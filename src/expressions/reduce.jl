# The Pythagorean relation `hi^2 == 1 + sign * lo^2`, used to reduce coefficients without
# calling the CAS. `hi` is eliminated in favour of `lo`, so the factor to remove comes first.
# `sign` is `-1` for a trigonometric pair (`cos^2 = 1 - sin^2`) and `+1` for a hyperbolic one.
struct ParamRelation
    hi::Num
    lo::Num
    sign::Int8
end

_as_num(x::Num) = x
_as_num(x) = Num(x)
ParamRelation(hi, lo, sign::Integer) =
    ParamRelation(_as_num(hi), _as_num(lo), Int8(sign))

# Whether a relation can fire at all: both members have to be symbolic, not folded
# constants. A numeric parameter collapses its members to literals, and the transient swap
# would then hand those to `substitute` and rewrite bare constants in every coefficient.
function _is_usable_rel(r::ParamRelation)
    return _is_sym_member(r.hi) && _is_sym_member(r.lo)
end

function _is_sym_member(x::Num)
    u = SymbolicUtils.unwrap(x)
    u isa SymbolicUtils.BasicSymbolic || return false
    return SymbolicUtils.iscall(u) || SymbolicUtils.issym(u)
end

@inline function _find_factor(m::Monomial, s)
    @inbounds for i in eachindex(m.syms)
        m.syms[i] === s && return i
    end
    return 0
end

# Rebuild with the exponents at `i` and `j` replaced, dropping factors that cancel.
# Index order is preserved, so `syms` stays sorted by `_fkey` as `_term_mul` and
# `_poly_add` require.
function _replace_exps(m::Monomial, i::Int, ei::Rational{Int}, j::Int, ej::Rational{Int})
    n = length(m.syms)
    syms = SymbolicUtils.BasicSymbolic[]
    exps = Rational{Int}[]
    sizehint!(syms, n)
    sizehint!(exps, n)
    @inbounds for k in 1:n
        e = k == i ? ei : (k == j ? ej : m.exps[k])
        iszero(e) && continue
        push!(syms, m.syms[k])
        push!(exps, e)
    end
    return Monomial(m.scalar, syms, exps)
end

@inline _reducible_power(m::Monomial, i::Int) =
    i != 0 && denominator(m.exps[i]) == 1 && numerator(m.exps[i]) >= 2

# A malformed or simply enormous symbolic power must not turn simplification into an
# accidental expansion bomb. This is a work bound, not a mathematical approximation: on
# refusal the original coefficient is returned byte-for-byte unchanged.
const _MAX_RELATION_TERMS = 4096

# Expanded rather than reached by `k` rounds of `_poly_mul`, which is quadratic: each round
# canonicalizes, so the accumulator grows by a term and is multiplied out again. `nothing`
# means that an integer binomial coefficient overflowed and the caller must abandon the
# whole rewrite.
function _binomial_power(lo, sign::Int8, k::Int)
    out = Vector{Monomial}(undef, k + 1)
    out[1] = Monomial(_ONE_C, _EMPTY_SYMS, _EMPTY_EXPS)
    c = 1
    for j in 1:k
        factor = k - j + 1
        c > typemax(Int) ÷ factor && return nothing
        c = c * factor ÷ j   # exact: the running binomial coefficient stays integral
        s = (sign < 0 && isodd(j)) ? -c : c
        out[j + 1] = Monomial(
            complex(Float64(s)), SymbolicUtils.BasicSymbolic[lo], Rational{Int}[2j],
        )
    end
    return out
end

# Project the raw number of terms before allocating a binomial. Canonicalization can make
# this much shorter (and often reduces it to one), so the projection is only a safety gate;
# the `gated` length decision is made on the canonical result afterwards.
function _projected_relation_terms(terms::Vector{Monomial}, hi)
    projected = 0
    for m in terms
        i = _find_factor(m, hi)
        add = if _reducible_power(m, i)
            numerator(m.exps[i]) ÷ 2 + 1
        else
            1
        end
        add > _MAX_RELATION_TERMS - projected && return nothing
        projected += add
    end
    return projected
end

# `hi^e -> (1 + sign*lo^2)^(e ÷ 2) * hi^(e % 2)`, term by term. Only integer exponents
# split; a radical of `hi` carries no such identity.
function _reduce_pythagorean(terms::Vector{Monomial}, hi, lo, sign::Int8)
    reducible = false
    for m in terms
        if _reducible_power(m, _find_factor(m, hi))
            reducible = true
            break
        end
    end
    reducible || return (terms, false)
    _projected_relation_terms(terms, hi) === nothing && return (terms, false)
    out = Monomial[]
    pow = Monomial[]
    last_k = -1   # terms of one expression almost always share an exponent
    for m in terms
        i = _find_factor(m, hi)
        if !_reducible_power(m, i)
            push!(out, m)
            continue
        end
        n = numerator(m.exps[i])
        k = n ÷ 2
        if k != last_k
            newpow = _binomial_power(lo, sign, k)
            newpow === nothing && return (terms, false)
            pow = newpow
            last_k = k
        end
        base = _replace_exps(m, i, Rational{Int}(n % 2), 0, 0 // 1)
        for p in pow
            push!(out, _term_mul(base, p))
        end
    end
    return (_canonical_terms!(out), true)
end

_apply_relation(terms::Vector{Monomial}, r::ParamRelation) = _reduce_pythagorean(
    terms, SymbolicUtils.unwrap(r.hi), SymbolicUtils.unwrap(r.lo), r.sign,
)

# Two relations that eliminate each other's retained factor would ping-pong forever;
# the cap is the backstop. Every relation used here converges in one pass.
const _MAX_REDUCE_PASSES = 8

# Apply every relation until nothing changes.
#
# `gated` keeps a rewrite only when it shortens the term list, which is what display
# wants but makes reduction non-monotone: the sole route to zero may pass through a
# longer intermediate. Zero-testing therefore runs ungated.
function _reduce_terms(terms::Vector{Monomial}, rels::Vector{ParamRelation}, gated::Bool)
    isempty(rels) && return terms
    cur = terms
    for _ in 1:_MAX_REDUCE_PASSES
        progress = false
        for r in rels
            new, changed = _apply_relation(cur, r)
            changed || continue
            (gated && length(new) >= length(cur)) && continue
            cur = new
            progress = true
        end
        progress || break
    end
    return cur
end

@inline function _strip_conj(s)
    (s isa SymbolicUtils.BasicSymbolic && _is_conj_call(s)) &&
        return (SymbolicUtils.arguments(s)[1], true)
    return (s, false)
end

# The sole argument of a one-argument call, else `nothing`.
@inline function _unary_arg(x)
    (x isa SymbolicUtils.BasicSymbolic && SymbolicUtils.iscall(x)) || return nothing
    args = SymbolicUtils.arguments(x)
    length(args) == 1 || return nothing
    return args[1]
end

struct _TrigKey
    arg_id::UInt
    wrapped::Bool
    family::UInt8
end

Base.isequal(a::_TrigKey, b::_TrigKey) =
    a.arg_id == b.arg_id && a.wrapped == b.wrapped && a.family == b.family
Base.hash(k::_TrigKey, h::UInt) = hash(k.family, hash(k.wrapped, hash(k.arg_id, h)))

struct _TrigHead{F}
    factor::F
    key::_TrigKey
    sign::Int8
end

# `cos^2 + sin^2 = 1` and `cosh^2 - sinh^2 = 1` hold unconditionally, so the pairs come
# from the coefficient's own factors. Both members must carry the same `conj` wrapping:
# `_conj_atom` wraps `Number`-symtype atoms, and `conj(cos u)^2 + conj(sin u)^2 = 1` too.
_trig_relations(terms::Vector{Monomial}) = _trig_relations!(ParamRelation[], terms)

# Appends to `rels`, which may already hold the declared relations: a discovered pair whose
# `hi` is already there would only be applied twice for the same fixpoint.
function _trig_relations!(rels::Vector{ParamRelation}, terms::Vector{Monomial})
    from = length(rels) + 1
    lows = Dict{_TrigKey, SymbolicUtils.BasicSymbolic}()
    heads = _TrigHead[]
    for m in terms, s in m.syms
        core, wrapped = _strip_conj(s)
        arg = _unary_arg(core)
        arg isa SymbolicUtils.BasicSymbolic || continue
        op = SymbolicUtils.operation(core)
        family = (op === cos || op === sin) ? UInt8(1) :
            ((op === cosh || op === sinh) ? UInt8(2) : UInt8(0))
        iszero(family) && continue
        key = _TrigKey(objectid(arg), wrapped, family)
        if op === sin || op === sinh
            lows[key] = s
        else
            push!(heads, _TrigHead(s, key, family == 1 ? Int8(-1) : Int8(1)))
        end
    end
    for head in heads
        duplicate = false
        for q in rels
            if SymbolicUtils.unwrap(q.hi) === head.factor
                duplicate = true
                break
            end
        end
        duplicate && continue
        lo = get(lows, head.key, nothing)
        lo === nothing && continue
        push!(rels, ParamRelation(Num(head.factor), Num(lo), head.sign))
    end
    # `_fkey` is `objectid`, so factor order within a monomial is address-derived. Sort by
    # the structural hash instead, keeping the result independent of allocation order. Only
    # the pairs found here are sorted, so declared relations keep their leading position.
    _sort_rels!(rels, from)
    return rels
end

function _sort_rels!(rels::Vector{ParamRelation}, from::Int)
    @inbounds for i in (from + 1):length(rels)
        r = rels[i]
        key = hash(SymbolicUtils.unwrap(r.hi))
        j = i - 1
        while j >= from && hash(SymbolicUtils.unwrap(rels[j].hi)) > key
            rels[j + 1] = rels[j]
            j -= 1
        end
        rels[j + 1] = r
    end
    return nothing
end

function _reduce_tail(t::Poly, rels::Vector{ParamRelation}, gated::Bool)
    out = _reduce_terms(t.terms, rels, gated)
    out === t.terms && return _poly_coeff(t)
    return _from_poly(out)   # empty -> `_CNUM_ZERO`; a `Poly` never reports itself zero
end

# Fold `cos^2 + sin^2` and `cosh^2 - sinh^2` factors of a coefficient.
function _reduce_trig(c::Coeff)
    t = c.tail
    t isa Native && return c
    if t isa Poly
        _has_trig_factor(t.terms) || return c   # common path allocates nothing
        rels = _trig_relations(t.terms)
        isempty(rels) && return c
        return _reduce_tail(t, rels, true)
    end
    _has_symbolic_trig(t.expr) || return c
    rels = ParamRelation[]
    _sym_trig_relations!(rels, c)
    isempty(rels) && return c
    return _reduce_via_transient(c, rels, true)
end

function _has_trig_factor(terms::Vector{Monomial})
    @inbounds for mi in eachindex(terms)
        syms = terms[mi].syms
        for si in eachindex(syms)
            s = syms[si]
            SymbolicUtils.iscall(s) || continue
            core = s
            if SymbolicUtils.operation(core) === conj
                args = SymbolicUtils.arguments(core)
                length(args) == 1 || continue
                inner = args[1]
                inner isa SymbolicUtils.BasicSymbolic || continue
                SymbolicUtils.iscall(inner) || continue
                core = inner
            end
            op = SymbolicUtils.operation(core)
            (op === cos || op === cosh) && return true
        end
    end
    return false
end

# Create a genuinely fresh stand-in and check it against every symbol already present in
# the coefficient. The explicit check is intentional: these atoms participate in global
# SymbolicUtils interning, so relying on a human-chosen prefix is never sound.
function _fresh_transient!(occupied::Vector{SymbolicUtils.BasicSymbolic})
    while true
        candidate = Symbolics.variable(gensym(:sqa_rel))
        raw = SymbolicUtils.unwrap(candidate)
        collision = false
        for s in occupied
            if isequal(s, raw)
                collision = true
                break
            end
        end
        collision && continue
        push!(occupied, raw)
        return candidate
    end
    return
end

# A relation on a composite argument (`cos(ω*t)`) is not an atom, so its coefficient sits
# on the raw symbolic tail, where the CAS folds degree 2 and nothing above it. Swap both
# members for plain symbols, which *are* atoms, reduce on the polynomial tier, swap back.
# The raw route keeps the complex expression intact; splitting it into real and imaginary
# parts here would materialize the representation on every coefficient.
function _reduce_raw_via_transient(
        c::Coeff, rels::Vector{ParamRelation}, gated::Bool,
    )
    tail = c.tail::RawSymbolicCoeff
    raw = tail.expr
    fwd = Dict{Num, Num}()
    back = Dict{Num, Num}()
    trels = ParamRelation[]
    occupied = SymbolicUtils.BasicSymbolic[]
    append!(occupied, Symbolics.get_variables(Num(raw)))
    for r in rels
        haskey(fwd, r.hi) && continue
        th = _fresh_transient!(occupied)
        tl = _fresh_transient!(occupied)
        fwd[r.hi] = th
        fwd[r.lo] = tl
        back[th] = r.hi
        back[tl] = r.lo
        push!(trels, ParamRelation(th, tl, r.sign))
    end
    substituted = SymbolicUtils.unwrap(Symbolics.substitute(Num(raw), fwd))
    sub = _recognize(substituted)
    sub.tail isa Poly || return c
    reduced = _reduce_tail(sub.tail, trels, gated)
    isequal(reduced, sub) && return c
    reduced.tail isa Native && return reduced
    result = SymbolicUtils.unwrap(Symbolics.substitute(Num(_raw_expression(reduced)), back))
    value = _const_value(result)
    value isa Number && return _to_cnum(value)
    result isa SymbolicUtils.BasicSymbolic || return c
    return _from_raw_arithmetic(result, tail.real_slot)
end

function _reduce_via_transient(c::Coeff, rels::Vector{ParamRelation}, gated::Bool)
    c.tail isa RawSymbolicCoeff && return _reduce_raw_via_transient(c, rels, gated)
    re, im = _raw_parts(c)
    return _reduce_via_transient(c, rels, gated, re, im)
end

function _reduce_via_transient(
        c::Coeff, rels::Vector{ParamRelation}, gated::Bool, re::Num, im::Num,
    )
    fwd = Dict{Num, Num}()
    back = Dict{Num, Num}()
    trels = ParamRelation[]
    occupied = SymbolicUtils.BasicSymbolic[]
    append!(occupied, Symbolics.get_variables(re))
    append!(occupied, Symbolics.get_variables(im))
    for r in rels
        haskey(fwd, r.hi) && continue   # a declared pair the discovery found again
        th = _fresh_transient!(occupied)
        tl = _fresh_transient!(occupied)
        fwd[r.hi] = th
        fwd[r.lo] = tl
        back[th] = r.hi
        back[tl] = r.lo
        push!(trels, ParamRelation(th, tl, r.sign))
    end
    sub = _cnum(Symbolics.substitute(re, fwd), Symbolics.substitute(im, fwd))
    st = sub.tail
    st isa Poly || return c
    red = _reduce_tail(st, trels, gated)
    isequal(red, sub) && return c
    rre, rim = _realimag(red)
    return _cnum(Symbolics.substitute(rre, back), Symbolics.substitute(rim, back))
end

# Reduce a coefficient by the given parameter relations.
function _reduce_cnum(c::Coeff, rels::Vector{ParamRelation}, gated::Bool)
    isempty(rels) && return c
    t = c.tail
    t isa Poly && return _reduce_tail(t, rels, gated)
    t isa Native && return c
    return _reduce_via_transient(c, rels, gated)
end

# On the polynomial tier the factors of a coefficient are explicit; on the raw symbolic
# tail they are buried in a SymbolicUtils tree, so walk it for the same pairs. Only the
# transform layer pays for this, never `simplify`.
function _has_symbolic_trig(x)
    x isa SymbolicUtils.BasicSymbolic || return false
    SymbolicUtils.iscall(x) || return false
    op = SymbolicUtils.operation(x)
    (op === cos || op === sin || op === cosh || op === sinh) && return true
    for arg in SymbolicUtils.arguments(x)
        _has_symbolic_trig(arg) && return true
    end
    return false
end

@inline _is_trig_relation(r::ParamRelation) =
    _has_symbolic_trig(SymbolicUtils.unwrap(r.hi)) ||
    _has_symbolic_trig(SymbolicUtils.unwrap(r.lo))

function _raw_uses_trig_relation(x, r::ParamRelation)
    for member in (r.hi, r.lo)
        raw = SymbolicUtils.unwrap(member)
        _has_symbolic_trig(raw) || continue
        args = SymbolicUtils.arguments(raw)
        isempty(args) || _raw_depends_on(x, only(args)) && return true
    end
    return false
end

function _collect_trig!(store::Vector{SymbolicUtils.BasicSymbolic}, x)
    x isa SymbolicUtils.BasicSymbolic || return nothing
    SymbolicUtils.iscall(x) || return nothing
    args = SymbolicUtils.arguments(x)
    op = SymbolicUtils.operation(x)
    if length(args) == 1 && (op === cos || op === sin || op === cosh || op === sinh)
        duplicate = false
        for y in store
            if y === x
                duplicate = true
                break
            end
        end
        duplicate || push!(store, x)
        return nothing
    end
    for arg in args
        _collect_trig!(store, arg)
    end
    return nothing
end

function _sym_trig_relations!(rels::Vector{ParamRelation}, c::Coeff)
    if c.tail isa RawSymbolicCoeff
        store = SymbolicUtils.BasicSymbolic[]
        _collect_trig!(store, c.tail.expr)
        return _sym_trig_relations!(rels, store)
    end
    re, im = _raw_parts(c)
    return _sym_trig_relations!(rels, c, re, im)
end

function _sym_trig_relations!(rels::Vector{ParamRelation}, c::Coeff, re::Num, im::Num)
    store = SymbolicUtils.BasicSymbolic[]
    _collect_trig!(store, SymbolicUtils.unwrap(re))
    _collect_trig!(store, SymbolicUtils.unwrap(im))
    return _sym_trig_relations!(rels, store)
end

function _sym_trig_relations!(rels::Vector{ParamRelation}, store::Vector{SymbolicUtils.BasicSymbolic})
    isempty(store) && return rels
    lows = Dict{_TrigKey, SymbolicUtils.BasicSymbolic}()
    heads = _TrigHead[]
    for factor in store
        op = SymbolicUtils.operation(factor)
        family = (op === cos || op === sin) ? UInt8(1) : UInt8(2)
        arg = only(SymbolicUtils.arguments(factor))
        arg isa SymbolicUtils.BasicSymbolic || continue
        key = _TrigKey(objectid(arg), false, family)
        if op === sin || op === sinh
            lows[key] = factor
        else
            push!(heads, _TrigHead(factor, key, family == 1 ? Int8(-1) : Int8(1)))
        end
    end
    for head in heads
        lo = get(lows, head.key, nothing)
        lo === nothing && continue
        push!(rels, ParamRelation(Num(head.factor), Num(lo), head.sign))
    end
    return rels
end

# Reduce by the declared relations together with the trigonometric pairs found in `c`, to a
# fixpoint over all of them.
function _reduce_all(
        c::Coeff, rels::Vector{ParamRelation}, gated::Bool,
        scratch::Vector{ParamRelation},
    )
    t = c.tail
    t isa Native && return c
    # Declared relations lead, then the pairs found in `c`. `scratch` is reused across the
    # terms of one expression, so a frame conjugation allocates no relation list per term.
    empty!(scratch)
    append!(scratch, rels)
    if t isa Poly
        _has_trig_factor(t.terms) && _trig_relations!(scratch, t.terms)
        isempty(scratch) && return c
        return _reduce_tail(t, scratch, gated)
    end
    if _has_symbolic_trig(t.expr)
        _sym_trig_relations!(scratch, c)
        isempty(scratch) && return c
        return _reduce_via_transient(c, scratch, gated)
    end
    # A raw coefficient without trig or phase factors cannot be affected by the
    # trigonometric relations attached to a unitary transform. In particular,
    # leave indexed parameters and time derivatives in their raw representation.
    if all(_is_trig_relation, rels)
        any(r -> _raw_uses_trig_relation(t.expr, r), rels) || return c
    end
    isempty(scratch) && return c
    return _reduce_via_transient(c, scratch, gated)
end
