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

# Expanded rather than reached by `k` rounds of `_poly_mul`, which is quadratic: each round
# canonicalizes, so the accumulator grows by a term and is multiplied out again.
function _binomial_power(lo, sign::Int8, k::Int)
    out = Vector{Monomial}(undef, k + 1)
    out[1] = Monomial(_ONE_C, _EMPTY_SYMS, _EMPTY_EXPS)
    c = 1
    for j in 1:k
        c = c * (k - j + 1) ÷ j   # exact: the running binomial coefficient stays integral
        s = (sign < 0 && isodd(j)) ? -c : c
        out[j + 1] = Monomial(
            complex(Float64(s)), SymbolicUtils.BasicSymbolic[lo], Rational{Int}[2j],
        )
    end
    return out
end

# `hi^e -> (1 + sign*lo^2)^(e ÷ 2) * hi^(e % 2)`, term by term. Only integer exponents
# split; a radical of `hi` carries no such identity.
function _reduce_pythagorean(terms::Vector{Monomial}, hi, lo, sign::Int8)
    any(m -> _reducible_power(m, _find_factor(m, hi)), terms) || return (terms, false)
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
            pow = _binomial_power(lo, sign, k)
            last_k = k
        end
        base = Monomial[_replace_exps(m, i, Rational{Int}(n % 2), 0, 0 // 1)]
        append!(out, _poly_mul(base, pow))
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

# The stored `sin`/`sinh` factor over the same argument and with the same `conj` wrapping.
# Matching a stored factor rather than building one keeps the `===` factor identity the
# reduction relies on. Stays a scan: it returns on the first match, so discovery measures
# linear, and an index would cost more than it saves at one or two pairs.
function _find_partner(terms::Vector{Monomial}, op, arg, wrapped::Bool)
    for m in terms, s in m.syms
        core, w = _strip_conj(s)
        w == wrapped || continue
        a = _unary_arg(core)
        (a === nothing || !(a === arg)) && continue
        SymbolicUtils.operation(core) === op || continue
        return s
    end
    return nothing
end

# The Pythagorean partner of a head, with the sign of its relation. One table: the polynomial
# and symbolic discovery paths differ in how they match factors, not in which heads pair up.
_trig_partner(op) = op === cos ? (sin, Int8(-1)) :
    (op === cosh ? (sinh, Int8(1)) : nothing)

# `cos^2 + sin^2 = 1` and `cosh^2 - sinh^2 = 1` hold unconditionally, so the pairs come
# from the coefficient's own factors. Both members must carry the same `conj` wrapping:
# `_conj_atom` wraps `Number`-symtype atoms, and `conj(cos u)^2 + conj(sin u)^2 = 1` too.
_trig_relations(terms::Vector{Monomial}) = _trig_relations!(ParamRelation[], terms)

# Appends to `rels`, which may already hold the declared relations: a discovered pair whose
# `hi` is already there would only be applied twice for the same fixpoint.
function _trig_relations!(rels::Vector{ParamRelation}, terms::Vector{Monomial})
    from = length(rels) + 1
    for m in terms, s in m.syms
        core, wrapped = _strip_conj(s)
        arg = _unary_arg(core)
        arg === nothing && continue
        partner = _trig_partner(SymbolicUtils.operation(core))
        partner === nothing && continue
        lop, sign = partner
        any(q -> SymbolicUtils.unwrap(q.hi) === s, rels) && continue
        lo = _find_partner(terms, lop, arg, wrapped)
        lo === nothing && continue
        push!(rels, ParamRelation(Num(s), Num(lo), sign))
    end
    # `_fkey` is `objectid`, so factor order within a monomial is address-derived. Sort by
    # the structural hash instead, keeping the result independent of allocation order. Only
    # the pairs found here are sorted, so declared relations keep their leading position.
    _sort_rels!(rels, from)
    return rels
end

function _sort_rels!(rels::Vector{ParamRelation}, from::Int)
    length(rels) - from >= 1 && _insertion_sort!(
        view(rels, from:length(rels)),
        (a, b) -> hash(SymbolicUtils.unwrap(a.hi)) < hash(SymbolicUtils.unwrap(b.hi)),
    )
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
    t isa Poly || return c
    _has_trig_factor(t.terms) || return c   # `simplify`'s common path allocates nothing
    rels = _trig_relations(t.terms)
    isempty(rels) && return c
    return _reduce_tail(t, rels, true)
end

function _has_trig_factor(terms::Vector{Monomial})
    for m in terms, s in m.syms
        core, _ = _strip_conj(s)
        _unary_arg(core) === nothing && continue
        op = SymbolicUtils.operation(core)
        (op === cos || op === cosh) && return true
    end
    return false
end

# Stand-ins for the transient swap below; the name makes a clash with a user parameter a
# non-concern. Cached because `Symbolics.variable` costs ~1.6 µs a call even though it
# returns the same interned object. Filled in `__init__` rather than at load: a symbol minted
# during precompilation is absent from the runtime intern table, so a later rebuild would
# mint a second one under a different `_fkey`.
const _N_TRANSIENT = 8
const _TRANSIENT_SYMS = Num[]

_new_transient_sym(k::Int) = Symbolics.variable(Symbol("__sqa_rel_", k))

_transient_sym(k::Int) =
    k <= length(_TRANSIENT_SYMS) ? @inbounds(_TRANSIENT_SYMS[k]) : _new_transient_sym(k)

# A relation on a composite argument (`cos(ω*t)`) is not an atom, so its coefficient sits
# on the `Complex{Num}` tail, where the CAS folds degree 2 and nothing above it. Swap both
# members for plain symbols, which *are* atoms, reduce on the polynomial tier, swap back.
# The stand-ins never leave this function.
function _reduce_via_transient(c::Coeff, rels::Vector{ParamRelation}, gated::Bool)
    fwd = Dict{Num, Num}()
    back = Dict{Num, Num}()
    trels = ParamRelation[]
    for r in rels
        haskey(fwd, r.hi) && continue   # a declared pair the discovery found again
        k = 2 * length(trels)
        th = _transient_sym(k + 1)
        tl = _transient_sym(k + 2)
        fwd[r.hi] = th
        fwd[r.lo] = tl
        back[th] = r.hi
        back[tl] = r.lo
        push!(trels, ParamRelation(th, tl, r.sign))
    end
    re, im = _realimag(c)
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

# On the polynomial tier the factors of a coefficient are explicit; on the `Complex{Num}`
# tail they are buried in a SymbolicUtils tree, so walk it for the same pairs. Only the
# transform layer pays for this, never `simplify`.
function _collect_trig!(store::Vector{SymbolicUtils.BasicSymbolic}, x)
    x isa SymbolicUtils.BasicSymbolic || return nothing
    SymbolicUtils.iscall(x) || return nothing
    args = SymbolicUtils.arguments(x)
    op = SymbolicUtils.operation(x)
    if length(args) == 1 && (op === cos || op === sin || op === cosh || op === sinh)
        any(y -> y === x, store) || push!(store, x)
        return nothing
    end
    for arg in args
        _collect_trig!(store, arg)
    end
    return nothing
end

function _pair_trig!(
        rels::Vector{ParamRelation}, store::Vector{SymbolicUtils.BasicSymbolic},
        hop::Function, lop::Function, sign::Int8,
    )
    for hi in store
        SymbolicUtils.operation(hi) === hop || continue
        harg = SymbolicUtils.arguments(hi)[1]
        for lo in store
            SymbolicUtils.operation(lo) === lop || continue
            isequal(SymbolicUtils.arguments(lo)[1], harg) || continue
            push!(rels, ParamRelation(Num(hi), Num(lo), sign))
            break
        end
    end
    return nothing
end

function _sym_trig_relations!(rels::Vector{ParamRelation}, c::Coeff)
    store = SymbolicUtils.BasicSymbolic[]
    _collect_trig!(store, SymbolicUtils.unwrap(real(c)))
    _collect_trig!(store, SymbolicUtils.unwrap(imag(c)))
    isempty(store) && return rels
    for (hop, (lop, sign)) in ((cos, _trig_partner(cos)), (cosh, _trig_partner(cosh)))
        _pair_trig!(rels, store, hop, lop, sign)
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
        _trig_relations!(scratch, t.terms)
        isempty(scratch) && return c
        return _reduce_tail(t, scratch, gated)
    end
    _sym_trig_relations!(scratch, c)
    isempty(scratch) && return c
    return _reduce_via_transient(c, scratch, gated)
end
