struct ParamRelation
    hi::Num
    lo::Num
    sign::Int8
end

as_num(x::Num) = x
as_num(x) = Num(x)
ParamRelation(hi, lo, sign::Integer) =
    ParamRelation(as_num(hi), as_num(lo), Int8(sign))

function is_usable_rel(r::ParamRelation)
    return is_sym_member(r.hi) && is_sym_member(r.lo)
end

function is_sym_member(x::Num)
    u = SymbolicUtils.unwrap(x)
    u isa SymbolicUtils.BasicSymbolic || return false
    return SymbolicUtils.iscall(u) || SymbolicUtils.issym(u)
end

@inline function find_factor(m::Monomial, s)
    @inbounds for i in eachindex(m.syms)
        m.syms[i] === s && return i
    end
    return 0
end

function replace_exps(m::Monomial{E}, i::Int, ei::Rational{Int}, j::Int, ej::Rational{Int}) where {E}
    n = length(m.syms)
    syms = RawExpression[]
    exps = Rational{Int}[]
    sizehint!(syms, n)
    sizehint!(exps, n)
    @inbounds for k in 1:n
        e = k == i ? ei : (k == j ? ej : m.exps[k])
        iszero(e) && continue
        push!(syms, m.syms[k])
        push!(exps, e)
    end
    return Monomial{E}(monomial_scalar(m), syms, exps)
end

@inline reducible_power(m::Monomial, i::Int) =
    i != 0 && denominator(m.exps[i]) == 1 && numerator(m.exps[i]) >= 2

const MAX_RELATION_TERMS = 4096

function binomial_power(
        ::Type{E}, lo::SymbolicUtils.BasicSymbolic, sign::Int8, k::Int,
    ) where {E <: ExactScalar}
    out = Vector{Monomial{E}}(undef, k + 1)
    out[1] = constant_term(E, one(E))
    c = one(E)
    for j in 1:k
        c = c * E(k - j + 1, 0) * inv(E(j, 0))
        s = (sign < 0 && isodd(j)) ? -c : c
        out[j + 1] = Monomial{E}(s, RawExpression[lo], Rational{Int}[2j])
    end
    return out
end

function exceeds_relation_budget(terms::Vector{<:Monomial}, hi)
    projected = 0
    for m in terms
        i = find_factor(m, hi)
        add = if reducible_power(m, i)
            numerator(m.exps[i]) ÷ 2 + 1
        else
            1
        end
        add > MAX_RELATION_TERMS - projected && return true
        projected += add
    end
    return false
end

function reduce_pythagorean(terms::Vector{Monomial{E}}, hi, lo, sign::Int8) where {E}
    any(m -> reducible_power(m, find_factor(m, hi)), terms) || return (terms, false)
    exceeds_relation_budget(terms, hi) && return (terms, false)
    out = Monomial{E}[]
    pow = Monomial{E}[]
    last_k = -1
    for m in terms
        i = find_factor(m, hi)
        if !reducible_power(m, i)
            push!(out, m)
            continue
        end
        n = numerator(m.exps[i])
        k = n ÷ 2
        if k != last_k
            pow = binomial_power(E, lo, sign, k)
            last_k = k
        end
        base = replace_exps(m, i, Rational{Int}(n % 2), 0, 0 // 1)
        for p in pow
            push!(out, term_mul(base, p))
        end
    end
    return (canonical_terms!(out), true)
end

const MAX_REDUCE_PASSES = 8

function reduce_terms(::Type{E}, terms::Vector{Monomial{E}}, rels::Vector{ParamRelation}, gated::Bool) where {E}
    isempty(rels) && return terms
    cur = terms
    for _ in 1:MAX_REDUCE_PASSES
        progress = false
        for r in rels
            new, changed = reduce_pythagorean(
                cur, SymbolicUtils.unwrap(r.hi), SymbolicUtils.unwrap(r.lo), r.sign,
            )
            changed || continue
            (gated && length(new) >= length(cur)) && continue
            cur = new
            progress = true
        end
        progress || break
    end
    return cur
end

@inline function strip_conj(s)
    (s isa SymbolicUtils.BasicSymbolic && is_call_of(s, conj)) &&
        return (SymbolicUtils.arguments(s)[1], true)
    return (s, false)
end

@inline is_unary_call(x) =
    x isa SymbolicUtils.BasicSymbolic && SymbolicUtils.iscall(x) &&
    length(SymbolicUtils.arguments(x)) == 1

function trig_relations!(rels::Vector{ParamRelation}, terms::Vector{<:Monomial})
    from = length(rels) + 1
    trig_relations!(rels, (s for m in terms for s in m.syms))
    return sort_rels!(rels, from)
end

const TrigKey = Tuple{UInt, Bool, UInt8}

function trig_relations!(rels::Vector{ParamRelation}, factors)
    lows = Dict{TrigKey, RawExpression}()
    heads = Tuple{RawExpression, TrigKey, Int8}[]
    for s in factors
        core, wrapped = strip_conj(s)
        is_unary_call(core) || continue
        op = SymbolicUtils.operation(core)
        family = trig_family(op)
        iszero(family) && continue
        key = (objectid(only(SymbolicUtils.arguments(core))), wrapped, family)
        if op === sin || op === sinh
            lows[key] = s
        else
            push!(heads, (s, key, family == 1 ? Int8(-1) : Int8(1)))
        end
    end
    for (head, key, sign) in heads
        haskey(lows, key) || continue
        any(q -> SymbolicUtils.unwrap(q.hi) === head, rels) && continue
        push!(rels, ParamRelation(Num(head), Num(lows[key]), sign))
    end
    return rels
end

trig_family(op) = (op === cos || op === sin) ? UInt8(1) :
    ((op === cosh || op === sinh) ? UInt8(2) : UInt8(0))
is_trig(op) = !iszero(trig_family(op))

function sort_rels!(rels::Vector{ParamRelation}, from::Int)
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
    return rels
end

reduce_tail(t::Poly, rels::Vector{ParamRelation}, gated::Bool) =
    tiered(reduce_terms, t, rels, gated)

function reduce_trig(c::Coeff)
    t = c.tail
    t isa Poly || return c
    has_trig_factor(t.terms) || return c
    rels = trig_relations!(ParamRelation[], t.terms)
    isempty(rels) && return c
    return reduce_tail(t, rels, true)
end

function has_trig_factor(terms::Vector{<:Monomial})
    for m in terms, s in m.syms
        core = first(strip_conj(s))
        is_unary_call(core) || continue
        op = SymbolicUtils.operation(core)
        (op === cos || op === cosh) && return true
    end
    return false
end

function fresh_transient!(occupied::Vector{RawExpression})
    while true
        raw = SymbolicUtils.Sym{SymbolicUtils.SymReal}(gensym(:sqa_rel); type = Real)
        any(s -> isequal(s, raw), occupied) && continue
        push!(occupied, raw)
        return Num(raw)
    end
    return
end

function reduce_raw_via_transient(
        c::Coeff, rels::Vector{ParamRelation}, gated::Bool,
    )
    tail = c.tail
    tail isa RawSymbolicCoeff || return c
    raw = tail.expr
    fwd = Dict{Num, Num}()
    back = Dict{Num, Num}()
    trels = ParamRelation[]
    occupied = RawExpression[]
    append!(occupied, Symbolics.get_variables(Num(raw)))
    for r in rels
        haskey(fwd, r.hi) && continue
        th = fresh_transient!(occupied)
        tl = fresh_transient!(occupied)
        fwd[r.hi] = th
        fwd[r.lo] = tl
        back[th] = r.hi
        back[tl] = r.lo
        push!(trels, ParamRelation(th, tl, r.sign))
    end
    substituted = SymbolicUtils.unwrap(Symbolics.substitute(Num(raw), fwd))
    sub = to_cnum(substituted)
    sub.tail isa Poly || return c
    reduced = reduce_tail(sub.tail, trels, gated)
    isequal(reduced, sub) && return c
    reduced.tail isa Native && return reduced
    result = SymbolicUtils.unwrap(Symbolics.substitute(Num(raw_tail(reduced)), back))
    result isa RawExpression || return c
    return raw_coeff(result, tail.known_real)
end

is_trig_call(x::SymbolicUtils.BasicSymbolic) = SymbolicUtils.iscall(x) && is_trig(SymbolicUtils.operation(x))
has_symbolic_trig(x) = has_subterm(is_trig_call, x)

@inline is_trig_relation(r::ParamRelation) =
    has_symbolic_trig(SymbolicUtils.unwrap(r.hi)) ||
    has_symbolic_trig(SymbolicUtils.unwrap(r.lo))

function raw_uses_trig_relation(x, r::ParamRelation)
    for member in (r.hi, r.lo), trig in collect_trig!(RawExpression[], SymbolicUtils.unwrap(member))
        raw_depends_on(x, only(SymbolicUtils.arguments(trig))) && return true
    end
    return false
end

function collect_trig!(store::Vector{RawExpression}, x)
    x isa SymbolicUtils.BasicSymbolic || return store
    SymbolicUtils.iscall(x) || return store
    args = SymbolicUtils.arguments(x)
    if length(args) == 1 && is_trig(SymbolicUtils.operation(x))
        args[1] isa SymbolicUtils.BasicSymbolic && !any(y -> y === x, store) && push!(store, x)
        return store
    end
    for arg in args
        collect_trig!(store, arg)
    end
    return store
end

function reduce_all(
        c::Coeff, rels::Vector{ParamRelation}, gated::Bool,
        scratch::Vector{ParamRelation},
    )
    t = c.tail
    t isa Native && return c
    empty!(scratch)
    append!(scratch, rels)
    if t isa Poly
        has_trig_factor(t.terms) && trig_relations!(scratch, t.terms)
        isempty(scratch) && return c
        return reduce_tail(t, scratch, gated)
    end
    if has_symbolic_trig(t.expr)
        trig_relations!(scratch, collect_trig!(RawExpression[], t.expr))
    elseif all(is_trig_relation, rels)
        any(r -> raw_uses_trig_relation(t.expr, r), rels) || return c
    end
    isempty(scratch) && return c
    return reduce_raw_via_transient(c, scratch, gated)
end
