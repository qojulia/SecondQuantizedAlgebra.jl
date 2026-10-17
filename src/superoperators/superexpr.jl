"""
    FormalSandwich

One sandwich `ρ ↦ left * ρ * right` of a [`SuperExpr`](@ref) with at least one formal side.
Both sides carry a unit coefficient; the scalar lives in the `SuperExpr`. A formal sandwich
carries neither a summation scope nor index constraints.
"""
struct FormalSandwich
    left::QExpr
    right::QExpr
    hash::UInt
    FormalSandwich(left::QExpr, right::QExpr) =
        new(left, right, hash(:FormalSandwich, hash(left, hash(right, zero(UInt)))))
end

Base.isequal(a::FormalSandwich, b::FormalSandwich) = a.hash == b.hash &&
    isequal(a.left, b.left) && isequal(a.right, b.right)
Base.:(==)(a::FormalSandwich, b::FormalSandwich) = isequal(a, b)
Base.hash(t::FormalSandwich, h::UInt) = hash(t.hash, h)

const FormalDict = Dict{FormalSandwich, CNum}

"""
    SuperExpr
    SuperExpr(L::Lindbladian)

A linear map on operators whose sides may contain formal operator functions.

A `SuperExpr` is a polynomial [`SuperOperator`](@ref) plus sandwiches `c * f ρ g` with a
[`QExpr`](@ref) side. It arises from [`sprepost`](@ref) with a formal side, from operations
on a `Lindbladian{QExpr}`, and from `SuperExpr(L)`. It supports the arithmetic,
composition, `transpose` and application of a `SuperOperator`, and keeps its type even when
no formal sandwich remains. `taylor(X, 0:n)`
lowers every formal side and returns a `SuperOperator`, and `SuperOperator(X)` converts an
`X` without formal sandwiches.

# Examples
```jldoctest
julia> using SecondQuantizedAlgebra: sprepost

julia> h = FockSpace(:c); a = Destroy(h, :a);

julia> X = sprepost(cos(a' * a), a)
cos(a' * a) ⋅ (·) ⋅ a

julia> X isa SecondQuantizedAlgebra.SuperExpr
true
```

See also [`SuperOperator`](@ref), [`taylor`](@ref).
"""
struct SuperExpr
    poly::SuperOperator
    formal::FormalDict
end

"""
    poly(X::SuperExpr)
    formal_terms(X::SuperExpr)

The polynomial part of a `SuperExpr`, and an iterator over its formal sandwiches with their
coefficients.

# Examples
```jldoctest
julia> using SecondQuantizedAlgebra: sprepost, spre, poly, formal_terms

julia> h = FockSpace(:c); a = Destroy(h, :a);

julia> X = sprepost(cos(a' * a), a) + spre(a);

julia> poly(X)
a ⋅ (·)

julia> length(collect(formal_terms(X)))
1
```

See also [`SuperExpr`](@ref).
"""
poly(X::SuperExpr) = X.poly
formal_terms(X::SuperExpr) = pairs(X.formal)
@doc (@doc poly) formal_terms

SuperExpr(S::SuperOperator) = SuperExpr(S, FormalDict())
SuperExpr(X::SuperExpr) = X
SuperExpr(L::Lindbladian{QAdd}) = SuperExpr(SuperOperator(L))

@noinline function formal_map_error()
    throw(
        ArgumentError(
            "the map has formal operator-function sides; lower it with `taylor(X, 0:n)` " *
                "before converting it to a `SuperOperator`",
        ),
    )
end

function SuperOperator(X::SuperExpr)
    isempty(X.formal) || formal_map_error()
    return X.poly
end

Base.zero(::Type{SuperExpr}) = SuperExpr(zero(SuperOperator))
Base.zero(::SuperExpr) = zero(SuperExpr)
Base.one(::Type{SuperExpr}) = SuperExpr(one(SuperOperator))
Base.one(::SuperExpr) = one(SuperExpr)
Base.iszero(X::SuperExpr) = iszero(X.poly) && isempty(X.formal)
Base.isone(X::SuperExpr) = isone(X.poly) && isempty(X.formal)

Base.isequal(a::SuperExpr, b::SuperExpr) = isequal(a.poly, b.poly) && isequal(a.formal, b.formal)
Base.:(==)(a::SuperExpr, b::SuperExpr) = isequal(a, b)
Base.hash(X::SuperExpr, h::UInt) = hash(:SuperExpr, hash(X.poly, hash(X.formal, h)))

const SideAtom = Tuple{CNum, QExprArg}

function side_atoms!(out::Vector{SideAtom}, c::CNum, x::QAdd)
    iszero(x) && return out
    isempty(x.indices) || formal_scope_error()
    for (t, v) in x.arguments
        push!(out, (c * v, side_qadd(t.ops, t.ne)))
    end
    return out
end

function side_atoms!(out::Vector{SideAtom}, c::CNum, x::QExpr)
    iszero(x) && return out
    if x.kind == QEXPR_ADD
        for arg in x.args
            side_atoms!(out, c, arg)
        end
        return out
    end
    qexpr_is_scalar(x) && return push!(out, (c * x.coeff, single_qadd(CNUM_ONE, Op[])))
    polynomial, body = polynomial_body(x)
    polynomial && return side_atoms!(out, c, body)
    return push!(out, (c * x.coeff, QExpr(x.kind, CNUM_ONE, copy(x.args))))
end

side_atoms(x::QExprArg) = side_atoms!(SideAtom[], CNUM_ONE, x)

atom_qexpr(x::QAdd) = qexpr_from_qadd(x)
atom_qexpr(x::QExpr) = x

struct SuperExprBuilder
    poly::SuperBuilder
    formal::FormalDict
end
SuperExprBuilder() = SuperExprBuilder(SuperBuilder(), FormalDict())

function add_atom_pair!(b::SuperExprBuilder, c::CNum, l::QExprArg, r::QExprArg)
    iszero(c) && return b
    if l isa QAdd && r isa QAdd
        tl, cl = first(l.arguments)
        tr, cr = first(r.arguments)
        canonical_sandwich!(b.poly.args, tl.ops, tr.ops, c * (cl * cr), merge_ne(tl.ne, tr.ne))
    else
        addto_key!(b.formal, FormalSandwich(atom_qexpr(l), atom_qexpr(r)), c)
    end
    return b
end

"""Add `c * left ρ right` with arbitrary polynomial or formal sides."""
function add_sandwich!(b::SuperExprBuilder, c::CNum, left::QExprArg, right::QExprArg)
    iszero(c) && return b
    for (cl, l) in side_atoms(left), (cr, r) in side_atoms(right)
        add_atom_pair!(b, c * (cl * cr), l, r)
    end
    return b
end

accumulate!(b::SuperExprBuilder, S::SuperOperator, scale::CNum = CNUM_ONE) =
    (accumulate!(b.poly, S, scale); b)

function accumulate!(b::SuperExprBuilder, X::SuperExpr, scale::CNum = CNUM_ONE)
    accumulate!(b.poly, X.poly, scale)
    iszero(scale) && return b
    for (t, c) in X.formal
        addto_key!(b.formal, t, c * scale)
    end
    return b
end

finish!(b::SuperExprBuilder) = SuperExpr(finish!(b.poly), b.formal)

sandwich_of(A::QAdd, B::QAdd) = sprepost_poly(A, B)
function sandwich_of(A::QExprArg, B::QExprArg)
    b = SuperExprBuilder()
    add_sandwich!(b, CNUM_ONE, A, B)
    return finish!(b)
end

function Base.:+(a::SuperExpr, b::SuperExpr)
    builder = SuperExprBuilder()
    accumulate!(builder, a)
    accumulate!(builder, b)
    return finish!(builder)
end
Base.:+(a::SuperExpr) = a

function scale(X::SuperExpr, c::CNum)
    builder = SuperExprBuilder()
    accumulate!(builder, X, c)
    return finish!(builder)
end
Base.:-(X::SuperExpr) = scale(X, CNUM_NEG1)
Base.:-(a::SuperExpr, b::SuperExpr) = a + (-b)
Base.:*(X::SuperExpr, c::Coefficient) = scale(X, to_cnum(c))
Base.:*(c::Coefficient, X::SuperExpr) = scale(X, to_cnum(c))

function lower_formal_hamiltonian!(b::SuperExprBuilder, H::QExpr)
    parts = H.kind == QEXPR_ADD ? H.args : QExprArg[H]
    for part in parts
        is_poly, p = part isa QAdd ? (true, part) : polynomial_qadd(part)
        if is_poly
            lower_hamiltonian!(b.poly, p)
        else
            one_side = single_qadd(CNUM_ONE, Op[])
            add_sandwich!(b, CNUM_NEG_IM, part, one_side)
            add_sandwich!(b, CNUM_IM, one_side, part)
        end
    end
    return b
end

function lower_formal_channel!(b::SuperExprBuilder, J::QExpr, K::QExpr, γ::CNum)
    one_side = single_qadd(CNUM_ONE, Op[])
    minus_half = CNUM_NEG1 * (CNUM_HALF * γ)
    for (cj, j) in side_atoms(J), (ck, k) in side_atoms(K)
        c = cj * conj(ck)
        kd = adjoint(k)
        add_sandwich!(b, γ * c, j, kd)
        kj = kd * j
        add_sandwich!(b, minus_half * c, kj, one_side)
        add_sandwich!(b, minus_half * c, one_side, kj)
    end
    return b
end

function SuperExpr(L::Lindbladian{QExpr})
    b = SuperExprBuilder()
    lower_formal_hamiltonian!(b, L.hamiltonian)
    for (ch, γ) in L.channels
        pj, J = polynomial_qadd(ch.jump)
        pk, K = polynomial_qadd(ch.partner)
        if pj && pk
            lower_channel!(b.poly, J, K, γ, ch.scope, ch.ne)
        else
            isempty(ch.scope) || formal_scope_error()
            lower_formal_channel!(b, ch.jump, ch.partner, γ)
        end
    end
    return finish!(b)
end

const AnyMap = Union{SuperOperator, SuperExpr, Lindbladian}

normal_to_symmetric(::AnyMap) = weyl_map_error()
symmetric_to_normal(::AnyMap) = weyl_map_error()

const PolyMap = Union{SuperOperator, Lindbladian{QAdd}}

lift(a::PolyMap, b::PolyMap) = (SuperOperator(a), SuperOperator(b))
lift(a::AnyMap, b::AnyMap) = (SuperExpr(a), SuperExpr(b))

lowered(S::Union{SuperOperator, SuperExpr}) = S
lowered(L::Lindbladian{QAdd}) = SuperOperator(L)
lowered(L::Lindbladian{QExpr}) = SuperExpr(L)

Base.:/(X::AnyMap, b::Number) =
    b isa Integer && !(b isa Bool) && !iszero(b) ? X * (1 // b) : X * inv(b)
Base.:/(X::AnyMap, b::Coefficient) = X * inv(b)
Base.://(X::AnyMap, b::Integer) = X * (1 // b)
Base.://(X::AnyMap, b::Coefficient) = X / b

function Base.:+(X::AnyMap, c::Coefficient)
    Y = lowered(X)
    return Y + scale(one(Y), to_cnum(c))
end
Base.:+(c::Coefficient, X::AnyMap) = X + c
Base.:-(X::AnyMap, c::Coefficient) = X + (-to_cnum(c))
Base.:-(c::Coefficient, X::AnyMap) = (-lowered(X)) + c

Base.transpose(L::Lindbladian) = transpose(lowered(L))
Base.adjoint(L::Lindbladian) = adjoint(lowered(L))
dagger(X::AnyMap) = adjoint(X)
qadjoint(X::AnyMap) = adjoint(X)
(L::Lindbladian)(ρ::Union{QField, Coefficient}) = lowered(L)(ρ)

Base.inv(::AnyMap) = negative_power_error()

function Base.:^(X::AnyMap, n::Integer)
    n >= 0 || negative_power_error()
    Y = lowered(X)
    n == 0 && return one(Y)
    n >= 2 && !isempty(get_sum_indices(Y)) && scoped_power_error()
    result = Y
    for _ in 2:n
        result = compose(result, Y)
    end
    return result
end

Base.:+(a::AnyMap, b::AnyMap) = +(lift(a, b)...)
Base.:-(a::AnyMap, b::AnyMap) = -(lift(a, b)...)
compose(a::AnyMap, b::AnyMap) = compose(lift(a, b)...)
Base.:*(a::AnyMap, b::AnyMap) = compose(a, b)
Base.:∘(a::AnyMap, b::AnyMap) = compose(a, b)
commutator(a::AnyMap, b::AnyMap) = compose(a, b) - compose(b, a)
anticommutator(a::AnyMap, b::AnyMap) = compose(a, b) + compose(b, a)


function poly_sandwiches(S::SuperOperator)
    out = Tuple{CNum, QExpr, QExpr}[]
    isempty(S.arguments) && return out
    isempty(S.indices) || formal_scope_error()
    for (t, c) in S.arguments
        push!(out, (c, qexpr_from_qadd(side_qadd(t.left, t.ne)), qexpr_from_qadd(side_qadd(t.right, t.ne))))
    end
    return out
end
formal_sandwiches(X::SuperExpr) = Tuple{CNum, QExpr, QExpr}[(c, t.left, t.right) for (t, c) in X.formal]

function compose_sandwiches!(
        b::SuperExprBuilder, as::Vector{Tuple{CNum, QExpr, QExpr}},
        bs::Vector{Tuple{CNum, QExpr, QExpr}},
    )
    for (ca, al, ar) in as, (cb, bl, br) in bs
        add_sandwich!(b, ca * cb, al * bl, br * ar)
    end
    return b
end

function compose(A::SuperExpr, B::SuperExpr)
    b = SuperExprBuilder()
    accumulate!(b, compose(A.poly, B.poly))
    af = formal_sandwiches(A)
    bf = formal_sandwiches(B)
    compose_sandwiches!(b, af, bf)
    isempty(af) || compose_sandwiches!(b, af, poly_sandwiches(B.poly))
    isempty(bf) || compose_sandwiches!(b, poly_sandwiches(A.poly), bf)
    return finish!(b)
end

function apply_formal(S::SuperOperator, ρ::QExpr)
    iszero(S) && return qexpr_zero()
    isempty(S.indices) || formal_scope_error()
    terms = QExprArg[]
    for (t, c) in S.arguments
        push!(terms, qexpr_product(QExprArg[side_qadd(t.left, t.ne), ρ, side_qadd(t.right, t.ne)], c))
    end
    return qexpr_sum(terms)
end

(S::SuperOperator)(ρ::QExpr) = apply_formal(S, ρ)

function apply_map(X::SuperExpr, ρ::QExprArg)
    terms = QExprArg[]
    iszero(X.poly) || push!(terms, ρ isa QAdd ? X.poly(ρ) : apply_formal(X.poly, ρ))
    for (t, c) in X.formal
        push!(terms, qexpr_product(QExprArg[t.left, ρ, t.right], c))
    end
    return qexpr_sum(terms)
end

(X::SuperExpr)(ρ::Union{QAdd, QExpr}) = apply_map(X, ρ)
(X::SuperExpr)(ρ::QSym) = apply_map(X, to_qadd(ρ))
(X::SuperExpr)(ρ::Coefficient) = apply_map(X, single_qadd(to_cnum(ρ), Op[]))

function Base.adjoint(X::SuperExpr)
    b = SuperExprBuilder()
    accumulate!(b, adjoint(X.poly))
    for (t, c) in X.formal
        add_sandwich!(b, conj(c), adjoint(t.left), adjoint(t.right))
    end
    return finish!(b)
end

function Base.transpose(X::SuperExpr)
    b = SuperExprBuilder()
    accumulate!(b, transpose(X.poly))
    for (t, c) in X.formal
        addto_key!(b.formal, FormalSandwich(t.right, t.left), c)
    end
    return finish!(b)
end

"""
    rewrite(f, fc, X)

`X` with `f` applied to every operator side and `fc` to every coefficient. A `Lindbladian`
renormalizes its channels; maps that change a channel's scope use `map_lindbladian`.
"""
function rewrite(f::F, fc::FC, X::SuperExpr) where {F, FC}
    b = SuperExprBuilder()
    accumulate!(b, f(X.poly))
    for (t, c) in X.formal
        add_sandwich!(b, fc(c), f(t.left), f(t.right))
    end
    return finish!(b)
end

rewrite(f::F, fc::FC, L::Lindbladian{O}) where {F, FC, O} =
    map_lindbladian(f, f, fc, same_scope, O, L)

const StructuredMap = Union{SuperExpr, Lindbladian}

normal_order(X::StructuredMap) = rewrite(normal_order, identity, X)
SymbolicUtils.simplify(X::StructuredMap; kwargs...) = rewrite(
    x -> SymbolicUtils.simplify(x; kwargs...), c -> simplify_prefactor(c; kwargs...), X,
)
Symbolics.expand(X::StructuredMap; kwargs...) =
    rewrite(x -> Symbolics.expand(x; kwargs...), c -> expand_prefactor(c; kwargs...), X)
expand_completeness(X::StructuredMap) = rewrite(expand_completeness, identity, X)
conjugate(X::StructuredMap, U::UnitaryTransform) = rewrite(x -> conjugate(x, U), identity, X)

function SymbolicUtils.substitute(X::StructuredMap, rules::AbstractDict; replace_adjoint::Bool = true)
    r = SubstitutionRules(rules, replace_adjoint)
    return rewrite(x -> substitute_with(x, r), c -> substitute_cnum(c, r), X)
end

assume_distinct_index(X::SuperExpr, pairs::Vector{Tuple{Index, Index}}) =
    rewrite(x -> assume_distinct_index(x, pairs), identity, X)
change_index(X::SuperExpr, m::IndexMap) =
    rewrite(x -> change_index(x, m), c -> change_index(c, m), X)

transform(X::SuperExpr, U::UnitaryTransform{StaticTime}) = conjugate(X, U)
transform(X::SuperExpr, U::UnitaryTransform{DynamicTime}) =
    conjugate(X, U) + SuperExpr(SuperOperator(hamiltonian_action(U.gauge)))

function get_indices(X::SuperExpr)
    out = get_indices(X.poly)
    for t in keys(X.formal)
        merge_indices!(out, get_indices(t.left))
        merge_indices!(out, get_indices(t.right))
    end
    return out
end

function foreach_operand(f::F, X::SuperExpr) where {F}
    f(X.poly)
    for t in keys(X.formal)
        f(t.left)
        f(t.right)
    end
    return nothing
end

function foreach_operand(f::F, L::Lindbladian) where {F}
    f(L.hamiltonian)
    for ch in keys(L.channels)
        f(ch.jump)
        f(ch.partner)
    end
    return nothing
end

rates(X::SuperExpr) = values(X.formal)
rates(L::Lindbladian) = values(L.channels)

function get_operators(X::StructuredMap)
    ops = Op[]
    foreach_operand(x -> append!(ops, get_operators(x)), X)
    unique!(ops)
    return sort!(ops; by = order_key)
end

function acts_on(X::StructuredMap)
    out = Int[]
    foreach_operand(x -> append!(out, acts_on(x)), X)
    return unique!(sort!(out))
end

function Symbolics.get_variables(X::StructuredMap; kwargs...)
    vars = SymbolicUtils.BasicSymbolic[]
    foreach_operand(x -> append!(vars, Symbolics.get_variables(x; kwargs...)), X)
    for c in rates(X)
        append!(vars, Symbolics.get_variables(real(c); kwargs...))
        append!(vars, Symbolics.get_variables(imag(c); kwargs...))
    end
    unique!(vars)
    sort!(vars; by = string)
    return vars
end

get_sum_indices(X::SuperExpr) = get_sum_indices(X.poly)

"""
    taylor(X::SuperExpr, ns)
    taylor(L::Lindbladian{QExpr}, ns)

Lower every formal side of a map to its Maclaurin polynomial through order `n`, for
`ns = 0:n`.

A `SuperExpr` becomes a `SuperOperator` and a `Lindbladian{QExpr}` a polynomial
`Lindbladian`. Each formal function is expanded as for an operator.
"""
function Symbolics.taylor(X::SuperExpr, ns::UnitRange{T}) where {T <: Integer}
    order = taylor_prefix_order(ns)
    b = SuperBuilder()
    accumulate!(b, X.poly)
    for (t, c) in X.formal
        accumulate!(b, c, lower_taylor_qexpr(t.left, order), lower_taylor_qexpr(t.right, order))
    end
    return finish!(b)
end

function Symbolics.taylor(L::Lindbladian{QExpr}, ns::UnitRange{T}) where {T <: Integer}
    order = taylor_prefix_order(ns)
    lowered_side(q) = lower_taylor_qexpr(q, order)
    return map_lindbladian(lowered_side, lowered_side, identity, same_scope, QAdd, L)
end
