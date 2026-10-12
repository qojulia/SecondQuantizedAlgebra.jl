# Cold-path superoperators with formal (QExpr) sides.

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

A superoperator with formal operator-function sides, `poly + Σ c * f ρ g`. The polynomial
part is a [`SuperOperator`](@ref); every sandwich with a [`QExpr`](@ref) side is a
`FormalSandwich` with its scalar moved into the coefficient. Sides are split into
their additive atoms, so `c * (f, g)`, `(c f, g)` and `(f, c g)` are one key, and a sandwich
whose sides are both polynomial lands in the polynomial part.

Operations keep the `SuperExpr` type even when the formal part is empty; call
`SuperOperator(X)` to convert, or `taylor(X, 0:n)` to lower every formal side.
"""
struct SuperExpr
    poly::SuperOperator
    formal::FormalDict
end

"""
    poly(X::SuperExpr) -> SuperOperator
    formal_terms(X::SuperExpr)

The polynomial part of `X`, and an iterator of its `FormalSandwich => coefficient` pairs.
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

## Atom decomposition

const SideAtom = Tuple{CNum, QExprArg}

function side_atoms!(out::Vector{SideAtom}, c::CNum, x::QAdd)
    iszero(x) && return out
    isempty(x.indices) || formal_scope_error()
    for (t, v) in x.arguments
        push!(out, (mul_cnum(c, v), QAdd(QTermDict(copy_key(t) => CNUM_ONE), Index[])))
    end
    return out
end

function side_atoms!(out::Vector{SideAtom}, c::CNum, x::QExpr)
    iszero(x) && return out
    if x.kind == QEXPR_ADD
        for arg in x.args
            side_atoms!(out, c, arg)
        end
    elseif qexpr_is_scalar(x)
        push!(out, (mul_cnum(c, x.coeff), single_qadd(CNUM_ONE, Op[])))
    elseif has_polynomial_body(x)
        side_atoms!(out, c, only(x.args))
    else
        push!(out, (mul_cnum(c, x.coeff), QExpr(x.kind, CNUM_ONE, qexpr_copy_args(x))))
    end
    return out
end

side_atoms(x::QExprArg) = side_atoms!(SideAtom[], CNUM_ONE, x)

atom_qexpr(x::QAdd) = qexpr_from_qadd(x)
atom_qexpr(x::QExpr) = x

struct SuperExprBuilder
    poly::SuperTermDict
    indices::Vector{Index}
    formal::FormalDict
end
SuperExprBuilder() = SuperExprBuilder(SuperTermDict(), Index[], FormalDict())

function add_formal!(d::FormalDict, key::FormalSandwich, c::CNum)
    iszero_cnum(c) && return d
    total = add_cnum(get(d, key, CNUM_ZERO), c)
    iszero_cnum(total) ? delete!(d, key) : (d[key] = total)
    return d
end

function add_atom_pair!(b::SuperExprBuilder, c::CNum, l::QExprArg, r::QExprArg)
    iszero_cnum(c) && return b
    if l isa QAdd && r isa QAdd
        tl, cl = first(l.arguments)
        tr, cr = first(r.arguments)
        canonical_sandwich!(b.poly, tl.ops, tr.ops, mul_cnum(c, mul_cnum(cl, cr)), merge_ne(tl.ne, tr.ne))
    else
        add_formal!(b.formal, FormalSandwich(atom_qexpr(l), atom_qexpr(r)), c)
    end
    return b
end

"""Add `c * left ρ right` with arbitrary polynomial or formal sides."""
function add_sandwich!(b::SuperExprBuilder, c::CNum, left::QExprArg, right::QExprArg)
    iszero_cnum(c) && return b
    for (cl, l) in side_atoms(left), (cr, r) in side_atoms(right)
        add_atom_pair!(b, mul_cnum(c, mul_cnum(cl, cr)), l, r)
    end
    return b
end

function accumulate!(b::SuperExprBuilder, S::SuperOperator, scale::CNum = CNUM_ONE)
    iszero_cnum(scale) && return b
    for (t, c) in S.arguments
        addto_key!(b.poly, copy_key(t), mul_cnum(c, scale))
    end
    merge_indices!(b.indices, S.indices)
    return b
end

function accumulate!(b::SuperExprBuilder, X::SuperExpr, scale::CNum = CNUM_ONE)
    accumulate!(b, X.poly, scale)
    iszero_cnum(scale) && return b
    for (t, c) in X.formal
        add_formal!(b.formal, t, mul_cnum(c, scale))
    end
    return b
end

build(b::SuperExprBuilder) = SuperExpr(superop_pruned(copy(b.poly), copy(b.indices)), copy(b.formal))

## Constructors

as_part(x::QExpr) = x
as_part(x::Union{QField, Coefficient}) = as_side(x)

sandwich_of(A::QAdd, B::QAdd) = sprepost_poly(A, B)
function sandwich_of(A::QExprArg, B::QExprArg)
    b = SuperExprBuilder()
    add_sandwich!(b, CNUM_ONE, A, B)
    return build(b)
end

## Linear arithmetic

function Base.:+(a::SuperExpr, b::SuperExpr)
    builder = SuperExprBuilder()
    accumulate!(builder, a)
    accumulate!(builder, b)
    return build(builder)
end
Base.:+(a::SuperExpr) = a

function scale(X::SuperExpr, c::CNum)
    builder = SuperExprBuilder()
    accumulate!(builder, X, c)
    return build(builder)
end
Base.:-(X::SuperExpr) = scale(X, CNUM_NEG1)
Base.:-(a::SuperExpr, b::SuperExpr) = a + (-b)
Base.:*(X::SuperExpr, c::Coefficient) = scale(X, to_cnum(c))
Base.:*(c::Coefficient, X::SuperExpr) = scale(X, to_cnum(c))
Base.:/(X::SuperExpr, b::Number) =
    b isa Integer && !(b isa Bool) && !iszero(b) ? X * (1 // b) : X * inv(b)
Base.:/(X::SuperExpr, b::Coefficient) = X * inv(b)
Base.://(X::SuperExpr, b::Integer) = X * (1 // b)
Base.://(X::SuperExpr, b::Coefficient) = X / b
Base.:+(X::SuperExpr, c::Coefficient) = X + scale(one(SuperExpr), to_cnum(c))
Base.:+(c::Coefficient, X::SuperExpr) = X + c
Base.:-(X::SuperExpr, c::Coefficient) = X + (-to_cnum(c))
Base.:-(c::Coefficient, X::SuperExpr) = (-X) + c

## Formal Lindbladians

function polynomial_qadd(q::QExpr)
    qexpr_is_scalar(q) && return (true, single_qadd(q.coeff, Op[]))
    iszero(q) && return (true, zero_qadd())
    has_polynomial_body(q) && return (true, only(q.args)::QAdd)
    return (false, zero_qadd())
end

function lower_formal_hamiltonian!(b::SuperExprBuilder, H::QExpr)
    parts = H.kind == QEXPR_ADD ? H.args : QExprArg[H]
    for part in parts
        is_poly, p = part isa QAdd ? (true, part) : polynomial_qadd(part)
        if is_poly
            merge_indices!(b.indices, lower_hamiltonian!(b.poly, p))
        else
            one_side = single_qadd(CNUM_ONE, Op[])
            add_sandwich!(b, CNUM_NEG_IM, part, one_side)
            add_sandwich!(b, CNUM_IM, one_side, part)
        end
    end
    return b
end

# `γ 𝒟[J, K]` over the additive atoms `J = Σ cⱼ jⱼ`, `K = Σ cₖ kₖ`: the loss terms use the
# same atom pairs as the gain term, so `transpose(L)(1)` cancels atom by atom.
function lower_formal_channel!(b::SuperExprBuilder, J::QExpr, K::QExpr, γ::CNum)
    one_side = single_qadd(CNUM_ONE, Op[])
    minus_half = mul_cnum(CNUM_NEG1, mul_cnum(CNUM_HALF, γ))
    for (cj, j) in side_atoms(J), (ck, k) in side_atoms(K)
        c = mul_cnum(cj, conj_cnum(ck))
        kd = adjoint(k)
        add_sandwich!(b, mul_cnum(γ, c), j, kd)
        kj = kd * j
        add_sandwich!(b, mul_cnum(minus_half, c), kj, one_side)
        add_sandwich!(b, mul_cnum(minus_half, c), one_side, kj)
    end
    return b
end

"""
    SuperExpr(L::Lindbladian{QExpr}) -> SuperExpr

Lower a formal Lindbladian to its sandwich form.
"""
function SuperExpr(L::Lindbladian{QExpr})
    b = SuperExprBuilder()
    lower_formal_hamiltonian!(b, L.hamiltonian)
    for (ch, γ) in L.channels
        pj, J = polynomial_qadd(ch.jump)
        pk, K = polynomial_qadd(ch.partner)
        if pj && pk
            lower_channel!(b.poly, J, K, γ, ch.scope, ch.ne)
            merge_indices!(b.indices, ch.scope)
        else
            isempty(ch.scope) || formal_scope_error()
            lower_formal_channel!(b, ch.jump, ch.partner, γ)
        end
    end
    return build(b)
end

## Promotion

const AnyMap = Union{SuperOperator, SuperExpr, Lindbladian}

normal_to_symmetric(::AnyMap) = weyl_map_error()
symmetric_to_normal(::AnyMap) = weyl_map_error()

as_super_expr(X::SuperExpr) = X
as_super_expr(S::SuperOperator) = SuperExpr(S)
as_super_expr(L::Lindbladian{QAdd}) = SuperExpr(SuperOperator(L))
as_super_expr(L::Lindbladian{QExpr}) = SuperExpr(L)

const PolyMap = Union{SuperOperator, Lindbladian{QAdd}}

# Mixed operands promote along `Lindbladian{QAdd} → SuperOperator → SuperExpr`; the
# result type depends only on the operand types.
lift(a::PolyMap, b::PolyMap) = (SuperOperator(a), SuperOperator(b))
lift(a::AnyMap, b::AnyMap) = (as_super_expr(a), as_super_expr(b))

Base.:+(a::AnyMap, b::AnyMap) = +(lift(a, b)...)
Base.:-(a::AnyMap, b::AnyMap) = -(lift(a, b)...)
compose(a::AnyMap, b::AnyMap) = compose(lift(a, b)...)
Base.:*(a::AnyMap, b::AnyMap) = compose(a, b)
Base.:∘(a::AnyMap, b::AnyMap) = compose(a, b)
commutator(a::AnyMap, b::AnyMap) = compose(a, b) - compose(b, a)
anticommutator(a::AnyMap, b::AnyMap) = compose(a, b) + compose(b, a)

Base.:+(L::Lindbladian{QExpr}, c::Coefficient) = SuperExpr(L) + c
Base.:+(c::Coefficient, L::Lindbladian{QExpr}) = SuperExpr(L) + c
Base.:-(L::Lindbladian{QExpr}, c::Coefficient) = SuperExpr(L) - c
Base.:-(c::Coefficient, L::Lindbladian{QExpr}) = c - SuperExpr(L)

## Composition, application, adjoint and transpose

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
        add_sandwich!(b, mul_cnum(ca, cb), al * bl, br * ar)
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
    return build(b)
end

const FormalOperand = Union{SuperExpr, Lindbladian{QExpr}}

Base.inv(::FormalOperand) = negative_power_error()

function Base.:^(X::FormalOperand, n::Integer)
    n >= 0 || negative_power_error()
    S = as_super_expr(X)
    n == 0 && return one(SuperExpr)
    result = S
    for _ in 2:n
        result = compose(result, S)
    end
    return result
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
(L::Lindbladian{QAdd})(ρ::QExpr) = apply_formal(SuperOperator(L), ρ)

function apply_map(X::SuperExpr, ρ::QExprArg)
    terms = QExprArg[]
    iszero(X.poly) || push!(terms, ρ isa QAdd ? X.poly(ρ) : apply_formal(X.poly, ρ))
    for (t, c) in X.formal
        push!(terms, qexpr_product(QExprArg[t.left, ρ, t.right], c))
    end
    return qexpr_sum(terms)
end

(X::SuperExpr)(ρ::Union{QAdd, QExpr}) = apply_map(X, ρ)
(X::SuperExpr)(ρ::QSym) = apply_map(X, single_qadd(CNUM_ONE, Op[ρ]))
(X::SuperExpr)(ρ::Coefficient) = apply_map(X, single_qadd(to_cnum(ρ), Op[]))
(L::Lindbladian{QExpr})(ρ::Union{QSym, QField, Coefficient}) = SuperExpr(L)(ρ)

function Base.adjoint(X::SuperExpr)
    b = SuperExprBuilder()
    accumulate!(b, adjoint(X.poly))
    for (t, c) in X.formal
        add_sandwich!(b, conj_cnum(c), adjoint(t.left), adjoint(t.right))
    end
    return build(b)
end

function Base.transpose(X::SuperExpr)
    b = SuperExprBuilder()
    accumulate!(b, transpose(X.poly))
    for (t, c) in X.formal
        add_formal!(b.formal, FormalSandwich(t.right, t.left), c)
    end
    return build(b)
end

Base.adjoint(L::Lindbladian{QExpr}) = adjoint(SuperExpr(L))
Base.transpose(L::Lindbladian{QExpr}) = transpose(SuperExpr(L))
dagger(X::SuperExpr) = adjoint(X)
qadjoint(X::SuperExpr) = adjoint(X)

## Structural operations

function map_formal(fp::FP, fs::FS, fc::FC, X::SuperExpr) where {FP, FS, FC}
    b = SuperExprBuilder()
    accumulate!(b, fp(X.poly))
    for (t, c) in X.formal
        add_sandwich!(b, fc(c), fs(t.left), fs(t.right))
    end
    return build(b)
end

normal_order(X::SuperExpr) = map_formal(normal_order, normal_order, identity, X)
SymbolicUtils.simplify(X::SuperExpr; kwargs...) = map_formal(
    S -> SymbolicUtils.simplify(S; kwargs...), q -> SymbolicUtils.simplify(q; kwargs...),
    c -> simplify_prefactor(c; kwargs...), X,
)
Symbolics.expand(X::SuperExpr; kwargs...) = map_formal(
    S -> Symbolics.expand(S; kwargs...), q -> Symbolics.expand(q; kwargs...),
    c -> expand_prefactor(c; kwargs...), X,
)
expand_completeness(X::SuperExpr) = map_formal(expand_completeness, expand_completeness, identity, X)
assume_distinct_index(X::SuperExpr, pairs::Vector{Tuple{Index, Index}}) = map_formal(
    S -> assume_distinct_index(S, pairs), q -> assume_distinct_index(q, pairs), identity, X,
)
change_index(X::SuperExpr, from::Index, to::Index) = map_formal(
    S -> change_index(S, from, to), q -> change_index(q, from, to), c -> change_index(c, from, to), X,
)
change_index(X::SuperExpr, pairs::AbstractDict{Index, Index}) = map_formal(
    S -> change_index(S, pairs), q -> change_index(q, pairs), c -> change_index(c, pairs), X,
)

function SymbolicUtils.substitute(X::SuperExpr, rules::AbstractDict; replace_adjoint::Bool = true)
    op_rules, scalar_rules = split_substitution_rules(rules, replace_adjoint)
    phase_rules = nonreal_phase_substitutions(scalar_rules)
    return map_formal(
        S -> iszero(S) ? S : substitute_split(S, op_rules, scalar_rules),
        q -> substitute_part(q, op_rules, scalar_rules, phase_rules),
        c -> substitute_cnum(c, scalar_rules, phase_rules), X,
    )
end

conjugate(X::SuperExpr, U::UnitaryTransform) =
    map_formal(S -> conjugate(S, U), q -> conjugate(q, U), identity, X)
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

function get_operators(X::SuperExpr)
    ops = copy(get_operators(X.poly))
    for t in keys(X.formal)
        append!(ops, get_operators(t.left))
        append!(ops, get_operators(t.right))
    end
    unique!(ops)
    sort!(ops; by = order_key)
    return ops
end

function acts_on(X::SuperExpr)
    out = copy(acts_on(X.poly))
    for t in keys(X.formal)
        append!(out, acts_on(t.left))
        append!(out, acts_on(t.right))
    end
    unique!(sort!(out))
    return out
end

function Symbolics.get_variables(X::SuperExpr; kwargs...)
    vars = copy(Symbolics.get_variables(X.poly; kwargs...))
    for (t, c) in X.formal
        append!(vars, Symbolics.get_variables(t.left; kwargs...))
        append!(vars, Symbolics.get_variables(t.right; kwargs...))
        append!(vars, Symbolics.get_variables(real(c); kwargs...))
        append!(vars, Symbolics.get_variables(imag(c); kwargs...))
    end
    unique!(vars)
    sort!(vars; by = string)
    return vars
end

get_sum_indices(X::SuperExpr) = get_sum_indices(X.poly)

## Taylor lowering

taylor_side(q::QExpr, order::Int) = lower_taylor_qexpr(q, order)

"""
    taylor(X::SuperExpr, ns) -> SuperOperator
    taylor(L::Lindbladian{QExpr}, ns) -> Lindbladian{QAdd}

Lower every formal side to its Maclaurin polynomial through order `n` (`ns = 0:n`, applied
per formal-function node as for `QExpr`).
"""
function Symbolics.taylor(X::SuperExpr, ns::UnitRange{T}) where {T <: Integer}
    order = taylor_prefix_order(ns)
    b = SuperBuilder()
    accumulate!(b, X.poly)
    for (t, c) in X.formal
        accumulate!(b, c, taylor_side(t.left, order), taylor_side(t.right, order))
    end
    return build(b)
end

function Symbolics.taylor(L::Lindbladian{QExpr}, ns::UnitRange{T}) where {T <: Integer}
    order = taylor_prefix_order(ns)
    b = LindbladianBuilder{QAdd}()
    add_hamiltonian!(b.hamiltonian, taylor_side(L.hamiltonian, order))
    for (ch, γ) in L.channels
        J = taylor_side(ch.jump, order)
        K = is_diagonal(ch) ? J : taylor_side(ch.partner, order)
        add_channel!(b, J, K, γ, copy(ch.scope), ch.ne)
    end
    return build(b)
end
