function map_sides(f::F, S::SuperOperator) where {F}
    out = SuperTermDict()
    for (t, c) in S.arguments
        lt = f(t.left, t.ne)
        rt = f(t.right, t.ne)
        temp = SuperTermDict()
        for (a, ca) in lt, (b, cb) in rt
            addto_key!(temp, SuperTerm(a.ops, b.ops, t.ne), c * (ca * cb))
        end
        settle_scope!(out, temp, bound_indices(t, c, S.indices))
    end
    return superop_pruned(out, copy(S.indices))
end

function restream_side(ops::Vector{Op}, ne::Vector{NonEqualPair})
    out = QTermDict()
    stream!(out, copy(ops), CNUM_ONE, ne)
    return out
end

normal_order(S::SuperOperator) = map_sides(restream_side, S)

map_coefficients(f::F, S::SuperOperator) where {F} =
    superop_pruned(mapped_terms(f, S.arguments, S.indices), copy(S.indices))

SymbolicUtils.simplify(S::SuperOperator; kwargs...) =
    map_coefficients(c -> simplify_prefactor(c; kwargs...), normal_order(S))
Symbolics.expand(S::SuperOperator; kwargs...) =
    map_coefficients(c -> expand_prefactor(c; kwargs...), S)
exponential_form(S::SuperOperator) = map_coefficients(exponential_form, S)
trigonometric_form(S::SuperOperator) = map_coefficients(trigonometric_form, S)

function SymbolicUtils.substitute(S::SuperOperator, rules::AbstractDict; replace_adjoint::Bool = true)
    return substitute_with(S, SubstitutionRules(rules, replace_adjoint))
end

function substitute_side(ops::Vector{Op}, ne::Vector{NonEqualPair}, r::SubstitutionRules)
    side = substitute_with(side_qadd(ops, ne), r)
    isempty(side.indices) || summed_side_error("an operator substitution inside a superoperator")
    return side
end

function substitute_with(S::SuperOperator, r::SubstitutionRules)
    (iszero(S) || isempty(r)) && return S
    out = SuperTermDict()
    for (t, c) in S.arguments
        nc = substitute_cnum(c, r)
        iszero(nc) && continue
        lt = substitute_side(t.left, t.ne, r)
        rt = substitute_side(t.right, t.ne, r)
        bound = bound_indices(t, c, S.indices)
        for (a, ca) in lt.arguments, (b, cb) in rt.arguments
            settle_sandwich!(
                out, a.ops, b.ops, nc * (ca * cb),
                merge_ne(t.ne, merge_ne(a.ne, b.ne)), bound,
            )
        end
    end
    return superop_pruned(out, copy(S.indices))
end

function Symbolics.get_variables(S::SuperOperator; kwargs...)
    vars = SymbolicUtils.BasicSymbolic[]
    for c in values(S.arguments)
        append!(vars, Symbolics.get_variables(real(c); kwargs...))
        append!(vars, Symbolics.get_variables(imag(c); kwargs...))
    end
    unique!(vars)
    sort!(vars; by = string)
    return vars
end

function get_operators(S::SuperOperator)
    ops = Op[]
    for t in keys(S.arguments)
        append!(ops, t.left)
        append!(ops, t.right)
    end
    unique!(ops)
    sort!(ops; by = order_key)
    return ops
end

@noinline function weyl_map_error()
    throw(
        ArgumentError(
            "Weyl ordering of a superoperator is not defined: relabelling each side as a " *
                "Weyl symbol changes the map. A phase-space generator needs Bopp operators.",
        ),
    )
end

same_scope(scope::Vector{Index}, ne::Vector{NonEqualPair}) = (true, copy(scope), ne)

"""
    map_lindbladian(fH, fJ, fγ, fscope, P, L) -> Lindbladian{P}

Rebuild `L` from transformed parts, renormalizing every channel. `fscope(scope, ne)` returns
whether the channel is kept, with its new scope and constraints.
"""
function map_lindbladian(
        fH::FH, fJ::FJ, fγ::FG, fscope::FS, ::Type{P}, L::Lindbladian,
    ) where {FH, FJ, FG, FS, P}
    b = LindbladianBuilder{P}()
    add_hamiltonian!(b.hamiltonian, fH(L.hamiltonian))
    for (ch, γ) in L.channels
        kept, scope, ne = fscope(ch.scope, ch.ne)
        kept || continue
        J = fJ(ch.jump)
        K = is_diagonal(ch) ? J : fJ(ch.partner)
        add_channel!(b, J, K, fγ(γ), scope, ne)
    end
    return finish!(b)
end

exponential_form(L::Lindbladian{QAdd}) = rewrite(exponential_form, exponential_form, L)
trigonometric_form(L::Lindbladian{QAdd}) = rewrite(trigonometric_form, trigonometric_form, L)
map_coefficients(f::F, L::Lindbladian{QAdd}) where {F} =
    rewrite(x -> map_coefficients(f, x), c -> to_cnum(f(c)), L)

function assume_distinct_index(L::Lindbladian{O}, pairs::Vector{Tuple{Index, Index}}) where {O}
    distinct(x) = assume_distinct_index(x, pairs)
    distinct_scope(scope, ne) = (true, copy(scope), merge_ne_pairs(ne, pairs))
    return map_lindbladian(distinct, distinct, identity, distinct_scope, O, L)
end

function change_index(L::Lindbladian{O}, m::IndexMap) where {O}
    renames_nothing(m) && return L
    renamed(x) = change_index(x, m)
    renamed_scope(scope, ne) = ne_becomes_contradictory(ne, m) ? (false, scope, ne) :
        (true, unique!(Index[reindex(m, idx) for idx in scope]), substitute_ne(ne, m))
    return map_lindbladian(renamed, renamed, renamed, renamed_scope, O, L)
end

function get_sum_indices(L::Lindbladian)
    out = copy(hamiltonian_indices(L.hamiltonian))
    for ch in keys(L.channels)
        merge_indices!(out, ch.scope)
    end
    return out
end
hamiltonian_indices(H::QAdd) = H.indices
hamiltonian_indices(::QExpr) = EMPTY_INDICES

function get_indices(L::Lindbladian)
    out = get_indices(L.hamiltonian)
    for ch in keys(L.channels)
        merge_indices!(out, ch.scope)
        merge_indices!(out, get_indices(ch.jump))
        merge_indices!(out, get_indices(ch.partner))
    end
    return out
end

function constraint_pairs(L::Lindbladian)
    out = L.hamiltonian isa QAdd ? constraint_pairs(L.hamiltonian) : NonEqualPair[]
    for ch in keys(L.channels), p in ch.ne
        push_ne_unique!(out, p)
    end
    return out
end

function accumulate!(b::SuperBuilder, S::SuperOperator, scale::CNum = CNUM_ONE)
    iszero(scale) && return b
    for (t, c) in S.arguments
        addto_key!(b.args, copy_key(t), scale === CNUM_ONE ? c : c * scale)
    end
    isempty(S.indices) || merge_indices!(b.indices, S.indices)
    return b
end
accumulate!(b::SuperBuilder, L::Lindbladian{QAdd}) = accumulate!(b, SuperOperator(L))

function accumulate!(b::SuperBuilder, c::Coefficient, A::Union{QAdd, QSym, Coefficient}, B::Union{QAdd, QSym, Coefficient})
    cc = to_cnum(c)
    iszero(cc) && return b
    add_sandwiches!(b.args, cc, as_side(A), as_side(B), "`accumulate!`")
    return b
end

build(b::SuperBuilder) = superop_pruned(copy(b.args), copy(b.indices))
finish!(b::SuperBuilder) = superop_pruned(b.args, b.indices)

function Base.sum(a::AbstractArray{SuperOperator})
    isempty(a) && return zero(SuperOperator)
    b = SuperBuilder()
    for x in a
        accumulate!(b, x)
    end
    return finish!(b)
end

function Base.sum(a::AbstractArray{Lindbladian{O}}) where {O}
    isempty(a) && return zero(Lindbladian{O})
    b = LindbladianBuilder{O}()
    for x in a
        accumulate!(b, x)
    end
    return finish!(b)
end

MA.operate(::typeof(sum), a::AbstractArray{SuperOperator}) = sum(a)
MA.operate(::typeof(sum), a::AbstractArray{Lindbladian{O}}) where {O} = sum(a)

"""
    superterm_order_key(t::SuperTerm)

Total ordering key for a sandwich: total length, then the left and right operators by
[`order_key`](@ref), then the canonical constraints.
"""
superterm_order_key(t::SuperTerm) = (
    length(t.left) + length(t.right), map(order_key, t.left), map(order_key, t.right),
    ne_sort_key(t.ne),
)

"""
    superoperator_order_key(S::SuperOperator) -> Vector

Total, reproducible ordering key: the sorted sandwich and coefficient keys. Two maps tie
exactly when they are `isequal` up to their summation scope.
"""
function superoperator_order_key(S::SuperOperator)
    pairs = [(superterm_order_key(t), coeff_key(c)) for (t, c) in S.arguments]
    sort!(pairs)
    return pairs
end

canonical_isless(a::SuperOperator, b::SuperOperator) = isless(
    (superoperator_order_key(a), indices_key(a.indices)), (superoperator_order_key(b), indices_key(b.indices)),
)

function canonical_isless(a::Channel{O}, b::Channel{O}) where {O}
    isequal(a.jump, b.jump) || return canonical_isless(a.jump, b.jump)
    isequal(a.partner, b.partner) || return canonical_isless(a.partner, b.partner)
    sa, sb = map(index_key, a.scope), map(index_key, b.scope)
    sa != sb && return isless(sa, sb)
    return isless(ne_sort_key(a.ne), ne_sort_key(b.ne))
end
