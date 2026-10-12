function map_sides(f::F, S::SuperOperator) where {F}
    out = SuperTermDict()
    for (t, c) in S.arguments
        lt = f(t.left, t.ne)
        rt = f(t.right, t.ne)
        temp = SuperTermDict()
        for (a, ca) in lt, (b, cb) in rt
            addto_key!(temp, SuperTerm(a.ops, b.ops, t.ne), mul_cnum(c, mul_cnum(ca, cb)))
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

function map_coefficients(f::F, S::SuperOperator) where {F}
    out = SuperTermDict()
    for (t, c) in S.arguments
        nc = f(c)
        iszero_cnum(nc) && continue
        addto_key!(out, copy_key(t), nc)
    end
    return superop_pruned(out, copy(S.indices))
end

SymbolicUtils.simplify(S::SuperOperator; kwargs...) =
    map_coefficients(c -> simplify_prefactor(c; kwargs...), normal_order(S))
Symbolics.expand(S::SuperOperator; kwargs...) =
    map_coefficients(c -> expand_prefactor(c; kwargs...), S)
exponential_form(S::SuperOperator) = map_coefficients(exponential_form, S)
trigonometric_form(S::SuperOperator) = map_coefficients(trigonometric_form, S)

function SymbolicUtils.substitute(S::SuperOperator, rules::AbstractDict; replace_adjoint::Bool = true)
    iszero(S) && return S
    op_rules, scalar_rules = split_substitution_rules(rules, replace_adjoint)
    return substitute_split(S, op_rules, scalar_rules)
end

function substitute_side(ops::Vector{Op}, ne::Vector{NonEqualPair}, op_rules::OpRules, scalar_rules::ScalarRules)
    side = substitute_split(side_qadd(ops, ne), op_rules, scalar_rules)
    isempty(side.indices) || summed_side_error("an operator substitution inside a superoperator")
    return side
end

function substitute_split(S::SuperOperator, op_rules::OpRules, scalar_rules::ScalarRules)
    phase_rules = nonreal_phase_substitutions(scalar_rules)
    out = SuperTermDict()
    for (t, c) in S.arguments
        nc = substitute_cnum(c, scalar_rules, phase_rules)
        iszero_cnum(nc) && continue
        lt = substitute_side(t.left, t.ne, op_rules, scalar_rules)
        rt = substitute_side(t.right, t.ne, op_rules, scalar_rules)
        bound = bound_indices(t, c, S.indices)
        for (a, ca) in lt.arguments, (b, cb) in rt.arguments
            settle_sandwich!(
                out, a.ops, b.ops, mul_cnum(nc, mul_cnum(ca, cb)),
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

"""
    map_lindbladian(fH, fJ, fγ, L) -> Lindbladian

Rebuild `L` from transformed parts, renormalizing every channel.
"""
function map_lindbladian(fH::FH, fJ::FJ, fγ::FG, L::Lindbladian{O}) where {FH, FJ, FG, O}
    b = LindbladianBuilder{O}()
    add_hamiltonian!(b.hamiltonian, fH(L.hamiltonian))
    for (ch, γ) in L.channels
        J = fJ(ch.jump)
        K = is_diagonal(ch) ? J : fJ(ch.partner)
        add_channel!(b, J, K, fγ(γ), copy(ch.scope), ch.ne)
    end
    return build(b)
end

normal_order(L::Lindbladian) = map_lindbladian(normal_order, normal_order, identity, L)
SymbolicUtils.simplify(L::Lindbladian; kwargs...) = map_lindbladian(
    x -> SymbolicUtils.simplify(x; kwargs...), x -> SymbolicUtils.simplify(x; kwargs...),
    c -> simplify_prefactor(c; kwargs...), L,
)
Symbolics.expand(L::Lindbladian; kwargs...) = map_lindbladian(
    x -> Symbolics.expand(x; kwargs...), x -> Symbolics.expand(x; kwargs...),
    c -> expand_prefactor(c; kwargs...), L,
)
expand_completeness(L::Lindbladian) =
    map_lindbladian(expand_completeness, expand_completeness, identity, L)
exponential_form(L::Lindbladian{QAdd}) =
    map_lindbladian(exponential_form, exponential_form, exponential_form, L)
trigonometric_form(L::Lindbladian{QAdd}) =
    map_lindbladian(trigonometric_form, trigonometric_form, trigonometric_form, L)

function assume_distinct_index(L::Lindbladian{O}, pairs::Vector{Tuple{Index, Index}}) where {O}
    b = LindbladianBuilder{O}()
    add_hamiltonian!(b.hamiltonian, assume_distinct_index(L.hamiltonian, pairs))
    for (ch, γ) in L.channels
        J = assume_distinct_index(ch.jump, pairs)
        K = is_diagonal(ch) ? J : assume_distinct_index(ch.partner, pairs)
        add_channel!(b, J, K, γ, copy(ch.scope), merge_ne_pairs(ch.ne, pairs))
    end
    return build(b)
end

function SymbolicUtils.substitute(L::Lindbladian{O}, rules::AbstractDict; replace_adjoint::Bool = true) where {O}
    op_rules, scalar_rules = split_substitution_rules(rules, replace_adjoint)
    phase_rules = nonreal_phase_substitutions(scalar_rules)
    return map_lindbladian(
        x -> substitute_part(x, op_rules, scalar_rules, phase_rules),
        x -> substitute_part(x, op_rules, scalar_rules, phase_rules),
        c -> substitute_cnum(c, scalar_rules, phase_rules), L,
    )
end

substitute_part(q::QAdd, op_rules::OpRules, scalar_rules::ScalarRules, ::PhaseRules) =
    iszero(q) ? q : substitute_split(q, op_rules, scalar_rules)
substitute_part(q::QExpr, op_rules::OpRules, scalar_rules::ScalarRules, phase_rules::PhaseRules) =
    iszero(q) ? q : substitute_qexpr_split(q, op_rules, scalar_rules, phase_rules)

change_index(L::Lindbladian, from::Index, to::Index) = change_index(L, Dict(from => to))

function change_index(L::Lindbladian{O}, pairs::AbstractDict{Index, Index}) where {O}
    isempty(pairs) && return L
    b = LindbladianBuilder{O}()
    add_hamiltonian!(b.hamiltonian, change_index(L.hamiltonian, pairs))
    for (ch, γ) in L.channels
        ne_becomes_contradictory(ch.ne, pairs) && continue
        J = change_index(ch.jump, pairs)
        K = is_diagonal(ch) ? J : change_index(ch.partner, pairs)
        sc = unique!(Index[get(pairs, idx, idx) for idx in ch.scope])
        add_channel!(b, J, K, change_index(γ, pairs), sc, substitute_ne(ch.ne, pairs))
    end
    return build(b)
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

function get_operators(L::Lindbladian)
    ops = copy(get_operators(L.hamiltonian))
    for ch in keys(L.channels)
        append!(ops, get_operators(ch.jump))
        append!(ops, get_operators(ch.partner))
    end
    unique!(ops)
    sort!(ops; by = order_key)
    return ops
end

function acts_on(L::Lindbladian)
    out = copy(acts_on(L.hamiltonian))
    for ch in keys(L.channels)
        append!(out, acts_on(ch.jump))
        append!(out, acts_on(ch.partner))
    end
    unique!(sort!(out))
    return out
end

function Symbolics.get_variables(L::Lindbladian; kwargs...)
    vars = SymbolicUtils.BasicSymbolic[]
    append!(vars, Symbolics.get_variables(L.hamiltonian; kwargs...))
    for (ch, γ) in L.channels
        append!(vars, Symbolics.get_variables(ch.jump; kwargs...))
        append!(vars, Symbolics.get_variables(ch.partner; kwargs...))
        append!(vars, Symbolics.get_variables(real(γ); kwargs...))
        append!(vars, Symbolics.get_variables(imag(γ); kwargs...))
    end
    unique!(vars)
    sort!(vars; by = string)
    return vars
end

function constraint_pairs(L::Lindbladian)
    out = L.hamiltonian isa QAdd ? constraint_pairs(L.hamiltonian) : NonEqualPair[]
    for ch in keys(L.channels), p in ch.ne
        push_ne_unique!(out, p)
    end
    return out
end

"""
    SuperBuilder()

Mutable accumulator for building a [`SuperOperator`](@ref) by repeated addition in linear
time. Add maps with `accumulate!(b, S)`, single sandwiches `c * A ρ B` with
`accumulate!(b, c, A, B)`, and finish with `build(b)`.
"""
struct SuperBuilder
    args::SuperTermDict
    indices::Vector{Index}
end
SuperBuilder() = SuperBuilder(SuperTermDict(), Index[])

function accumulate!(b::SuperBuilder, S::SuperOperator)
    for (t, c) in S.arguments
        addto_key!(b.args, copy_key(t), c)
    end
    isempty(S.indices) || merge_indices!(b.indices, S.indices)
    return b
end
accumulate!(b::SuperBuilder, L::Lindbladian{QAdd}) = accumulate!(b, SuperOperator(L))

function accumulate!(b::SuperBuilder, c::Coefficient, A::Union{QField, Coefficient}, B::Union{QField, Coefficient})
    cc = to_cnum(c)
    iszero_cnum(cc) && return b
    Aq = as_side(A)
    Bq = as_side(B)
    isempty(Aq.indices) || summed_side_error("`accumulate!`")
    isempty(Bq.indices) || summed_side_error("`accumulate!`")
    for (ta, ca) in Aq.arguments, (tb, cb) in Bq.arguments
        canonical_sandwich!(b.args, ta.ops, tb.ops, mul_cnum(cc, mul_cnum(ca, cb)), merge_ne(ta.ne, tb.ne))
    end
    return b
end

build(b::SuperBuilder) = superop_pruned(copy(b.args), copy(b.indices))

function Base.sum(a::AbstractArray{SuperOperator})
    isempty(a) && return zero(SuperOperator)
    b = SuperBuilder()
    for x in a
        accumulate!(b, x)
    end
    return build(b)
end

function Base.sum(a::AbstractArray{Lindbladian{O}}) where {O}
    isempty(a) && return zero(Lindbladian{O})
    b = LindbladianBuilder{O}()
    for x in a
        accumulate!(b, x)
    end
    return build(b)
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

"""
    channel_order_key(c::Channel)

Total ordering key for a channel: its jump and partner, then its scope and constraints.
"""
channel_order_key(c::Channel{QAdd}) = (
    qadd_order_key(c.jump), qadd_order_key(c.partner),
    [(Int(i.space_index), name_rank(i.name_id)) for i in c.scope], ne_sort_key(c.ne),
)
channel_order_key(c::Channel{QExpr}) = (
    sprint(show, c.jump), sprint(show, c.partner),
    [(Int(i.space_index), name_rank(i.name_id)) for i in c.scope], ne_sort_key(c.ne),
)
