"""
    to_numeric(expr, basis, sites::Dict{Int, Vector{Int}}; parameter, operators)
    to_numeric(expr, basis, sites::Dict{Int, Vector{Int}}[, operators[, parameter]])

Matrix of an expression with sums over indices, on a basis with one factor per site.

Each sum is written out over the sites. `sites` maps the position of each subsystem in the
symbolic space to the positions of its copies in `basis`; a subsystem without an index has a
single copy. Pairs of indices that must differ skip their coinciding sites. `parameter`
gives the values of all parameters. A site-dependent parameter ``g_i`` is keyed by site, with
the value of ``g_1`` under `substitute(g, Dict(index_sym(i) => 1))`. `operators` maps
operators to matrices as for the other forms. This form is available for QuantumOpticsBase
only.

# Examples
```jldoctest
julia> using QuantumOpticsBase

julia> h = FockSpace(:cavity) ⊗ NLevelSpace(:atom, 2);

julia> a = Destroy(h, :a, 1); i = Index(h, :i, 2, 2);

julia> σ(k, l) = IndexedOperator(Transition(h, :σ, k, l, 2), i);

julia> g = IndexedVariable(:g, i); H = Σ(g * (a' * σ(1, 2) + a * σ(2, 1)), i);

julia> b = FockBasis(3) ⊗ NLevelBasis(2) ⊗ NLevelBasis(2);

julia> gs = Dict(substitute(g, Dict(index_sym(i) => k)) => 0.1k for k in 1:2);

julia> size(to_numeric(H, b, Dict(1 => [1], 2 => [2, 3]); parameter = gs))
(16, 16)
```

See also [`Σ`](@ref), [`index_sym`](@ref).
"""
function to_numeric(
        q::QAdd, b::Basis,
        sites::AbstractDict{Int, Vector{Int}},
        d::AbstractDict{<:QSym} = NO_SUBS,
        scalar_subs::AbstractDict = NO_SCALAR_SUBS;
        operators::AbstractDict{<:QSym} = d,
        parameter::AbstractDict = scalar_subs,
    )
    return to_numeric_indexed(QuantumOpticsBackend(), b, q, sites, operators, parameter)
end

function to_numeric_indexed(
        be::NumericBackend, basis, q::QAdd,
        sites::AbstractDict{Int, Vector{Int}},
        d::AbstractDict{<:QSym}, scalar_subs::AbstractDict,
    )
    isempty(q.indices) &&
        return to_numeric(q, basis; operators = d, adjoint_ops = false, parameter = scalar_subs)
    validate_sites(be, basis, sites)
    ctx = NumericContext(be, basis, d, Dict{Int, Vector{Int}}(sites))
    terms = OpTerm[(to_complex(c), ops) for (c, ops) in unrolled_terms(q, ctx.sites, split_scalar_subs(scalar_subs))]
    return static_operator(terms, ctx, EagerOperator())
end

function validate_sites(be::NumericBackend, basis, sites)
    nsubsystems = numeric_num_subsystems(be, basis)
    seen = Set{Int}()
    for (space, slots) in sites
        for slot in slots
            1 <= slot <= nsubsystems || throw(
                ArgumentError(
                    "sites[$space] contains slot $slot, but the numeric basis has " *
                        "$nsubsystems subsystem$(nsubsystems == 1 ? "" : "s")",
                ),
            )
            slot in seen && throw(
                ArgumentError(
                    "numeric subsystem slot $slot occurs more than once in the sites layout",
                ),
            )
            push!(seen, slot)
        end
    end
    return sites
end

function for_each_assignment(f::F, dep::Vector{Index}, ne::Vector{NonEqualPair}, sites::Dict{Int, Vector{Int}}) where {F}
    if isempty(dep)
        f(Dict{Index, Index}(), Dict{Index, Index}())
        return nothing
    end
    lens = Int[length(get(sites, Int(idx.space_index), Int[])) for idx in dep]
    any(iszero, lens) && throw(
        ArgumentError("no sites supplied for the space of a bound summation index"),
    )
    sub_op = Dict{Index, Index}()
    sub_coef = Dict{Index, Index}()
    for combo in 1:prod(lens)
        empty!(sub_op)
        empty!(sub_coef)
        rem = combo - 1
        for k in eachindex(dep)
            kpos = (rem % lens[k]) + 1
            rem ÷= lens[k]
            idx = dep[k]
            sub_op[idx] = Index(idx.name_id, idx.range_id, idx.space_index, Int32(kpos))
            sub_coef[idx] = Index(Int32(0), idx.range_id, idx.space_index, Int32(kpos))
        end
        violates_ne(ne, sub_op) && continue
        f(sub_op, sub_coef)
    end
    return nothing
end

function foreach_unrolled(f::F, X::Union{QAdd, SuperOperator}, sites::Dict{Int, Vector{Int}}, s::SplitScalarRules) where {F}
    for (t, c) in X.arguments
        for_each_assignment(bound_indices(t, c, X.indices), t.ne, sites) do sub_op, sub_coef
            f(t, apply_scalar_subs(change_index(c, sub_coef), s), sub_op)
        end
    end
    return nothing
end

relabel(ops::Vector{Op}, sub::Dict{Index, Index}) = Op[change_index(o, sub) for o in ops]

function unrolled_terms(q::QAdd, sites::Dict{Int, Vector{Int}}, s::SplitScalarRules)
    terms = CoeffTerm[]
    foreach_unrolled(q, sites, s) do t, c, sub_op
        push!(terms, (c, relabel(t.ops, sub_op)))
    end
    return terms
end

function violates_ne(ne::Vector{NonEqualPair}, sub::Dict{Index, Index})
    for (a, b) in ne
        ra = get(sub, a, a)
        rb = get(sub, b, b)
        va = index_slot(ra)
        vb = index_slot(rb)
        iszero(va) && continue
        iszero(vb) && continue
        ra.space_index == rb.space_index || continue
        va == vb && return true
    end
    return false
end

function resolve_slot(op::QSym, sites::AbstractDict{Int, Vector{Int}})
    si = Int(op.space_index)
    slots = get(sites, si, Int[])
    if isempty(slots)
        return si
    end
    if has_index(op.index)
        v = index_slot(op.index)
        if 1 <= v <= length(slots)
            return slots[v]
        end
    end
    length(slots) == 1 && return slots[1]
    throw(
        ArgumentError(
            "cannot resolve slot for operator $(op): index is not a concrete integer site, " *
                "and sites[$si] has $(length(slots)) candidates",
        ),
    )
end
