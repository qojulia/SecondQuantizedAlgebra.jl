"""
    to_numeric(q::QAdd, basis, sites::Dict{Int, Vector{Int}}[, d[, scalar_subs]])

Convert a `QAdd` that carries bound summation indices (`q.indices`) to a numeric operator
on a `basis` whose layout replicates one or more abstract subspaces. `sites[orig_space_index]`
is the list of slots that realise that subspace; non-indexed subspaces map to a single-slot
entry.

For each term, every bound index the term depends on is unrolled over the length of its
`sites` vector, with `term.ne` constraints filtering out diagonal combinations. Concrete-site
operators are routed to `sites[op.space_index][k]`, where `k` is the substituted integer
site. Terms that do not depend on any bound index are emitted once as written, matching the
`Σ` convention that i-independent residuals already carry their range factor.

`d` substitutes individual `QSym` operators with custom numeric operators. `scalar_subs`
substitutes symbolic scalar parameters (e.g. `Δ` or `g(k)` for the integer `k` produced by
index unrolling) with numeric values.
"""
function to_numeric(
        q::QAdd, b::Basis,
        sites::AbstractDict{Int, Vector{Int}},
        d::AbstractDict{<:QSym} = NO_SUBS,
        scalar_subs::AbstractDict = NO_SCALAR_SUBS,
    )
    return to_numeric_indexed(QuantumOpticsBackend(), b, q, sites, d, scalar_subs)
end

function to_numeric_indexed(
        be::NumericBackend, basis, q::QAdd,
        sites::AbstractDict{Int, Vector{Int}},
        d::AbstractDict{<:QSym}, scalar_subs::AbstractDict,
    )
    isempty(q.indices) && return to_numeric(q, basis, d)
    validate_sites(be, basis, sites)
    ctx = NumericContext(be, basis, d, scalar_subs, Dict{Int, Vector{Int}}(sites))
    sub_re, sub_im, has_imag = split_scalar_subs(scalar_subs)
    terms = OpTerm[]
    for (term, c) in q.arguments
        accumulate_indexed_term!(terms, term, c, q.indices, ctx, sub_re, sub_im, has_imag)
    end
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

function push_indexed_combo!(terms::Vector{OpTerm}, ops::Vector{Op}, c::CNum)
    push!(terms, (to_complex(c), ops))
    return terms
end

function accumulate_indexed_term!(
        terms::Vector{OpTerm}, term::QTerm, c::CNum, indices::Vector{Index},
        ctx::NumericContext, sub_re::ScalarRules, sub_im::ScalarRules, has_imag::Bool,
    )
    dep_indices = bound_indices(term, c, indices)
    if isempty(dep_indices)
        c_resolved = apply_scalar_subs(c, sub_re, sub_im, has_imag)
        return push_indexed_combo!(terms, term.ops, c_resolved)
    end
    lens = Int[length(ctx.sites[Int(idx.space_index)]) for idx in dep_indices]
    total = prod(lens)
    sub_op = Dict{Index, Index}()
    sub_coef = Dict{Index, Index}()
    for combo in 1:total
        empty!(sub_op)
        empty!(sub_coef)
        rem = combo - 1
        for k in 1:length(dep_indices)
            kpos = (rem % lens[k]) + 1
            rem ÷= lens[k]
            idx = dep_indices[k]
            sub_op[idx] = Index(idx.name_id, idx.range_id, idx.space_index, Int32(kpos))
            sub_coef[idx] = Index(Int32(0), idx.range_id, idx.space_index, Int32(kpos))
        end
        violates_ne(term.ne, sub_op) && continue
        new_ops = Op[change_index(op, sub_op) for op in term.ops]
        new_c = change_index(c, sub_coef)
        new_c = apply_scalar_subs(new_c, sub_re, sub_im, has_imag)
        push_indexed_combo!(terms, new_ops, new_c)
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
