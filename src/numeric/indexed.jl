# Indexed numeric conversion. The unroll logic (which bound indices a term depends on,
# the `ne` diagonal filter, the `sub_op`/`sub_coef` dual substitution) is backend-neutral
# and kept verbatim; only the leaf emission is routed through the backend hooks at the
# resolved site, so every combination becomes one entry of the same term list the static
# path assembles.

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
        d::AbstractDict{<:QSym} = _NO_SUBS,
        scalar_subs::AbstractDict = _NO_SCALAR_SUBS,
    )
    isempty(q.indices) && return to_numeric(q, b, d)
    be = QuantumOpticsBackend()
    ctx = NumericContext(be, b, d, scalar_subs, Dict{Int, Vector{Int}}(sites))
    sub_re, sub_im, has_imag = _split_scalar_subs(scalar_subs)
    terms = Tuple{ComplexF64, Vector{Any}}[]
    for (term, c) in q.arguments
        _accumulate_indexed_term!(terms, term, c, q.indices, ctx, sub_re, sub_im, has_imag)
    end
    return numeric_assemble(be, b, terms)
end

function _numeric_leaf_indexed(op::Op, ctx::NumericContext)
    haskey(ctx.op_subs, op) && return ctx.op_subs[op]
    slot = _resolve_slot(op, ctx.sites)
    sub = numeric_subbasis(ctx.backend, ctx.basis, slot)
    leaf = numeric_operator(ctx.backend, op, sub)
    return numeric_embed(ctx.backend, ctx.basis, slot, leaf)
end

function _push_indexed_combo!(terms, ops::Vector{Op}, c::CNum, ctx::NumericContext)
    factors = Any[_numeric_leaf_indexed(op, ctx) for op in ops]
    push!(terms, (_to_complex(c), factors))
    return terms
end

function _accumulate_indexed_term!(
        terms, term::QTerm, c::CNum, indices::Vector{Index},
        ctx::NumericContext, sub_re::Dict, sub_im::Dict, has_imag::Bool,
    )
    dep_indices = Index[idx for idx in indices if _depends_on_index_term(c, term.ops, idx)]
    if isempty(dep_indices)
        c_resolved = _apply_scalar_subs(c, sub_re, sub_im, has_imag)
        return _push_indexed_combo!(terms, term.ops, c_resolved, ctx)
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
            # Operators keep the index name (slot kpos drives routing / ne checks and
            # lets a resolved op still match a user `d` key); coefficients use the
            # anonymous name_id-0 form so `index_sym` is Num(kpos), resolving g(i)→g(k).
            sub_op[idx] = Index(idx.name_id, idx.range_id, idx.space_index, Int32(kpos))
            sub_coef[idx] = Index(Int32(0), idx.range_id, idx.space_index, Int32(kpos))
        end
        _violates_ne(term.ne, sub_op) && continue
        new_ops = Op[change_index(op, sub_op) for op in term.ops]
        new_c = change_index(c, sub_coef)
        new_c = _apply_scalar_subs(new_c, sub_re, sub_im, has_imag)
        _push_indexed_combo!(terms, new_ops, new_c, ctx)
    end
    return terms
end

function _violates_ne(ne::Vector{NonEqualPair}, sub::Dict{Index, Index})
    for (a, b) in ne
        ra = get(sub, a, a)
        rb = get(sub, b, b)
        va = index_slot(ra)
        vb = index_slot(rb)
        va === nothing && continue
        vb === nothing && continue
        ra.space_index == rb.space_index || continue
        va == vb && return true
    end
    return false
end

function _resolve_slot(op::QSym, sites::AbstractDict{Int, Vector{Int}})
    si = Int(op.space_index)
    slots = get(sites, si, Int[])
    if isempty(slots)
        return si
    end
    if has_index(op.index)
        v = index_slot(op.index)
        if v !== nothing && 1 <= v <= length(slots)
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
