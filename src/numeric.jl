const QuantumState = Union{QuantumOpticsBase.StateVector, QuantumOpticsBase.AbstractOperator}

# `Union{}` value keeps the haskey-fallback branch concrete-typed.
const _NO_SUBS = Dict{QSym, Union{}}()

"""
    to_numeric(op, basis [, d::AbstractDict{<:QSym}])
    to_numeric(op, state [, d::AbstractDict{<:QSym}])
    to_numeric(op, basis; parameter=Dict(), time_parameter=Dict(), operators=Dict(),
               adjoint_ops=true, op_type=sparse)
    to_numeric(ops::AbstractVector, basis; kwargs...)

Convert a symbolic operator expression to a numeric QuantumOpticsBase operator.
`d` substitutes individual `QSym`s with custom numeric operators. Throws
`ArgumentError` if any prefactor cannot be reduced to a concrete `ComplexF64`.

The keyword form first substitutes scalar `parameter`s, then translates using
custom numeric `operators`. Missing adjoint entries are added to `operators`
when `adjoint_ops=true`. If `time_parameter` is non-empty, values may be numbers
or functions of time and the result is a closure `t -> op(t)`.

The result representation is controlled by `op_type` and does not depend on the
shape of the expression: a bare operator, a product, and a sum all return the
same operator type. `op_type` is applied once to the assembled operator:

  - `op_type=sparse` (default): a sparse `Operator`. A good general default: it
    is memory-scalable, type-stable, and the form the QuantumOptics solvers
    consume most efficiently.
  - `op_type=dense`: a dense `Operator`. Prefer for small Hilbert spaces.
  - `op_type=identity`: the natural lazy representation (`LazyTensor` /
    `LazyProduct` / `LazySum`). Opt in for large tensor-product spaces where an
    operator is local to a few subsystems, so the full Kronecker product is
    never materialised. (In the time-dependent path with more than one term,
    terms are still materialised to keep the returned closure type-stable.)

The positional forms `to_numeric(op, basis[, d])` always return `sparse`; use
the keyword form to select `dense` or `identity`.

# Examples

```jldoctest
julia> using QuantumOpticsBase

julia> h = FockSpace(:f);

julia> @qnumbers a::Destroy(h);

julia> b = FockBasis(5); ψ = fockstate(b, 2);

julia> real(QuantumOpticsBase.expect(to_numeric(a' * a, ψ), ψ)) ≈ 2
true
```

See also [`numeric_average`](@ref).
"""
function to_numeric end

function to_numeric(op::QField, b::QuantumOpticsBase.Basis; kwargs...)
    return _to_numeric_kw(op, b; kwargs...)
end

function to_numeric(x::Union{Number, SymbolicUtils.BasicSymbolic}, b::QuantumOpticsBase.Basis; kwargs...)
    return _to_numeric_kw(x, b; kwargs...)
end

function to_numeric(op, state::QuantumState; kwargs...)
    return to_numeric(op, QuantumOpticsBase.basis(state); kwargs...)
end

function to_numeric(ops::AbstractVector, b::QuantumOpticsBase.Basis; kwargs...)
    return [to_numeric(op, b; kwargs...) for op in ops]
end

# Per-basis kind-chain on the single concrete `Op`.
function to_numeric(op::Op, b::QuantumOpticsBase.FockBasis)
    is_destroy(op)  && return QuantumOpticsBase.destroy(b)
    is_create(op)   && return QuantumOpticsBase.create(b)
    is_position(op) && return (QuantumOpticsBase.destroy(b) + QuantumOpticsBase.create(b)) / sqrt(2)
    is_momentum(op) && return im * (QuantumOpticsBase.create(b) - QuantumOpticsBase.destroy(b)) / sqrt(2)
    throw(ArgumentError("Op kind $(op.kind) does not act on a FockBasis"))
end
function to_numeric(op::Op, b::QuantumOpticsBase.NLevelBasis)
    is_transition(op) && return QuantumOpticsBase.transition(b, Int(op.l1), Int(op.l2))
    throw(ArgumentError("Op kind $(op.kind) does not act on an NLevelBasis"))
end
function to_numeric(op::Op, b::QuantumOpticsBase.SpinBasis)
    if is_pauli(op)
        b.spinnumber == 1 // 2 || throw(ArgumentError("Pauli operators require SpinBasis(1//2), got SpinBasis($(b.spinnumber))"))
        op.l1 == 1 && return QuantumOpticsBase.sigmax(b)
        op.l1 == 2 && return QuantumOpticsBase.sigmay(b)
        return QuantumOpticsBase.sigmaz(b)
    end
    if is_spin(op)
        op.l1 == 1 && return 0.5 * QuantumOpticsBase.sigmax(b)
        op.l1 == 2 && return 0.5 * QuantumOpticsBase.sigmay(b)
        return 0.5 * QuantumOpticsBase.sigmaz(b)
    end
    throw(ArgumentError("Op kind $(op.kind) does not act on a SpinBasis"))
end
# ── Internal lazy embedding ─────────────────────────────────────────────────
# `_embed` builds the natural lazy representation of a single operator on `b`
# (a `LazyTensor` on a `CompositeBasis`). It is the composable building block for
# products and sums; the public `to_numeric` methods materialise the assembled
# result via `op_type` (sparse by default). Keeping the embedding lazy lets a
# whole term be assembled as lazy factors and converted to a concrete operator
# exactly once, at the top.
_embed(op::Op, b::QuantumOpticsBase.Basis) = to_numeric(op, b)
function _embed(op::Op, b::QuantumOpticsBase.CompositeBasis)
    si = Int(op.space_index)
    return QuantumOpticsBase.LazyTensor(b, [si], (_embed(op, b.bases[si]),))
end
function _embed(op::Op, b::QuantumOpticsBase.Basis, d::AbstractDict{<:QSym})
    haskey(d, op) && return d[op]
    return _embed(op, b)
end

# Public positional forms return a concrete sparse operator, independent of
# expression shape. To choose the representation (`dense`, or `identity` for the
# lazy form) use the keyword form `to_numeric(op, b; op_type=...)`.
to_numeric(op::Op, b::QuantumOpticsBase.CompositeBasis) =
    QuantumOpticsBase.sparse(_embed(op, b))
to_numeric(op::Op, b::QuantumOpticsBase.Basis, d::AbstractDict{<:QSym}) =
    QuantumOpticsBase.sparse(_embed(op, b, d))
to_numeric(s::QAdd, b::QuantumOpticsBase.Basis, d::AbstractDict{<:QSym} = _NO_SUBS) =
    QuantumOpticsBase.sparse(_assemble(s, b, d))

# Natural (lazy) assembly of an operator sum; the caller materialises it.
# Scalar-term and operator-product branches return different operator types
# (Diagonal identity vs SparseMatrixCSC), so the conditional stays in-loop.
function _assemble(s::QAdd, b::QuantumOpticsBase.Basis, d::AbstractDict{<:QSym} = _NO_SUBS)
    iter = iterate(s.arguments)
    iter === nothing && return _to_complex(_CNUM_ZERO) * _lazy_one(b)
    (term, c), st = iter
    result = isempty(term.ops) ? _to_complex(c) * _lazy_one(b) : _to_complex(c) * _product(term.ops, b, d)
    while true
        next = iterate(s.arguments, st)
        next === nothing && return result
        (term, c), st = next
        result += isempty(term.ops) ? _to_complex(c) * _lazy_one(b) : _to_complex(c) * _product(term.ops, b, d)
    end
    return
end

function _product(ops::Vector{Op}, b::QuantumOpticsBase.Basis, d::AbstractDict{<:QSym} = _NO_SUBS)
    acc = _embed(first(ops), b, d)
    for i in 2:length(ops)
        acc *= _embed(ops[i], b, d)
    end
    return acc
end

# Lazy (un-materialised) numeric form. Used where the consumer contracts the
# operator with a state immediately and therefore never needs a concrete matrix,
# e.g. `expect` against a `LazyKet` (which has no `expect` method for a sparse
# operator, and whose whole purpose is to avoid materialising the Kronecker
# product). Mirrors the positional `to_numeric(op, b, d)` dispatch, minus the
# final `sparse`.
_to_numeric_lazy(op::Op, b::QuantumOpticsBase.Basis, d::AbstractDict{<:QSym}) = _embed(op, b, d)
_to_numeric_lazy(s::QAdd, b::QuantumOpticsBase.Basis, d::AbstractDict{<:QSym}) = _assemble(s, b, d)

function _to_numeric_kw(
        op,
        b::QuantumOpticsBase.Basis;
        parameter = Dict(),
        time_parameter = Dict(),
        operators = Dict{QSym, Any}(),
        adjoint_ops = true,
        op_type = QuantumOpticsBase.sparse,
    )
    param = _expand_parameter(parameter)
    tp = _normalize_time_parameter(time_parameter)
    ops = _numeric_operator_dict(operators, adjoint_ops)
    return _to_numeric_translated(
        op, b; parameter = param, time_parameter = tp, operators = ops, op_type,
    )
end

function _to_numeric_translated(
        op::QSym,
        b::QuantumOpticsBase.Basis;
        parameter,
        time_parameter,
        operators,
        op_type,
    )
    return _to_numeric_translated(
        _single_qadd(_CNUM_ONE, Op[op]), b;
        parameter, time_parameter, operators, op_type,
    )
end

function _to_numeric_translated(
        op::QAdd,
        b::QuantumOpticsBase.Basis;
        parameter,
        time_parameter,
        operators,
        op_type,
    )
    if isempty(time_parameter)
        return _to_numeric_static(op, b; parameter, operators, op_type)
    end

    op_ = substitute(op, parameter)
    if iszero(op_)
        z = op_type(_to_complex(_CNUM_ZERO) * _lazy_one(b))
        return t -> z
    end

    pairs = collect(op_)
    if length(pairs) == 1
        term, c = pairs[1]
        return _translate_term(
            term.ops, to_num(c), b, time_parameter, operators, op_type,
        )
    end

    # Use a concrete operator representation for every term so the tuple of
    # wrapped closures has one return type, even for scalar-plus-operator sums.
    term_op_type = op_type === identity ? QuantumOpticsBase.sparse : op_type
    first_res = _translate_term(
        pairs[1].first.ops, to_num(pairs[1].second), b,
        time_parameter, operators, term_op_type,
    )
    OpType = typeof(first_res(0.0))
    FW = FunctionWrapper{OpType, Tuple{Float64}}
    wrapped = ntuple(length(pairs)) do k
        res = k == 1 ? first_res : _translate_term(
                pairs[k].first.ops, to_num(pairs[k].second), b,
                time_parameter, operators, term_op_type,
            )
        FW(res)
    end

    return t -> begin
        tt = Float64(t)
        result = wrapped[1](tt)
        for i in 2:length(wrapped)
            result = result + wrapped[i](tt)
        end
        result
    end
end

function _to_numeric_translated(
        arg,
        b::QuantumOpticsBase.Basis;
        parameter,
        time_parameter,
        operators,
        op_type,
    )
    arg_sub = substitute(arg, parameter)
    one_b = _lazy_one(b)
    c = _as_cnum(arg_sub)

    if isempty(time_parameter)
        _coeff_is_const(c) || throw(
            ArgumentError(
                "cannot translate symbolic scalar `$arg_sub` without a value: supply it via " *
                    "`parameter` or `time_parameter`.",
            ),
        )
        return op_type(_const_coeff(c) * one_b)
    end

    if _coeff_is_const(c)
        val = _const_coeff(c)
        op = op_type(val * one_b)
        return t -> op
    end
    basevars, valuefuncs = _time_basis(time_parameter)
    _check_time_variables(c, basevars)
    pref = _compile_coeff(c, basevars...)
    return t -> op_type(pref(map(f -> f(t), valuefuncs)...) * one_b)
end

function _to_numeric_static(
        op::QAdd,
        b::QuantumOpticsBase.Basis;
        parameter,
        operators,
        op_type,
    )
    op_ = substitute(op, parameter)
    iszero(op_) && return op_type(_to_complex(_CNUM_ZERO) * _lazy_one(b))

    result = nothing
    for (term, c_) in op_
        c = to_num(c_)
        _coeff_is_const(c) || throw(
            ArgumentError(
                "cannot translate `$op` to a static operator: coefficient `$c` still " *
                    "depends on a free variable; supply it via `parameter` or `time_parameter`.",
            ),
        )
        contrib = _const_coeff(c) * _numeric_product(term.ops, b, operators, identity)
        result = result === nothing ? contrib : result + contrib
    end
    return op_type(result)
end

_to_numeric_with_ops(op, b::QuantumOpticsBase.Basis, operators::AbstractDict{<:QSym}) =
    isempty(operators) ? _embed(op, b) : _embed(op, b, operators)

function _numeric_product(
        ops::Vector{Op},
        b::QuantumOpticsBase.Basis,
        operators::AbstractDict{<:QSym},
        op_type,
    )
    if isempty(ops)
        return op_type(_lazy_one(b))
    end
    acc = _to_numeric_with_ops(ops[1], b, operators)
    for i in 2:length(ops)
        acc *= _to_numeric_with_ops(ops[i], b, operators)
    end
    return op_type(acc)
end

function _numeric_operator_dict(operators, adjoint_ops::Bool)
    out = Dict{QSym, Any}()
    for (k, v) in operators
        k isa QSym || throw(ArgumentError("operator substitution key `$k` is not a QSym"))
        out[k] = v
    end
    if adjoint_ops
        for (k, v) in operators
            k isa QSym || continue
            k_adj = Base.adjoint(k)
            haskey(out, k_adj) || (out[k_adj] = Base.adjoint(v))
        end
    end
    return out
end

function _expand_parameter(parameter)
    isempty(parameter) && return parameter
    out = Dict{Any, Any}()
    for (k, v) in parameter
        if k isa Complex
            out[real(k)] = real(v)
            ik = imag(k)
            SymbolicUtils.unwrap(ik) isa Number || (out[ik] = imag(v))
        else
            out[k] = v
        end
    end
    return out
end

function _normalize_time_parameter(time_parameter)
    isempty(time_parameter) && return time_parameter
    out = Dict{Any, Any}()
    for (k, v) in time_parameter
        out[k] = v isa Number ? (t -> v) : v
    end
    return out
end

function _time_basis(time_parameter)
    basevars = Any[]
    valuefuncs = Any[]
    for (k, f) in time_parameter
        uk = SymbolicUtils.unwrap(k)
        if SymbolicUtils.issym(uk)
            push!(basevars, k)
            push!(valuefuncs, f)
            continue
        end
        vs = collect(Symbolics.get_variables(k))
        length(vs) == 1 || throw(
            ArgumentError(
                "time_parameter key `$k` must depend on exactly one variable, got $(length(vs)).",
            ),
        )
        wv = Symbolics.wrap(vs[1])
        if isequal(uk, SymbolicUtils.unwrap(conj(wv)))
            push!(basevars, wv)
            push!(valuefuncs, t -> conj(f(t)))
        else
            throw(
                ArgumentError(
                    "unsupported time_parameter key `$k`; only a bare variable `v` or `conj(v)` are supported.",
                ),
            )
        end
    end
    return basevars, Tuple(valuefuncs)
end

_coeff_is_const(c::Complex{Num}) =
    isempty(Symbolics.get_variables(real(c))) && isempty(Symbolics.get_variables(imag(c)))

_const_coeff(c::Complex{Num}) = _to_complex(c)

_as_cnum(x::Complex) = Complex{Num}(Num(real(x)), Num(imag(x)))
_as_cnum(x) = Complex{Num}(Num(x), Num(false))

function _coefficient_variables(c::Complex{Num})
    vars = Any[]
    append!(vars, Symbolics.get_variables(real(c)))
    append!(vars, Symbolics.get_variables(imag(c)))
    unique!(vars)
    return vars
end

function _check_time_variables(c::Complex{Num}, basevars)
    vars = _coefficient_variables(c)
    base_unwrapped = SymbolicUtils.unwrap.(basevars)
    missing = Any[]
    for v in vars
        any(b -> isequal(v, b), base_unwrapped) || push!(missing, v)
    end
    isempty(missing) && return nothing
    throw(
        ArgumentError(
            "time-dependent coefficient `$c` depends on variables without time values: " *
                join(string.(missing), ", "),
        ),
    )
end

function _compile_coeff(c::Complex{Num}, vars...)
    f_re = build_function(real(c), vars...; expression = Val(false))
    f_im = build_function(imag(c), vars...; expression = Val(false))
    g_re = f_re isa Tuple ? first(f_re) : f_re
    g_im = f_im isa Tuple ? first(f_im) : f_im
    return (vals...) -> g_re(vals...) + im * g_im(vals...)
end

function _translate_term(
        ops::Vector{Op},
        c::Complex{Num},
        b::QuantumOpticsBase.Basis,
        time_parameter,
        operators,
        op_type,
    )
    prodop = _numeric_product(ops, b, operators, op_type)
    if _coeff_is_const(c)
        op = _const_coeff(c) * prodop
        return t -> op
    end
    basevars, valuefuncs = _time_basis(time_parameter)
    _check_time_variables(c, basevars)
    pref = _compile_coeff(c, basevars...)
    return t -> pref(map(f -> f(t), valuefuncs)...) * prodop
end

const _NO_SCALAR_SUBS = Dict{Num, Any}()

"""
    to_numeric(q::QAdd, b::CompositeBasis, sites::Dict{Int, Vector{Int}}[, d[, scalar_subs]])

Convert a `QAdd` that carries bound summation indices (`q.indices`) to a numeric
operator on a `CompositeBasis` whose layout replicates one or more abstract
subspaces. `sites[orig_space_index]` is the list of slots into `b.bases` that
realize that subspace; non-indexed subspaces map to a single-slot entry.

For each term, every bound index that the term depends on is unrolled over the
length of its `sites` vector, with `term.ne` constraints filtering out diagonal
combinations. Concrete-site operators are routed to `sites[op.space_index][k]`,
where `k` is the substituted integer site. Terms that do not depend on any
bound index are emitted once as written, matching the `Σ` convention that
i-independent residuals already carry their range factor.

`d` substitutes individual `QSym` operators with custom numeric operators (same
semantics as the single-basis-slot path). `scalar_subs` substitutes symbolic
scalar parameters (e.g. `Δ` or `g(k)` for the integer `k` produced by index
unrolling) with numeric values; supply this to resolve indexed parameters that
only become fully concrete after the symbolic sum is unrolled.
"""
function to_numeric(
        q::QAdd, b::QuantumOpticsBase.CompositeBasis,
        sites::AbstractDict{Int, Vector{Int}},
        d::AbstractDict{<:QSym} = _NO_SUBS,
        scalar_subs::AbstractDict = _NO_SCALAR_SUBS,
    )
    if isempty(q.indices)
        return to_numeric(q, b, d)
    end
    sub_re, sub_im, has_imag = _split_scalar_subs(scalar_subs)
    result = _to_complex(_CNUM_ZERO) * _lazy_one(b)
    for (term, c) in q.arguments
        result = _accumulate_indexed_term!(
            result, term, c, q.indices, b, sites, d, sub_re, sub_im, has_imag,
        )
    end
    return QuantumOpticsBase.sparse(result)
end

function _accumulate_indexed_term!(
        acc, term::QTerm, c::CNum, indices::Vector{Index},
        b::QuantumOpticsBase.CompositeBasis,
        sites::AbstractDict{Int, Vector{Int}},
        d::AbstractDict{<:QSym},
        sub_re::Dict, sub_im::Dict, has_imag::Bool,
    )
    dep_indices = Index[idx for idx in indices if _depends_on_index_term(c, term.ops, idx)]
    if isempty(dep_indices)
        c_resolved = _apply_scalar_subs(c, sub_re, sub_im, has_imag)
        return acc + _emit_indexed_combo(term.ops, c_resolved, b, sites, d)
    end
    lens = Int[length(sites[Int(idx.space_index)]) for idx in dep_indices]
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
        if _violates_ne(term.ne, sub_op)
            continue
        end
        new_ops = Op[change_index(op, sub_op) for op in term.ops]
        new_c = change_index(c, sub_coef)
        new_c = _apply_scalar_subs(new_c, sub_re, sub_im, has_imag)
        acc = acc + _emit_indexed_combo(new_ops, new_c, b, sites, d)
    end
    return acc
end

# Split user-supplied scalar substitutions into real and imag legs once per
# `to_numeric` call. A complex RHS `re + im*ii` propagates as separate real
# and imaginary substitutions, preserving the `Complex{Num}` invariant.
function _split_scalar_subs(scalar_subs::AbstractDict)
    sub_re = Dict{Any, Any}()
    sub_im = Dict{Any, Any}()
    has_imag = false
    for (k, v) in scalar_subs
        kraw = SymbolicUtils.unwrap(k)
        if v isa Complex
            sub_re[kraw] = real(v)
            sub_im[kraw] = imag(v)
            iszero(imag(v)) || (has_imag = true)
        else
            sub_re[kraw] = v
            sub_im[kraw] = 0
        end
    end
    return sub_re, sub_im, has_imag
end

function _apply_scalar_subs(c::CNum, sub_re::Dict, sub_im::Dict, has_imag::Bool)
    (isempty(sub_re) || _is_native(c)) && return c
    re_part = SymbolicUtils.unwrap(real(c))
    im_part = SymbolicUtils.unwrap(imag(c))
    if !has_imag
        return _cnum(
            Num(Symbolics.substitute(re_part, sub_re)),
            Num(Symbolics.substitute(im_part, sub_re)),
        )
    end
    # (re + i*im) -> (re|sub_re - im|sub_im) + i*(re|sub_im + im|sub_re)
    new_re = Symbolics.substitute(re_part, sub_re) -
        Symbolics.substitute(im_part, sub_im)
    new_im = Symbolics.substitute(re_part, sub_im) +
        Symbolics.substitute(im_part, sub_re)
    return _cnum(Num(new_re), Num(new_im))
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

function _emit_indexed_combo(
        ops::Vector{Op}, c::CNum,
        b::QuantumOpticsBase.CompositeBasis,
        sites::AbstractDict{Int, Vector{Int}},
        d::AbstractDict{<:QSym},
    )
    if isempty(ops)
        return _to_complex(c) * _lazy_one(b)
    end
    acc = _site_routed_op(ops[1], b, sites, d)
    for k in 2:length(ops)
        acc = acc * _site_routed_op(ops[k], b, sites, d)
    end
    return _to_complex(c) * acc
end

function _site_routed_op(
        op::QSym, b::QuantumOpticsBase.CompositeBasis,
        sites::AbstractDict{Int, Vector{Int}},
        d::AbstractDict{<:QSym},
    )
    slot = _resolve_slot(op, sites)
    haskey(d, op) && return d[op]
    return QuantumOpticsBase.LazyTensor(b, [slot], (to_numeric(op, b.bases[slot]),))
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
        )
    )
end

to_numeric(x::Number, b::QuantumOpticsBase.Basis, ::AbstractDict{<:QSym} = _NO_SUBS) = QuantumOpticsBase.sparse(_to_complex(x) * _lazy_one(b))
to_numeric(op::QField, state::QuantumState, d::AbstractDict{<:QSym} = _NO_SUBS) = to_numeric(op, QuantumOpticsBase.basis(state), d)
to_numeric(x::Number, state::QuantumState) = to_numeric(x, QuantumOpticsBase.basis(state))

_lazy_one(b::QuantumOpticsBase.Basis) = one(b)
_lazy_one(b::QuantumOpticsBase.CompositeBasis) = QuantumOpticsBase.LazyTensor(b, collect(1:length(b.bases)), Tuple(one(bi) for bi in b.bases))

function _reduce_const(n::Num)::ComplexF64
    x = Symbolics.value(n)
    try
        return _fold_const(x)
    catch err
        err isa ArgumentError || rethrow()
        isempty(Symbolics.get_variables(n)) || rethrow()
        return _compile_const(n)
    end
end

_compile_const(n::Num)::ComplexF64 = ComplexF64(symbolic_to_float(n))

function _fold_const(x)::ComplexF64
    x isa Number && return x
    if x isa SymbolicUtils.BasicSymbolic
        if SymbolicUtils.iscall(x)
            op = SymbolicUtils.operation(x)
            args = SymbolicUtils.arguments(x)
            if op === (+)
                acc = zero(ComplexF64)
                for a in args
                    acc += _fold_const(a)
                end
                return acc
            elseif op === (*)
                acc = one(ComplexF64)
                for a in args
                    acc *= _fold_const(a)
                end
                return acc
            elseif op === conj
                return conj(_fold_const(first(args)))
            elseif op === real
                return real(_fold_const(first(args)))
            elseif op === imag
                return imag(_fold_const(first(args)))
            elseif op === (/)
                return _fold_const(first(args)) / _fold_const(last(args))
            elseif op === (^)
                return _fold_const(first(args))^_fold_const(last(args))
            elseif op === (-)
                return length(args) == 1 ? -_fold_const(first(args)) :
                    _fold_const(first(args)) - _fold_const(last(args))
            end
        elseif SymbolicUtils.isconst(x)
            return x.val
        end
    end
    throw(ArgumentError("cannot reduce symbolic expression $x to a concrete number"))
end

_to_complex(c::Coeff) = _is_native(c) ? c.z : _to_complex(to_num(c))

# One method (union-split budget) routing every input through `convert ∘ Complex`,
# the only pattern that infers to ComplexF64 from `Any` after `Symbolics.value`.
function _to_complex(x)
    x isa ComplexF64   && return x
    x isa Complex{Num} && return _reduce_const(real(x)) + im * _reduce_const(imag(x))
    x isa Complex      && return convert(ComplexF64, x)
    x isa Num          && return _reduce_const(x)
    if x isa SymbolicUtils.BasicSymbolic
        SymbolicUtils.isconst(x) || throw(ArgumentError("cannot reduce symbolic expression $x to a concrete number"))
        return convert(ComplexF64, Complex(x.val, false))
    end
    return convert(ComplexF64, Complex(x, false))
end

"""
    numeric_average(op, state[, d::AbstractDict{<:QSym}]) -> ComplexF64
    numeric_average(op, states::AbstractVector[, d::AbstractDict{<:QSym}]) -> Vector

Expectation value ``\\langle \\psi | \\hat{O} | \\psi \\rangle`` of a symbolic
operator expression. Averaged `BasicSymbolic` expressions (from
[`average`](@ref)) are unwrapped automatically. `d` substitutes symbolic
operators with custom numeric operators before evaluation.

# Examples

```jldoctest
julia> using QuantumOpticsBase

julia> h = FockSpace(:f);

julia> @qnumbers a::Destroy(h);

julia> b = FockBasis(5); ψ = fockstate(b, 2);

julia> real(numeric_average(a' * a, ψ)) ≈ 2
true
```

See also [`to_numeric`](@ref), [`average`](@ref).
"""
function numeric_average end

numeric_average(op, state::QuantumState, d::AbstractDict{<:QSym} = _NO_SUBS) = _numeric_average(op, state, d)
function numeric_average(op, states::AbstractVector, d::AbstractDict{<:QSym} = _NO_SUBS)
    isempty(states) && throw(ArgumentError("numeric_average: states vector is empty"))
    return _numeric_average_vec(op, states, d)
end

_numeric_average(op::QField, state::QuantumState, d::AbstractDict{<:QSym}) = ComplexF64(QuantumOpticsBase.expect(_to_numeric_lazy(op, QuantumOpticsBase.basis(state), d), state))
_numeric_average(x::Number, ::QuantumState, ::AbstractDict{<:QSym}) = _to_complex(x)
_numeric_average(x::Num, state::QuantumState, d::AbstractDict{<:QSym}) = _numeric_average(SymbolicUtils.unwrap(x), state, d)

# `arguments(::BasicSymbolic)` is statically `Any`; the `_to_complex` wraps
# reseal each recursive result to `ComplexF64`.
function _numeric_average(avg::SymbolicUtils.BasicSymbolic, state::QuantumState, d::AbstractDict{<:QSym})
    if SymbolicUtils.hasmetadata(avg, AverageOperator) ||
            (SymbolicUtils.iscall(avg) && SymbolicUtils.operation(avg) isa AvgFunc) ||
            is_indexed_sum(avg)
        return _to_complex(_numeric_average(undo_average(avg), state, d))
    end
    SymbolicUtils.isconst(avg) && return _to_complex(_numeric_average(avg.val, state, d))
    SymbolicUtils.iscall(avg) || throw(ArgumentError("numeric_average not implemented for $(typeof(avg))"))
    f, args = SymbolicUtils.operation(avg), SymbolicUtils.arguments(avg)
    f === (^) && return _to_complex(_numeric_average(args[1], state, d))^_to_complex(args[2])
    f === (+) && return _sum_avg(args, state, d)
    f === (*) && return _prod_avg(args, state, d)
    throw(ArgumentError("numeric_average not implemented for operation $f"))
end

function _sum_avg(args, state::QuantumState, d::AbstractDict{<:QSym})
    r = ComplexF64(0)
    for a in args
        r += _to_complex(_numeric_average(a, state, d))
    end
    return r
end
function _prod_avg(args, state::QuantumState, d::AbstractDict{<:QSym})
    r = ComplexF64(1)
    for a in args
        r *= _to_complex(_numeric_average(a, state, d))
    end
    return r
end

_numeric_average_vec(op::QField, states::AbstractVector, d::AbstractDict{<:QSym}) = QuantumOpticsBase.expect(QuantumOpticsBase.sparse(to_numeric(op, first(states), d)), states)
_numeric_average_vec(op, states::AbstractVector, d::AbstractDict{<:QSym}) = [_numeric_average(op, ψ, d) for ψ in states]

"""
    expect(op::QField, state[, d::AbstractDict{<:QSym}])

Alias for [`numeric_average`](@ref) on operator expressions. For symbolic
scalar expressions such as `average(op)`, call [`numeric_average`](@ref)
directly.
"""
expect(op::QField, state, d::AbstractDict{<:QSym} = _NO_SUBS) = numeric_average(op, state, d)
