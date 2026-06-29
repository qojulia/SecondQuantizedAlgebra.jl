# Backend-neutral core: walk a canonical `QAdd` into a backend-neutral term list and hand
# it to the backend's lazy assembler. The core never does operator arithmetic and never
# names a concrete numeric type; everything numeric goes through the hooks in `backend.jl`.

# `Union{}` value keeps the haskey-fallback branch concrete-typed (a real op_subs dict
# widens the leaf to `Any`, which only the kwargs/dict paths use and which are not on the
# `@inferred` contract).
const _NO_SUBS = Dict{QSym, Union{}}()
const _EMPTY_SITES = Dict{Int, Vector{Int}}()

"""
    NumericContext(backend, basis, op_subs, scalar_subs, sites)

Immutable, all-concrete-fields bundle threaded through the numeric core: the backend
singleton, the opaque per-backend `basis`, the operator-substitution dict `op_subs`, the
indexed `scalar_subs`, and the indexed `sites` map (empty for the non-indexed path). All
fields are concrete via the type parameters, so the context never widens inference.
"""
struct NumericContext{BE <: NumericBackend, B, SUBS <: AbstractDict, SS <: AbstractDict}
    backend::BE
    basis::B
    op_subs::SUBS
    scalar_subs::SS
    sites::Dict{Int, Vector{Int}}
end

NumericContext(be::NumericBackend, basis, op_subs::AbstractDict = _NO_SUBS) =
    NumericContext(be, basis, op_subs, _NO_SCALAR_SUBS, _EMPTY_SITES)

# --- backend resolution for QuantumOptics states (QI types, so it lives here) -----------
# A `Basis` is inherently QuantumOptics, so the Basis-form `to_numeric` hardcodes
# `QuantumOpticsBackend()` rather than going through `_backend_of` (keeping `_backend_of`
# strictly a *state* query: a `Basis` is never a state, so it must not appear as one). The
# QuantumToolbox extension adds the `QuantumObject` state methods.
_backend_of(::Union{StateVector, AbstractOperator}) = QuantumOpticsBackend()
_basis_of(s::Union{StateVector, AbstractOperator}) = basis(s)

# --- leaf / term construction ----------------------------------------------------------

# The numeric leaf for one operator: a custom override from `op_subs`, else the backend
# operator embedded at its slot. With `op_subs = _NO_SUBS` the override branch has value
# type `Union{}` and contributes nothing, so the inferred type is exactly the embed type.
function _numeric_leaf(op::Op, ctx::NumericContext)
    haskey(ctx.op_subs, op) && return ctx.op_subs[op]
    slot = Int(op.space_index)
    sub = numeric_subbasis(ctx.backend, ctx.basis, slot)
    leaf = numeric_operator(ctx.backend, op, sub)
    return numeric_embed(ctx.backend, ctx.basis, slot, leaf)
end

_term_factors(ops::Vector{Op}, ctx::NumericContext) = [_numeric_leaf(op, ctx) for op in ops]

# --- static assembly -------------------------------------------------------------------

# Build the (ComplexF64, factors) term list and assemble. Coefficients must already be
# concrete (the kwargs path substitutes `parameter` before calling this); `_to_complex`
# throws `ArgumentError` for a still-symbolic coefficient.
function _to_numeric_static(q::QAdd, ctx::NumericContext)
    terms = [(_to_complex(c), _term_factors(term.ops, ctx)) for (term, c) in q.arguments]
    return numeric_assemble(ctx.backend, ctx.basis, terms)
end

# --- time-dependent assembly -----------------------------------------------------------

# Each coefficient becomes a `ComplexF64` (constant term) or a `t -> ComplexF64` closure
# (genuinely time-dependent), reusing the kept coefficient machinery. Factors are stored
# behind `Vector{Any}` because the TD path is not on the `@inferred` contract and the
# backend assembler only needs to iterate them.
const _TDTerm = Tuple{Union{ComplexF64, Function}, Vector{Any}}

function _td_coeff(c::Complex{Num}, basevars, valuefuncs)
    if _coeff_is_const(c)
        return _const_coeff(c)
    end
    _check_time_variables(c, basevars)
    pref = _compile_coeff(c, basevars...)
    return let pref = pref, valuefuncs = valuefuncs
        t -> pref(map(f -> f(t), valuefuncs)...)
    end
end

function _to_numeric_td(q::QAdd, ctx::NumericContext, time_parameter)
    basevars, valuefuncs = _time_basis(time_parameter)
    td_terms = _TDTerm[]
    for (term, c_) in q.arguments
        coef = _td_coeff(to_num(c_), basevars, valuefuncs)
        push!(td_terms, (coef, collect(Any, _term_factors(term.ops, ctx))))
    end
    isempty(td_terms) && push!(td_terms, (0.0im, Any[]))
    return numeric_assemble_td(ctx.backend, ctx.basis, td_terms)
end

# --- keyword translation boundary ------------------------------------------------------

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

function _to_numeric_kw(be::NumericBackend, op, b; parameter, time_parameter, operators, adjoint_ops)
    param = _expand_parameter(parameter)
    tp = _normalize_time_parameter(time_parameter)
    ops = _numeric_operator_dict(operators, adjoint_ops)
    ctx = NumericContext(be, b, ops)
    return _to_numeric_translated(op, ctx, param, tp)
end

# QSym: wrap as a one-term QAdd and reuse the QAdd path.
_to_numeric_translated(op::QSym, ctx::NumericContext, parameter, time_parameter) =
    _to_numeric_translated(_single_qadd(_CNUM_ONE, Op[op]), ctx, parameter, time_parameter)

function _to_numeric_translated(op::QAdd, ctx::NumericContext, parameter, time_parameter)
    op_ = substitute(op, parameter)
    isempty(time_parameter) && return _to_numeric_static(op_, ctx)
    return _to_numeric_td(op_, ctx, time_parameter)
end

# Bare scalar: static -> scaled identity; time-dependent -> native TD over identity.
function _to_numeric_translated(arg, ctx::NumericContext, parameter, time_parameter)
    arg_sub = substitute(arg, parameter)
    c = _as_cnum(arg_sub)
    if isempty(time_parameter)
        _coeff_is_const(c) || throw(
            ArgumentError(
                "cannot translate symbolic scalar `$arg_sub` without a value: supply it via " *
                    "`parameter` or `time_parameter`.",
            ),
        )
        return _const_coeff(c) * numeric_identity(ctx.backend, ctx.basis)
    end
    basevars, valuefuncs = _time_basis(time_parameter)
    coef = _td_coeff(c, basevars, valuefuncs)
    return numeric_assemble_td(ctx.backend, ctx.basis, _TDTerm[(coef, Any[])])
end
