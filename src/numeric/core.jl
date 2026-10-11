const NO_SUBS = Dict{QSym, Union{}}()
const EMPTY_SITES = Dict{Int, Vector{Int}}()

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

NumericContext(be::NumericBackend, basis, op_subs::AbstractDict = NO_SUBS) =
    NumericContext(be, basis, op_subs, NO_SCALAR_SUBS, EMPTY_SITES)

numeric_backend(::Union{StateVector, AbstractOperator}) = QuantumOpticsBackend()

function numeric_leaf(op::Op, ctx::NumericContext)
    haskey(ctx.op_subs, op) && return ctx.op_subs[op]
    slot = Int(op.space_index)
    sub = numeric_subbasis(ctx.backend, ctx.basis, slot)
    leaf = numeric_operator(ctx.backend, op, sub)
    return numeric_embed(ctx.backend, ctx.basis, slot, leaf)
end

term_factors(ops::Vector{Op}, ctx::NumericContext) = [numeric_leaf(op, ctx) for op in ops]

function to_numeric_static(q::QAdd, ctx::NumericContext)
    terms = [(to_complex(c), term_factors(term.ops, ctx)) for (term, c) in q.arguments]
    return numeric_assemble(ctx.backend, ctx.basis, terms)
end

to_numeric_lazy(op::Op, ctx::NumericContext) = numeric_leaf(op, ctx)
to_numeric_lazy(q::QAdd, ctx::NumericContext) = to_numeric_static(q, ctx)

const TDTerm = Tuple{Union{ComplexF64, Function, PControlCoeff}, Vector{Any}}

function td_coeff(c::Complex{Num}, basevars, valuefuncs, p_aware::Bool)
    if coeff_is_const(c)
        return const_coeff(c)
    end
    check_time_variables(c, basevars)
    pref = compile_coeff(c, basevars...)
    p_aware && return PControlCoeff(
        let pref = pref, valuefuncs = valuefuncs
            (p, t) -> pref(map(f -> f(p, t), valuefuncs)...)
        end
    )
    return let pref = pref, valuefuncs = valuefuncs
        t -> pref(map(f -> f(t), valuefuncs)...)
    end
end

function to_numeric_td(q::QAdd, ctx::NumericContext, time_parameter)
    basevars, valuefuncs, p_aware = time_basis(time_parameter)
    td_terms = TDTerm[]
    for (term, c_) in q.arguments
        coef = td_coeff(to_num(c_), basevars, valuefuncs, p_aware)
        push!(td_terms, (coef, collect(Any, term_factors(term.ops, ctx))))
    end
    isempty(td_terms) && push!(td_terms, (0.0im, Any[]))
    return numeric_assemble_td(ctx.backend, ctx.basis, td_terms)
end

function numeric_operator_dict(operators, adjoint_ops::Bool)
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

function to_numeric_kw(be::NumericBackend, op, b; parameter, time_parameter, operators, adjoint_ops, op_type)
    param = expand_parameter(parameter)
    tp = normalize_time_parameter(time_parameter)
    ops = numeric_operator_dict(operators, adjoint_ops)
    ctx = NumericContext(be, b, ops)
    return to_numeric_translated(op, ctx, param, tp, op_type)
end

to_numeric_translated(op::QSym, ctx::NumericContext, parameter, time_parameter, op_type) =
    to_numeric_translated(single_qadd(CNUM_ONE, Op[op]), ctx, parameter, time_parameter, op_type)

function to_numeric_translated(op::QAdd, ctx::NumericContext, parameter, time_parameter, op_type)
    op_ = substitute(op, parameter)
    isempty(time_parameter) && return numeric_materialize(ctx.backend, to_numeric_static(op_, ctx), op_type)
    check_td_op_type(op_type)
    return to_numeric_td(op_, ctx, time_parameter)
end

function to_numeric_translated(arg, ctx::NumericContext, parameter, time_parameter, op_type)
    arg_sub = substitute(arg, parameter)
    c = as_cnum(arg_sub)
    if isempty(time_parameter)
        coeff_is_const(c) || throw(
            ArgumentError(
                "cannot translate symbolic scalar `$arg_sub` without a value: supply it via " *
                    "`parameter` or `time_parameter`.",
            ),
        )
        return numeric_materialize(ctx.backend, const_coeff(c) * numeric_identity(ctx.backend, ctx.basis), op_type)
    end
    check_td_op_type(op_type)
    basevars, valuefuncs, p_aware = time_basis(time_parameter)
    coef = td_coeff(c, basevars, valuefuncs, p_aware)
    return numeric_assemble_td(ctx.backend, ctx.basis, TDTerm[(coef, Any[])])
end

function check_td_op_type(op_type)
    (op_type isa EagerOperator || op_type === identity) && return op_type
    throw(
        ArgumentError(
            "op_type is not supported when time_parameter is non-empty; native " *
                "time-dependent operators cannot be materialized once at conversion time",
        ),
    )
end
