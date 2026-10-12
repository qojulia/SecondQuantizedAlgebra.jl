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
    slot = resolve_slot(op, ctx.sites)
    sub = numeric_subbasis(ctx.backend, ctx.basis, slot)
    leaf = numeric_operator(ctx.backend, op, sub)
    return numeric_embed(ctx.backend, ctx.basis, slot, leaf)
end

function term_factors(ops::Vector{Op}, ctx::C) where {C <: NumericContext}
    factors = Vector{Base.promote_op(numeric_leaf, Op, C)}(undef, length(ops))
    for (k, op) in enumerate(ops)
        factors[k] = numeric_leaf(op, ctx)
    end
    return factors
end

const OpTerm = Tuple{ComplexF64, Vector{Op}}

op_terms(q::QAdd) = OpTerm[(to_complex(c), term.ops) for (term, c) in q.arguments]

function assemble_op_terms(terms::Vector{OpTerm}, ctx::NumericContext)
    numeric_terms = [(c, term_factors(ops, ctx)) for (c, ops) in terms]
    isempty(numeric_terms) && push!(numeric_terms, (0.0im, term_factors(Op[], ctx)))
    return numeric_assemble(ctx.backend, ctx.basis, numeric_terms)
end

to_numeric_static(q::QAdd, ctx::NumericContext) = assemble_op_terms(op_terms(q), ctx)

static_operator(terms::Vector{OpTerm}, ctx::NumericContext, op_type) =
    numeric_materialize(ctx.backend, assemble_op_terms(terms, ctx), op_type)
static_operator(q::QAdd, ctx::NumericContext, op_type) = static_operator(op_terms(q), ctx, op_type)

to_numeric_lazy(op::Op, ctx::NumericContext) = numeric_leaf(op, ctx)
to_numeric_lazy(q::QAdd, ctx::NumericContext) = to_numeric_static(q, ctx)

const TDCoeff = Union{ComplexF64, Function}

function td_coeff(c::Complex{Num}, basevars::Vector{Num}, @nospecialize(valuefuncs), p_aware::Bool)
    if coeff_is_const(c)
        return const_coeff(c)
    end
    check_time_variables(c, basevars)
    pref = compile_coeff(c, basevars...)
    p_aware && return let pref = pref, valuefuncs = valuefuncs
        (p, t) -> pref(map(f -> f(p, t), valuefuncs)...)
    end
    return let pref = pref, valuefuncs = valuefuncs
        t -> pref(map(f -> f(t), valuefuncs)...)
    end
end

function to_numeric_td(q::QAdd, ctx::NumericContext, time_parameter)
    basevars, valuefuncs, p_aware = time_basis(time_parameter)
    no_factors = term_factors(Op[], ctx)
    td_terms = Tuple{TDCoeff, typeof(no_factors)}[]
    for (term, c_) in q.arguments
        coef = td_coeff(to_num(c_), basevars, valuefuncs, p_aware)
        push!(td_terms, (coef, term_factors(term.ops, ctx)))
    end
    isempty(td_terms) && push!(td_terms, (0.0im, no_factors))
    return numeric_assemble_td(ctx.backend, ctx.basis, td_terms)
end

numeric_operator_dict(::NumericBackend, operators::Dict{QSym, Union{}}, ::Bool) = operators
function numeric_operator_dict(be::NumericBackend, operators::AbstractDict, adjoint_ops::Bool)
    explicit = Dict(operator_key(k) => v for (k, v) in operators)
    adjoints = Dict(Base.adjoint(k) => operator_adjoint(be, v) for (k, v) in explicit if adjoint_ops)
    out = Dict{Op, Union{valtype(explicit), valtype(adjoints)}}()
    return merge!(out, adjoints, explicit)
end

operator_adjoint(::NumericBackend, v) = Base.adjoint(v)

operator_key(k::Op) = k
operator_key(k) = throw(ArgumentError("operator substitution key `$k` is not a QSym"))

function to_numeric_kw(be::NumericBackend, op, b; parameter, time_parameter, operators, adjoint_ops, op_type)
    rules = expand_parameter(parameter)
    tp = normalize_time_parameter(time_parameter)
    ops = numeric_operator_dict(be, operators, adjoint_ops)
    ctx = NumericContext(be, b, ops)
    return to_numeric_translated(op, ctx, rules, tp, op_type)
end

apply_parameter(q::QAdd, (op_rules, scalar_rules)::Tuple{OpRules, ScalarRules}) =
    isempty(op_rules) && isempty(scalar_rules) ? q : substitute_split(q, op_rules, scalar_rules)

apply_parameter(x::Complex{Num}, (_, scalar_rules)::Tuple{OpRules, ScalarRules}) =
    isempty(scalar_rules) ? x : Complex{Num}(
        Num(SymbolicUtils.substitute(raw_rule(real(x)), scalar_rules)),
        Num(SymbolicUtils.substitute(raw_rule(imag(x)), scalar_rules)),
    )
apply_parameter(x::Number, (_, scalar_rules)::Tuple{OpRules, ScalarRules}) =
    isempty(scalar_rules) ? x : SymbolicUtils.substitute(raw_rule(x), scalar_rules)
apply_parameter(x::SymbolicUtils.BasicSymbolic, (_, scalar_rules)::Tuple{OpRules, ScalarRules}) =
    isempty(scalar_rules) ? x : SymbolicUtils.substitute(x, scalar_rules)

to_numeric_translated(op::QSym, ctx::NumericContext, parameter, time_parameter, op_type) =
    to_numeric_translated(single_qadd(CNUM_ONE, Op[op]), ctx, parameter, time_parameter, op_type)

function to_numeric_translated(op::QAdd, ctx::NumericContext, parameter, time_parameter, op_type)
    op_ = apply_parameter(op, parameter)
    isempty(time_parameter) && return static_operator(op_, ctx, op_type)
    check_td_op_type(op_type)
    return to_numeric_td(op_, ctx, time_parameter)
end

function to_numeric_translated(
        arg::Union{Number, SymbolicUtils.BasicSymbolic}, ctx::NumericContext, parameter, time_parameter, op_type,
    )
    arg_sub = apply_parameter(arg, parameter)
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
    return to_numeric_td(single_qadd(to_cnum(c), Op[]), ctx, time_parameter)
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
