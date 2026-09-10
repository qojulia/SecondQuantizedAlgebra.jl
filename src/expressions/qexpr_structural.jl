# Structural operations for the cold-path formal expression layer.

function qexpr_rebuild(kind::QExprKind, coeff::CNum, args::Vector{QExprArg})
    if kind == QEXPR_ADD
        return qexpr_scale(qexpr_sum(args), coeff)
    elseif kind == QEXPR_MUL
        return qexpr_product(args, coeff)
    end
    length(args) == 1 || error("formal function node must have exactly one argument")
    return qexpr_scale(qexpr_call(kind, only(args)), coeff)
end

function substitute_qexpr_split(
        q::QExpr,
        op_rules::AbstractDict,
        scalar_rules::AbstractDict,
        phase_rules::Vector{Pair{Any, Any}},
    )
    new_coeff = substitute_cnum(q.coeff, scalar_rules, phase_rules)
    iszero_cnum(new_coeff) && return qexpr_zero()

    args = QExprArg[]
    sizehint!(args, length(q.args))
    for arg in q.args
        if arg isa QAdd
            push!(args, substitute_split(arg, op_rules, scalar_rules))
        else
            push!(args, substitute_qexpr_split(arg, op_rules, scalar_rules, phase_rules))
        end
    end

    if q.kind == QEXPR_EXPIM
        argument = only(args)
        qexpr_provably_hermitian(argument) || throw(
            ArgumentError(
                "substitution made the argument of `expim(::QField)` non-Hermitian",
            ),
        )
    end
    return qexpr_rebuild(q.kind, new_coeff, args)
end

function SymbolicUtils.substitute(q::QExpr, rules::AbstractDict; replace_adjoint = true)
    iszero(q) && return q
    op_rules, scalar_rules = split_substitution_rules(rules, replace_adjoint)
    phase_rules = nonreal_phase_substitutions(scalar_rules)
    return substitute_qexpr_split(q, op_rules, scalar_rules, phase_rules)
end

function change_index(q::QExpr, from::Index, to::Index)
    coeff = change_index(q.coeff, from, to)
    args = QExprArg[]
    sizehint!(args, length(q.args))
    for arg in q.args
        push!(args, change_index(arg, from, to))
    end
    return qexpr_rebuild(q.kind, coeff, args)
end

function change_index(q::QExpr, pairs::AbstractDict{Index, Index})
    isempty(pairs) && return q
    coeff = change_index(q.coeff, pairs)
    args = QExprArg[]
    sizehint!(args, length(q.args))
    for arg in q.args
        push!(args, change_index(arg, pairs))
    end
    return qexpr_rebuild(q.kind, coeff, args)
end

function get_indices(q::QExpr)
    indices = Index[]
    for arg in q.args
        for idx in get_indices(arg)
            idx in indices || push!(indices, idx)
        end
    end
    sort!(indices)
    return indices
end

function get_operators(q::QExpr)
    operators = Op[]
    for arg in q.args
        append!(operators, get_operators(arg))
    end
    unique!(operators)
    sort!(operators; by = order_key)
    return operators
end

function Symbolics.get_variables(q::QExpr; kwargs...)::Vector{SymbolicUtils.BasicSymbolic}
    vars = SymbolicUtils.BasicSymbolic[]
    append!(vars, Symbolics.get_variables(real(q.coeff); kwargs...))
    append!(vars, Symbolics.get_variables(imag(q.coeff); kwargs...))
    for arg in q.args
        append!(vars, Symbolics.get_variables(arg; kwargs...))
    end
    unique!(vars)
    sort!(vars; by = string)
    return vars
end

function Symbolics.get_variables(
        q::QExpr,
        varlist;
        is_atomic = SymbolicUtils.default_is_atomic,
        kwargs...,
    )::Vector{SymbolicUtils.BasicSymbolic}
    vars = SymbolicUtils.BasicSymbolic[]
    append!(
        vars,
        Symbolics.get_variables(real(q.coeff), varlist; is_atomic = is_atomic, kwargs...),
    )
    append!(
        vars,
        Symbolics.get_variables(imag(q.coeff), varlist; is_atomic = is_atomic, kwargs...),
    )
    for arg in q.args
        append!(
            vars,
            Symbolics.get_variables(arg, varlist; is_atomic = is_atomic, kwargs...),
        )
    end
    unique!(vars)
    sort!(vars; by = string)
    return vars
end

function Symbolics.get_variables!(buffer, q::QExpr; kwargs...)
    Symbolics.get_variables!(buffer, real(q.coeff); kwargs...)
    Symbolics.get_variables!(buffer, imag(q.coeff); kwargs...)
    for arg in q.args
        Symbolics.get_variables!(buffer, arg; kwargs...)
    end
    return buffer
end

function Symbolics.get_variables!(
        buffer,
        q::QExpr,
        varlist;
        is_atomic = SymbolicUtils.default_is_atomic,
        kwargs...,
    )
    Symbolics.get_variables!(
        buffer, real(q.coeff), varlist; is_atomic = is_atomic, kwargs...
    )
    Symbolics.get_variables!(
        buffer, imag(q.coeff), varlist; is_atomic = is_atomic, kwargs...
    )
    for arg in q.args
        Symbolics.get_variables!(buffer, arg, varlist; is_atomic = is_atomic, kwargs...)
    end
    return buffer
end

function acts_on(q::QExpr)
    spaces = Int[]
    for arg in q.args
        append!(spaces, acts_on(arg))
    end
    unique!(sort!(spaces))
    return spaces
end

function normal_order(q::QExpr)
    args = QExprArg[]
    sizehint!(args, length(q.args))
    for arg in q.args
        push!(args, normal_order(arg))
    end
    return qexpr_rebuild(q.kind, q.coeff, args)
end

function SymbolicUtils.simplify(q::QExpr; kwargs...)
    coeff = simplify_prefactor(q.coeff; kwargs...)
    args = QExprArg[]
    sizehint!(args, length(q.args))
    for arg in q.args
        push!(args, SymbolicUtils.simplify(arg; kwargs...))
    end
    return qexpr_rebuild(q.kind, coeff, args)
end

function Symbolics.expand(q::QExpr; kwargs...)
    coeff = expand_prefactor(q.coeff; kwargs...)
    args = QExprArg[]
    sizehint!(args, length(q.args))
    for arg in q.args
        push!(args, Symbolics.expand(arg; kwargs...))
    end
    return qexpr_rebuild(q.kind, coeff, args)
end

function expand_completeness(q::QExpr)
    args = QExprArg[]
    sizehint!(args, length(q.args))
    for arg in q.args
        push!(args, expand_completeness(arg))
    end
    return qexpr_rebuild(q.kind, q.coeff, args)
end

function assume_distinct_index(q::QExpr, pairs::Vector{Tuple{Index, Index}})
    args = QExprArg[]
    sizehint!(args, length(q.args))
    for arg in q.args
        push!(args, assume_distinct_index(arg, pairs))
    end
    return qexpr_rebuild(q.kind, q.coeff, args)
end

commutator(a::QExpr, b::QExpr) = a * b - b * a
commutator(a::QExpr, b::Union{QSym, QAdd}) = a * b - b * a
commutator(a::Union{QSym, QAdd}, b::QExpr) = a * b - b * a
commutator(::QExpr, ::Number) = qexpr_zero()
commutator(::Number, ::QExpr) = qexpr_zero()
