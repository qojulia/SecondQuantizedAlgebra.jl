function map_args(f::F, q::QExpr, coeff::CNum = q.coeff) where {F}
    q.kind == QEXPR_ADD || q.kind == QEXPR_MUL ||
        return qexpr_scale(qexpr_call(q.kind, f(only(q.args))), coeff)
    return qexpr_rebuild(q.kind, coeff, QExprArg[f(arg) for arg in q.args])
end

function substitute_with(q::QExpr, r::SubstitutionRules)
    (iszero(q) || isempty(r)) && return q
    new_coeff = substitute_cnum(q.coeff, r)
    iszero(new_coeff) && return qexpr_zero()
    return map_args(q, new_coeff) do arg
        new_arg = substitute_with(arg, r)
        q.kind == QEXPR_EXPIM && !qexpr_provably_hermitian(new_arg) && throw(
            ArgumentError("substitution made the argument of `expim(::QField)` non-Hermitian"),
        )
        new_arg
    end
end

SymbolicUtils.substitute(q::QExpr, rules::AbstractDict; replace_adjoint::Bool = true) =
    substitute_with(q, SubstitutionRules(rules, replace_adjoint))

function change_index(q::QExpr, m::IndexMap)
    renames_nothing(m) && return q
    return map_args(arg -> change_index(arg, m), q, change_index(q.coeff, m))
end

get_indices(q::QExpr) = sort!(unique!(mapreduce(get_indices, vcat, q.args; init = Index[])))

get_operators(q::QExpr) =
    sort!(unique!(mapreduce(get_operators, vcat, q.args; init = Op[])); by = order_key)

function acts_on(q::QExpr)
    spaces = Int[]
    for arg in q.args
        append!(spaces, acts_on(arg))
    end
    return unique!(sort!(spaces))
end

function Symbolics.get_variables(q::QExpr, varlist...; kwargs...)
    vars = Symbolics.get_variables!(SymbolicUtils.BasicSymbolic[], q, varlist...; kwargs...)
    unique!(vars)
    sort!(vars; by = string)
    return vars
end

function Symbolics.get_variables!(buffer, q::QExpr, varlist...; kwargs...)
    Symbolics.get_variables!(buffer, real(q.coeff), varlist...; kwargs...)
    Symbolics.get_variables!(buffer, imag(q.coeff), varlist...; kwargs...)
    for arg in q.args
        Symbolics.get_variables!(buffer, arg, varlist...; kwargs...)
    end
    return buffer
end

normal_order(q::QExpr) = map_args(normal_order, q)

SymbolicUtils.simplify(q::QExpr; kwargs...) =
    map_args(arg -> SymbolicUtils.simplify(arg; kwargs...), q, simplify_prefactor(q.coeff; kwargs...))

Symbolics.expand(q::QExpr; kwargs...) =
    map_args(arg -> Symbolics.expand(arg; kwargs...), q, expand_prefactor(q.coeff; kwargs...))

expand_completeness(q::QExpr) = map_args(expand_completeness, q)

assume_distinct_index(q::QExpr, pairs::Vector{Tuple{Index, Index}}) =
    map_args(arg -> assume_distinct_index(arg, pairs), q)

commutator(a::QExpr, b::Union{QSym, QAdd, QExpr}) = a * b - b * a
commutator(a::Union{QSym, QAdd}, b::QExpr) = a * b - b * a
commutator(::QExpr, ::Number) = qexpr_zero()
commutator(::Number, ::QExpr) = qexpr_zero()
