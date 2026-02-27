function average end

"""
    AvgSym <: CNumber

Symbolic number representing the average over an operator.
See also: [`average`](@ref)
"""
struct AvgSym <: CNumber end

const Average = SQABasicSymbolic

function is_average(x)
    x isa SQABasicSymbolic || return false
    is_symtype(x, AvgSym) && return true
    if SymbolicUtils.iscall(x)
        f = SymbolicUtils.operation(x)
        if f === (+) || f === (-)
            return all(arg -> is_average(unwrap_const(arg)), SymbolicUtils.arguments(x))
        end
    end
    return false
end

const sym_average = begin # Symbolic function for averages
    T = SymbolicUtils.FnType{Tuple{QNumber},AvgSym,Nothing}
    SymbolicUtils.Sym{SQA_VARTYPE}(:avg; type=T)
end

# Type promotion -- average(::QNumber)::Number
SymbolicUtils.promote_symtype(::typeof(sym_average), ::Type{<:QNumber}) = AvgSym
SymbolicUtils.promote_symtype(::typeof(*), ::Type{AvgSym}, ::Type{AvgSym}) = CNumber
SymbolicUtils.promote_symtype(::typeof(*), ::Type{AvgSym}, ::Type{<:CNumber}) = CNumber
SymbolicUtils.promote_symtype(::typeof(*), ::Type{<:CNumber}, ::Type{AvgSym}) = CNumber
SymbolicUtils.promote_symtype(::typeof(+), ::Type{AvgSym}, ::Type{AvgSym}) = CNumber
SymbolicUtils.promote_symtype(::typeof(+), ::Type{AvgSym}, ::Type{<:CNumber}) = CNumber
SymbolicUtils.promote_symtype(::typeof(+), ::Type{<:CNumber}, ::Type{AvgSym}) = CNumber
SymbolicUtils.promote_symtype(::typeof(-), ::Type{AvgSym}, ::Type{AvgSym}) = CNumber
SymbolicUtils.promote_symtype(::typeof(-), ::Type{AvgSym}, ::Type{<:CNumber}) = CNumber
SymbolicUtils.promote_symtype(::typeof(-), ::Type{<:CNumber}, ::Type{AvgSym}) = CNumber

# Direct construction of average symbolic expression
function _average(operator)
    return SymbolicUtils.Term{SQA_VARTYPE}(sym_average, [operator]; type=AvgSym)
end

# ensure that BasicSymbolic{<:AvgSym} are only single averages
# function *(a::Average, b::Average)
#     if isequal(a, b)
#         return SymbolicUtils.Mul(CNumber, 1, Dict(a=>2))
#     end
#     return SymbolicUtils.Mul(CNumber, 1, Dict(a=>1, b=>1))
# end
# function +(a::Average, b::Average)
#     if isequal(a, b)
#         return SymbolicUtils.Add(CNumber, 0, Dict(a=>2))
#     end
#     return SymbolicUtils.Add(CNumber, 0, Dict(a=>1, b=>1))
# end
# ∨ https://github.com/qojulia/SecondQuantizedAlgebra.jl/issues/28
function SymbolicUtils.Add(::Type{AvgSym}, coeff, dict; kw...)
    SymbolicUtils.Add(SQA_VARTYPE, coeff, dict; type=CNumber, kw...)
end
function SymbolicUtils.Mul(::Type{AvgSym}, coeff, dict; kw...)
    SymbolicUtils.Mul(SQA_VARTYPE, coeff, dict; type=CNumber, kw...)
end

function acts_on(s::SymbolicUtils.BasicSymbolic)
    if SymbolicUtils.iscall(s)
        f = SymbolicUtils.operation(s)
        if f === sym_average
            return acts_on(unwrap_const(SymbolicUtils.arguments(s)[1]))
        else
            aon = []
            for arg in SymbolicUtils.arguments(s)
                append!(aon, acts_on(unwrap_const(arg)))
            end
            unique!(aon)
            sort!(aon)
            return aon
        end
    else
        return Int[]
    end
end

"""
    average(::QNumber)

Compute the average of an operator.
"""
average(op::QSym) = _average(op)
function average(op::QTerm)
    f = SymbolicUtils.operation(op)
    if f===(+) || f===(-) # linearity
        args = map(average, SymbolicUtils.arguments(op))
        return f(args...)
    elseif f === (*)
        # Move constants out of average
        c = op.arg_c
        op_ = QMul(1, op.args_nc)
        return c*_average(op_)
    else
        error("Unknown function $f")
    end
end

average(x::SNuN) = x

function undo_average(t)
    if SymbolicUtils.iscall(t)
        f = SymbolicUtils.operation(t)
        if isequal(f, sym_average) # "===" results in false sometimes in Symbolics version > 5
            return unwrap_const(SymbolicUtils.arguments(t)[1])
        else
            args = map(arg -> undo_average(unwrap_const(arg)), SymbolicUtils.arguments(t))
            return f(args...)
        end
    else
        return t
    end
end

function undo_average(eq::Symbolics.Equation)
    lhs = undo_average(eq.lhs)
    rhs = undo_average(eq.rhs)
    return Symbolics.Equation(lhs, rhs)
end
