## TermInterface / SymbolicUtils integration

TermInterface.head(::QField) = :call

# QSym — leaves (not callable)
SymbolicUtils.iscall(::QSym) = false
TermInterface.metadata(::QSym) = nothing

# QAdd — sums (the only compound expression type)
SymbolicUtils.iscall(::QAdd) = true
SymbolicUtils.iscall(::Type{QAdd}) = true
SymbolicUtils.operation(::QAdd) = (+)
SymbolicUtils.arguments(a::QAdd) = sorted_arguments(a)
TermInterface.metadata(::QAdd) = nothing

function TermInterface.maketerm(::Type{QAdd}, ::typeof(+), args, metadata)
    result = args[1]
    for i in 2:length(args)
        result = result + args[i]
    end
    return result
end

# Type promotion
SymbolicUtils.symtype(x::T) where {T <: QField} = T

for f in SymbolicUtils.basic_diadic
    @eval SymbolicUtils.promote_symtype(::$(typeof(f)), Ts::Type{<:QField}...) = promote_type(Ts...)
    @eval SymbolicUtils.promote_symtype(::$(typeof(f)), T::Type{<:QField}, Ts...) = T
    @eval SymbolicUtils.promote_symtype(::$(typeof(f)), T::Type{<:QField}, S::Type{<:Number}) = T
    @eval SymbolicUtils.promote_symtype(::$(typeof(f)), T::Type{<:Number}, S::Type{<:QField}) = S
    @eval SymbolicUtils.promote_symtype(::$(typeof(f)), T::Type{<:QField}, S::Type{<:QField}) = promote_type(T, S)
end

Base.one(::T) where {T <: QField} = one(T)
Base.one(::Type{<:QField}) = 1
Base.isone(::QField) = false
Base.zero(::T) where {T <: QField} = zero(T)
Base.zero(::Type{<:QField}) = 0
