## TermInterface / SymbolicUtils integration

# Head
TermInterface.head(::QField) = :call

# QSym — leaves (not callable)
SymbolicUtils.iscall(::QSym) = false
TermInterface.metadata(::QSym) = nothing

# QMul — products
SymbolicUtils.iscall(::QMul) = true
SymbolicUtils.iscall(::Type{<:QMul}) = true
SymbolicUtils.operation(::QMul) = (*)
SymbolicUtils.arguments(a::QMul) = Any[a.arg_c, a.args_nc...]
TermInterface.metadata(::QMul) = nothing

function TermInterface.maketerm(::Type{<:QMul}, ::typeof(*), args, metadata)
    args_c = filter(x -> !(x isa QField), args)
    args_nc = filter(x -> x isa QField, args)
    arg_c = isempty(args_c) ? 1 : *(args_c...)
    return QMul(arg_c, QSym[args_nc...])
end

# QAdd — sums
SymbolicUtils.iscall(::QAdd) = true
SymbolicUtils.iscall(::Type{<:QAdd}) = true
SymbolicUtils.operation(::QAdd) = (+)
SymbolicUtils.arguments(a::QAdd) = a.arguments
TermInterface.metadata(::QAdd) = nothing

function TermInterface.maketerm(::Type{<:QAdd}, ::typeof(+), args, metadata)
    # TermInterface may pass Vector{Any}; convert to concrete QMul vector
    muls = QMul[x isa QMul ? x : QMul(1, QSym[x]) for x in args]
    TT = promote_type((typeof(m.arg_c) for m in muls)...)
    return QAdd(QMul{TT}[convert(QMul{TT}, m) for m in muls])
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

# one / zero / isone / iszero for QField
Base.one(::T) where {T <: QField} = one(T)
Base.one(::Type{<:QField}) = 1
Base.isone(::QField) = false
Base.zero(::T) where {T <: QField} = zero(T)
Base.zero(::Type{<:QField}) = 0
