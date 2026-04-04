"""
    QAdd{T<:Number} <: QTerm

Sum of [`QMul{T}`](@ref) terms.

Fields:
- `arguments::Vector{QMul{T}}` — terms of the sum
"""
struct QAdd{T<:Number} <: QTerm
    arguments::Vector{QMul{T}}
    function QAdd(args::Vector{QMul{T}}) where {T<:Number}
        return new{T}(args)
    end
end

Base.length(a::QAdd) = length(a.arguments)

# Equality and hashing
function Base.isequal(a::QAdd, b::QAdd)
    length(a.arguments) == length(b.arguments) || return false
    for (x, y) in zip(a.arguments, b.arguments)
        isequal(x, y) || return false
    end
    return true
end
Base.hash(q::QAdd, h::UInt) = hash(:QAdd, hash(q.arguments, h))

# Adjoint
function Base.adjoint(q::QAdd{T}) where {T}
    return QAdd(QMul{T}[adjoint(t) for t in q.arguments])
end

# Promote
Base.promote_rule(::Type{QAdd{S}}, ::Type{QAdd{T}}) where {S,T} = QAdd{promote_type(S, T)}
function Base.convert(::Type{QAdd{T}}, x::QAdd{S}) where {T<:Number,S<:Number}
    return QAdd(QMul{T}[convert(QMul{T}, m) for m in x.arguments])
end

# Helpers: wrap QSym/scalar as QMul
_to_qmul(a::QSym, ::Type{T}) where {T} = QMul(one(T), QSym[a])
_to_qmul(a::QMul{S}, ::Type{T}) where {S,T} = convert(QMul{promote_type(S, T)}, a)
_scalar_qmul(x::T) where {T<:Number} = QMul(x, QSym[])

## Addition — always returns QAdd

# QSym + QSym
function Base.:+(a::QSym, b::QSym)
    return QAdd(QMul{Int}[_to_qmul(a, Int), _to_qmul(b, Int)])
end

# QMul + QMul
function Base.:+(a::QMul{S}, b::QMul{T}) where {S,T}
    TT = promote_type(S, T)
    return QAdd(QMul{TT}[convert(QMul{TT}, a), convert(QMul{TT}, b)])
end

# QMul + QSym
function Base.:+(a::QMul{T}, b::QSym) where {T}
    return QAdd(QMul{T}[a, _to_qmul(b, T)])
end
Base.:+(a::QSym, b::QMul{T}) where {T} = b + a

# QAdd + QMul
function Base.:+(a::QAdd{S}, b::QMul{T}) where {S,T}
    TT = promote_type(S, T)
    args = QMul{TT}[convert(QMul{TT}, x) for x in a.arguments]
    push!(args, convert(QMul{TT}, b))
    return QAdd(args)
end
Base.:+(a::QMul, b::QAdd) = b + a

# QAdd + QSym
function Base.:+(a::QAdd{T}, b::QSym) where {T}
    return a + _to_qmul(b, T)
end
Base.:+(a::QSym, b::QAdd) = b + a

# QAdd + QAdd
function Base.:+(a::QAdd{S}, b::QAdd{T}) where {S,T}
    TT = promote_type(S, T)
    args = QMul{TT}[convert(QMul{TT}, x) for x in a.arguments]
    for x in b.arguments
        push!(args, convert(QMul{TT}, x))
    end
    return QAdd(args)
end

# QField + Number
function Base.:+(a::QSym, b::T) where {T<:Number}
    return QAdd(QMul{T}[_to_qmul(a, T), _scalar_qmul(b)])
end
Base.:+(a::Number, b::QSym) = b + a

function Base.:+(a::QMul{S}, b::T) where {S,T<:Number}
    TT = promote_type(S, T)
    return QAdd(QMul{TT}[convert(QMul{TT}, a), QMul(convert(TT, b), QSym[])])
end
Base.:+(a::Number, b::QMul) = b + a

function Base.:+(a::QAdd{S}, b::T) where {S,T<:Number}
    TT = promote_type(S, T)
    args = QMul{TT}[convert(QMul{TT}, x) for x in a.arguments]
    push!(args, QMul(convert(TT, b), QSym[]))
    return QAdd(args)
end
Base.:+(a::Number, b::QAdd) = b + a

# Subtraction
Base.:-(a::QAdd) = QAdd(QMul[-t for t in a.arguments])
Base.:-(a::QField, b::QField) = a + (-b)
Base.:-(a::QField, b::Number) = a + (-b)
Base.:-(a::Number, b::QField) = a + (-b)

## QAdd * ... (distributive)

# QAdd * Number
function Base.:*(a::QAdd{S}, b::T) where {S,T<:Number}
    TT = promote_type(S, T)
    return QAdd(QMul{TT}[convert(QMul{TT}, t) * b for t in a.arguments])
end
Base.:*(a::Number, b::QAdd) = b * a

# QAdd * QSym
function Base.:*(a::QAdd{T}, b::QSym) where {T}
    return QAdd(QMul{T}[t * b for t in a.arguments])
end
function Base.:*(a::QSym, b::QAdd{T}) where {T}
    return QAdd(QMul{T}[a * t for t in b.arguments])
end

# QAdd * QMul
function Base.:*(a::QAdd{S}, b::QMul{T}) where {S,T}
    TT = promote_type(S, T)
    return QAdd(QMul{TT}[t * b for t in a.arguments])
end
function Base.:*(a::QMul{S}, b::QAdd{T}) where {S,T}
    TT = promote_type(S, T)
    return QAdd(QMul{TT}[a * t for t in b.arguments])
end

# QAdd * QAdd
function Base.:*(a::QAdd{S}, b::QAdd{T}) where {S,T}
    TT = promote_type(S, T)
    args = QMul{TT}[]
    for ai in a.arguments, bi in b.arguments
        push!(args, ai * bi)
    end
    return QAdd(args)
end

# QAdd / Number
Base.:/(a::QAdd, b::Number) = a * (1 // b)
