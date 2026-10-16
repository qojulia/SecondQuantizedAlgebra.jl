Base.:*(a::QSym, b::Coefficient) = single_qadd(to_cnum(b), Op[a])
Base.:*(b::Coefficient, a::QSym) = a * b

function Base.:*(a::QAdd, b::Coefficient)
    b isa Number && isone(b) && return a
    cb = to_cnum(b)
    d = QTermDict()
    for (term, c) in a.arguments
        new_c = c * cb
        iszero(new_c) && continue
        d[copy_key(term)] = new_c
    end
    return QAdd(d, drop_unused_indices(d, a.indices))
end
Base.:*(a::Coefficient, b::QAdd) = b * a

function Base.:+(a::QSym, b::QSym)
    d = QTermDict()
    addto!(d, Op[a], CNUM_ONE)
    addto!(d, Op[b], CNUM_ONE)
    return QAdd(d, EMPTY_INDICES)
end

function Base.:+(a::QAdd, b::QSym)
    d = copy_args(a.arguments)
    addto!(d, Op[b], CNUM_ONE)
    return QAdd(d, drop_unused_indices(d, a.indices))
end
Base.:+(a::QSym, b::QAdd) = b + a

function Base.:+(a::QAdd, b::QAdd)
    d = copy_args(a.arguments)
    merge_into!(d, b)
    return QAdd(d, drop_unused_indices(d, merge_unique(a.indices, b.indices)))
end

function Base.:+(a::QSym, b::Coefficient)
    d = QTermDict()
    addto!(d, Op[a], CNUM_ONE)
    addto!(d, EMPTY_OPS, to_cnum(b))
    return QAdd(d, EMPTY_INDICES)
end
Base.:+(a::Coefficient, b::QSym) = b + a

function Base.:+(a::QAdd, b::Coefficient)
    b isa Number && iszero(b) && return a
    d = copy_args(a.arguments)
    addto!(d, EMPTY_OPS, to_cnum(b))
    return QAdd(d, copy(a.indices))
end
Base.:+(a::Coefficient, b::QAdd) = b + a

const ZERO_QADD = QAdd(QTermDict(), Index[])
zero_qadd() = ZERO_QADD
Base.zero(::Type{QAdd}) = zero_qadd()
Base.zero(::QAdd) = zero_qadd()

Base.:+(a::QSym) = to_qadd(a)
Base.:+(a::QAdd) = a

Base.:-(a::QSym) = single_qadd(CNUM_NEG1, Op[a])

function Base.:-(a::QAdd)
    d = QTermDict()
    for (term, c) in a.arguments
        d[copy_key(term)] = -c
    end
    return QAdd(d, copy(a.indices))
end

Base.:-(a::QField, b::QField) = a + (-b)
Base.:-(a::QField, b::Coefficient) = a + (-b)
Base.:-(a::Coefficient, b::QField) = a + (-b)

Base.:/(a::Union{QSym, QAdd}, b::Number) =
    b isa Integer && !(b isa Bool) && !iszero(b) ? a * (1 // b) : a * inv(b)
Base.:/(a::Union{QSym, QAdd}, b::Coefficient) = a * inv(b)

Base.://(a::Union{QSym, QAdd}, b::Integer) = a * (1 // b)
Base.://(a::Union{QSym, QAdd}, b::Coefficient) = a / b

function Base.:^(a::QSym, n::Integer)
    n >= 0 || throw(ArgumentError("Negative powers not supported"))
    n == 0 && return single_qadd(CNUM_ONE, EMPTY_OPS)
    out = QTermDict()
    ops = Op[a for _ in 1:n]
    canonicalize!(out, ops, CNUM_ONE, EMPTY_NE)
    return QAdd(out, EMPTY_INDICES)
end

function Base.:^(a::QAdd, n::Integer)
    n >= 0 || throw(ArgumentError("Negative powers not supported"))
    n == 0 && return single_qadd(CNUM_ONE, EMPTY_OPS)
    result = a
    for _ in 2:n
        result = result * a
    end
    return result
end

function merge_into!(d::QTermDict, r::QAdd)
    for (term, c) in r.arguments
        addto_key!(d, copy_key(term), c)
    end
    return r.indices
end
