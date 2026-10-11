struct QAddBuilder
    args::QTermDict
    indices::Vector{Index}
end
QAddBuilder() = QAddBuilder(QTermDict(), Index[])

function merge_indices!(dest::Vector{Index}, src::Vector{Index})
    for x in src
        x ∉ dest && push!(dest, x)
    end
    return dest
end

function accumulate!(b::QAddBuilder, x::QAdd)
    for (term, c) in x.arguments
        addto_key!(b.args, copy_key(term), c)
    end
    isempty(x.indices) || merge_indices!(b.indices, x.indices)
    return b
end
accumulate!(b::QAddBuilder, x::QSym) = (addto!(b.args, Op[x], CNUM_ONE); b)
function accumulate!(b::QAddBuilder, x::Coefficient)
    x isa Number && iszero(x) && return b
    addto!(b.args, EMPTY_OPS, to_cnum(x))
    return b
end

function build(b::QAddBuilder)
    isempty(b.args) && return zero_qadd()
    return QAdd(b.args, drop_unused_indices(b.args, b.indices))
end

function MA.operate(
        ::typeof(sum), a::AbstractArray{<:QAdd}; init = zero(QAdd),
    )
    isempty(a) && return init
    if init isa QAdd || init isa QSym || init isa Coefficient
        b = QAddBuilder()
        accumulate!(b, init)
        for x in a
            accumulate!(b, x)
        end
        return build(b)
    end
    return mapreduce(identity, MA.add!!, a; init)
end
MA.promote_operation(::typeof(sum), ::Type{<:AbstractArray{<:QAdd}}) = QAdd

function Base.sum(a::AbstractArray{<:QAdd})
    isempty(a) && return zero(QAdd)
    b = QAddBuilder()
    for x in a
        accumulate!(b, x)
    end
    return build(b)
end
function Base.sum(f, a::AbstractArray{<:QAdd})
    isempty(a) && return zero(QAdd)
    b = QAddBuilder()
    for x in a
        accumulate!(b, f(x))
    end
    return build(b)
end
Base.reduce(::typeof(+), a::Vector{<:QAdd}) = sum(a)
