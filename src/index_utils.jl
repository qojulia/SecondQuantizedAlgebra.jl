# get_indices functions
function get_indices(term::QMul)
    args_nc = filter(x -> x isa IndexedOperator, term.args_nc)
    return unique(vcat(get_indices(args_nc), get_indices(term.arg_c)))
end
get_indices(a::IndexedOperator) = [a.ind]
get_indices(vec::AbstractVector) = unique(vcat(get_indices.(vec)...))
get_indices(term::SpecialIndexedTerm) = get_indices(term.term)
get_indices(term::IndexedAverageSum) = get_indices(term.term)
get_indices(term::IndexedAverageDoubleSum) = get_indices(term.innerSum)
function get_indices(term::SQABasicSymbolic)
    if is_symtype(term, IndexedAverageSum) && SymbolicUtils.hasmetadata(term, IndexedAverageSum)
        return get_indices(TermInterface.metadata(term)[IndexedAverageSum])
    elseif is_symtype(term, IndexedAverageDoubleSum) &&
        SymbolicUtils.hasmetadata(term, IndexedAverageDoubleSum)
        return get_indices(TermInterface.metadata(term)[IndexedAverageDoubleSum])
    elseif is_symtype(term, SpecialIndexedAverage) &&
        SymbolicUtils.hasmetadata(term, SpecialIndexedAverage)
        return get_indices(TermInterface.metadata(term)[SpecialIndexedAverage])
    end
    if is_symtype(term, DoubleIndexedVariable) &&
       SymbolicUtils.hasmetadata(term, DoubleIndexedVariable)
        meta = TermInterface.metadata(term)[DoubleIndexedVariable]
        return unique([meta.ind1, meta.ind2])
    elseif is_symtype(term, IndexedVariable)
        return [TermInterface.metadata(term)[IndexedVariable].ind]
    elseif is_average(term)
        return get_indices(unwrap_const(arguments(term)[1]))
    elseif iscall(term)
        return get_indices(map(unwrap_const, arguments(term)))
    end
    return Index[]
end

const Sums = Union{SingleSum,DoubleSum}
# get_indices(x::Sums) = unique(get_indices(arguments(x)))
get_indices(x::SingleSum) = get_indices(x.term)
get_indices(x::DoubleSum) = get_indices(x.innerSum.term)
get_indices(x::Number) = []
function get_indices(term)
    return iscall(term) ? get_indices(map(unwrap_const, arguments(term))) : []
end

#Usability functions:
Σ(a, b) = DoubleSum(a, b)  #Double-Sum here, because if variable a is not a single sum it will create a single sum anyway
Σ(a, b, c; kwargs...) = DoubleSum(a, b, c; kwargs...)
∑(args...; kwargs...) = Σ(args...; kwargs...)

IndexedOperator(x::IndexableOps, numb::Int64) = NumberedOperator(x, numb)
function IndexedOperator(x::QTerm, numb::Int64) # σ(1,1,2)
    f = SymbolicUtils.operation(x)
    args = SymbolicUtils.arguments(x)
    f([NumberedOperator(arg, numb) for arg in args]...)
end

#Numeric Conversion of NumberedOperators
function to_numeric(
    op::NumberedOperator,
    b::QuantumOpticsBase.CompositeBasis;
    ranges::Vector{Int64}=Int64[],
    kwargs...,
)
    if isempty(ranges)
        error(
            "When calling to_numeric for indexed Operators, specification of the \"ranges\" keyword is needed! This keyword requires a vector of Integers, which specify the maximum range of the index for each hilbertspace.",
        )
    end
    h = hilbert(op)
    if h isa ProductSpace
        if length(h.spaces) != length(ranges)
            error("Unequal length of hilbertspaces and ranges!")
        end
    else
        if length(ranges) != 1
            error("Wrong number of entries in ranges!")
        end
    end
    start = 0
    if h !== nothing #this is fine here since there are assertions above
        aon_ = acts_on(op)
        for i in 1:(aon_ - 1)
            start = start + ranges[i]
        end
    end
    aon = 0
    if start == 0
        aon = getNumber(op)[1] - 1
    else
        aon = op.numb + start
    end
    op_ = _to_numeric(op.op, b.bases[aon]; kwargs...)
    return QuantumOpticsBase.LazyTensor(b, [aon], (op_,))
end

function inadjoint(q::QMul)
    qad = adjoint(q)
    inorder!(qad)
    return qad
end
inadjoint(op::QNumber) = adjoint(op)
inadjoint(s::SQABasicSymbolic) = _adjoint(s)
inadjoint(x) = adjoint(x)

function inorder!(q::QMul)
    sort!(q.args_nc; by=get_numbers)
    sort!(q.args_nc; by=getIndName)
    sort!(q.args_nc; by=acts_on)
    return merge_commutators(q.arg_c, q.args_nc)
end
function inorder!(v::T) where {T<:SymbolicUtils.BasicSymbolic}
    if SymbolicUtils.iscall(v)
        f = SymbolicUtils.operation(v)
        args = map(arg -> inorder!(unwrap_const(arg)), SymbolicUtils.arguments(v))
        return TermInterface.maketerm(T, f, args, TermInterface.metadata(v))
    end
    return v
end
inorder!(x) = x

function inorder!(v::Average)
    f = operation(v)
    if f == conj
        return conj(inorder!(unwrap_const(arguments(v)[1])))
    end
    return average(inorder!(unwrap_const(arguments(v)[1])))
end

get_numbers(term::QMul) = unique(vcat(get_numbers.(term.args_nc)...))
get_numbers(term::Average) = get_numbers(unwrap_const(arguments(term)[1]))
get_numbers(x::NumberedOperator) = [x.numb]
get_numbers(x::Vector) = unique(vcat(get_numbers.(x)...))
get_numbers(x) = []

getNumber(x::NumberedOperator) = [acts_on(x) + x.numb]
getNumber(x::QMul) = acts_on(x)
getNumber(x) = [acts_on(x)] # this is so that, any other operator still behaves the same as before
