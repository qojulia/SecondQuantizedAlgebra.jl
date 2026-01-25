#Base file for defining DoubleIndexedSums

"""
    DoubleSum <: QTerm

Defines a symbolic summation over another [`SingleSum`](@ref), using one [`Index`](@ref) entity. This corresponds to a double-summation over a multiplication of terms.

Fields:
======

* innerSum: A [`SingleSum`](@ref) entity.
* sum_index: The index, for which the (outer) summation will go over.
* NEI: (optional) A vector of indices, for which the (outer) summation-index can not be equal with.

"""
struct DoubleSum{M} <: QTerm
    innerSum::SingleSum
    sum_index::Index
    NEI::Vector{Index}
    metadata::M
    function DoubleSum(innerSum::SingleSum, sum_index::Index, NEI::Vector, metadata)
        try
            return new{typeof(metadata)}(innerSum, sum_index, NEI, metadata)
        catch e
            println(
                "Could not create DoubleSum with input: term= $(innerSum) ; sum_index=c$(sum_index) ; NEI= $(NEI) ; metadata= $(metadata)",
            )
            rethrow(e)
        end
    end
end
function DoubleSum(innerSum::SingleSum, sum_index::Index, NEI::Vector; metadata=NO_METADATA)
    if innerSum.sum_index == sum_index
        error("summation index is the same as the index of the inner sum")
    else
        extraterm = 0
        NEI_ = copy(NEI)
        if sum_index in innerSum.non_equal_indices
            for index in get_indices(innerSum.term)
                if isequal(index, innerSum.sum_index)
                    (innerSum.sum_index ∉ NEI_) && push!(NEI_, index)
                    continue
                end
                if index != sum_index && index ∉ NEI && isequal(index.aon, sum_index.aon)
                    extraterm =
                        extraterm + SingleSum(
                            change_index(innerSum.term, sum_index, index),
                            innerSum.sum_index,
                            innerSum.non_equal_indices,
                        )
                    push!(NEI_, index)
                end
            end
        end
        if innerSum.term isa QMul
            # put terms of the outer index in front
            indicesToOrder = sort([innerSum.sum_index, sum_index]; by=getIndName)
            newargs = order_by_index(innerSum.term.args_nc, indicesToOrder)
            qmul = 0
            if length(newargs) == 1
                qmul = *(innerSum.term.arg_c, newargs[1])
            else
                qmul = *(innerSum.term.arg_c, newargs...)
            end
            sort!(NEI_)
            innerSum_ = SingleSum(qmul, innerSum.sum_index, innerSum.non_equal_indices)
            if innerSum_ isa SingleSum
                if isequal(extraterm, 0) || SymbolicUtils._iszero(extraterm)
                    return DoubleSum(innerSum_, sum_index, NEI_, metadata)
                end
                return DoubleSum(innerSum_, sum_index, NEI_, metadata) + extraterm
            else
                return DoubleSum(innerSum_, sum_index, NEI_; metadata=metadata) + extraterm
            end
        else # this case is only reachable if innersum has only one operator ->  no extraterm anyway
            sort!(NEI)
            return DoubleSum(innerSum, sum_index, NEI, metadata)
        end
    end
end
function DoubleSum(
    innerSum::IndexedAdd, sum_index::Index, NEI::Vector; metadata=NO_METADATA
)
    if iscall(innerSum)
        op = operation(innerSum)
        if op === +
            sums = [
                DoubleSum(arg, sum_index, NEI; metadata=metadata) for
                arg in map(unwrap_const, arguments(innerSum))
            ]
            isempty(sums) && return 0
            length(sums) == 1 && return sums[1]
            return +(sums...)
        end
    end
    return SingleSum(innerSum, sum_index, NEI; metadata=metadata)
end
function DoubleSum(
    innerSum::SQABasicSymbolic, sum_index::Index, NEI::Vector; metadata=NO_METADATA
)
    if iscall(innerSum)
        op = operation(innerSum)
        if op === +
            sums = [
                DoubleSum(arg, sum_index, NEI; metadata=metadata) for
                arg in map(unwrap_const, arguments(innerSum))
            ]
            isempty(sums) && return 0
            length(sums) == 1 && return sums[1]
            return +(sums...)
        end
    end
    return SingleSum(innerSum, sum_index, NEI; metadata=metadata)
end
DoubleSum(x, ind::Index, NEI::Vector; metadata=NO_METADATA) = SingleSum(x, ind, NEI)
DoubleSum(x, ind::Index; metadata=NO_METADATA) = DoubleSum(x, ind, Index[])

#In this constructor the NEI is considered so, that all indices given in ind are unequal to any of the NEI
function DoubleSum(term::QMul, ind::Vector{Index}, NEI::Vector{Index}; metadata=NO_METADATA)
    if length(ind) != 2
        error("Can only create Double-Sum with 2 indices!")
    end
    return DoubleSum(SingleSum(term, ind[1], NEI), ind[2], NEI; metadata=metadata)
end
function DoubleSum(
    term, innerInd::Index, outerInd::Index; non_equal::Bool=false, metadata=NO_METADATA
)
    if non_equal
        innerSum = SingleSum(term, innerInd, [outerInd])
        return DoubleSum(innerSum, outerInd, []; metadata=metadata)
    else
        innerSum = SingleSum(term, innerInd, [])
        return DoubleSum(innerSum, outerInd, []; metadata=metadata)
    end
end

hilbert(elem::DoubleSum) = hilbert(elem.sum_index)
#multiplications
*(elem::Number, sum::DoubleSum) = DoubleSum(elem*sum.innerSum, sum.sum_index, sum.NEI)
*(sum::DoubleSum, elem::Number) = DoubleSum(sum.innerSum*elem, sum.sum_index, sum.NEI)
*(sum::DoubleSum, qmul::QMul) = qmul.arg_c*(*(sum, qmul.args_nc...))
function *(qmul::QMul, sum::DoubleSum)
    sum_ = sum
    for i in length(qmul.args_nc):-1:1
        sum_ = qmul.args_nc[i] * sum_
    end
    return qmul.arg_c*sum_
end

function _mul_elem_sum(elem, sum::DoubleSum)
    sum_ = SymbolicUtils.simplify(sum)
    if !(sum_ isa DoubleSum) # issue QuantumCumulants 223
        return elem*sum_
    end
    NEI = copy(sum.NEI)
    aon_sum = sum.sum_index.aon
    elem_ind = _get_index(elem)
    aon_elem = elem_ind.aon
    if aon_sum ≠ aon_elem # issue QuantumCumulants 256
        return DoubleSum(elem*sum.innerSum, sum.sum_index, NEI)
    end
    if elem_ind != sum.sum_index && elem_ind ∉ NEI
        if (sum.sum_index.aon != sum.innerSum.sum_index.aon) # indices for different ops
            if isequal(elem_ind.aon, sum.sum_index.aon)
                push!(NEI, elem_ind)
                addterm = SingleSum(
                    elem*change_index(sum.innerSum.term, sum.sum_index, elem_ind),
                    sum.innerSum.sum_index,
                    sum.innerSum.non_equal_indices,
                )
                return DoubleSum(elem*sum.innerSum, sum.sum_index, NEI) + addterm
            end
            return DoubleSum(elem*sum.innerSum, sum.sum_index, NEI)
        end
        NEI_ = [NEI..., elem_ind] # issue QuantumCumulants #169 (scaling of double sum)
        ds_term = DoubleSum(
            SingleSum(
                elem*sum.innerSum.term,
                sum.innerSum.sum_index,
                [sum.innerSum.non_equal_indices..., elem_ind],
            ),
            sum.sum_index,
            NEI_,
        )
        new_non_equal_indices1 = replace(sum.NEI, sum.innerSum.sum_index => elem_ind)
        ss_term1 = SingleSum(
            elem*change_index(sum.innerSum.term, sum.innerSum.sum_index, elem_ind),
            sum.sum_index,
            new_non_equal_indices1,
        )
        new_non_equal_indices2 = unique([
            replace(sum.innerSum.non_equal_indices, sum.sum_index => elem_ind)..., elem_ind
        ]) #issue QuantumCumulants #223
        ss_term2 = SingleSum(
            elem*change_index(sum.innerSum.term, sum.sum_index, elem_ind),
            sum.innerSum.sum_index,
            new_non_equal_indices2,
        )
        return ds_term + ss_term1 + ss_term2
    end
    return DoubleSum(elem*sum.innerSum, sum.sum_index, NEI)
end
function *(elem::IndexedOperator, sum::DoubleSum)
    return _mul_elem_sum(elem, sum)
end
function *(elem::SQABasicSymbolic, sum::DoubleSum)
    _is_indexed_symbolic(elem) ||
        return DoubleSum(elem*sum.innerSum, sum.sum_index, sum.NEI)
    return _mul_elem_sum(elem, sum)
end
function _mul_sum_elem(sum::DoubleSum, elem)
    sum_ = SymbolicUtils.simplify(sum)
    if !(sum_ isa DoubleSum) # issue QuantumCumulants 223
        return sum_*elem
    end
    NEI = copy(sum.NEI)
    aon_sum = sum.sum_index.aon
    elem_ind = _get_index(elem)
    aon_elem = elem_ind.aon
    if aon_sum ≠ aon_elem # issue QuantumCumulants 256
        return DoubleSum(sum.innerSum*elem, sum.sum_index, NEI)
    end
    if elem_ind != sum.sum_index && elem_ind ∉ NEI
        if (sum.sum_index.aon != sum.innerSum.sum_index.aon) # indices for different ops
            if isequal(elem_ind.aon, sum.sum_index.aon)
                push!(NEI, elem_ind)
                addterm = SingleSum(
                    change_index(sum.innerSum.term, sum.sum_index, elem_ind)*elem,
                    sum.innerSum.sum_index,
                    sum.innerSum.non_equal_indices,
                )
                return DoubleSum(sum.innerSum*elem, sum.sum_index, NEI) + addterm
            end
            return DoubleSum(sum.innerSum*elem, sum.sum_index, NEI)
        end
        function *(sum::DoubleSum, elem::IndexedOperator)
            return _mul_sum_elem(sum, elem)
        end
        function *(sum::DoubleSum, elem::SQABasicSymbolic)
            _is_indexed_symbolic(elem) ||
                return DoubleSum(sum.innerSum*elem, sum.sum_index, sum.NEI)
            return _mul_sum_elem(sum, elem)
        end
        NEI_ = [NEI..., elem_ind] # issue QuantumCumulants #169 (scaling of double sum)
        ds_term = DoubleSum(
            SingleSum(
                sum.innerSum.term*elem,
                sum.innerSum.sum_index,
                [sum.innerSum.non_equal_indices..., elem_ind],
            ),
            sum.sum_index,
            NEI_,
        )
        new_non_equal_indices1 = replace(sum.NEI, sum.innerSum.sum_index => elem_ind)
        ss_term1 = SingleSum(
            change_index(sum.innerSum.term, sum.innerSum.sum_index, elem_ind)*elem,
            sum.sum_index,
            new_non_equal_indices1,
        )
        # new_non_equal_indices2 = replace(sum.innerSum.non_equal_indices, sum.sum_index => elem.ind)
        new_non_equal_indices2 = unique([
            replace(sum.innerSum.non_equal_indices, sum.sum_index => elem_ind)..., elem_ind
        ]) #issue QuantumCumulants #223
        ss_term2 = SingleSum(
            change_index(sum.innerSum.term, sum.sum_index, elem_ind)*elem,
            sum.innerSum.sum_index,
            new_non_equal_indices2,
        )
        return ds_term + ss_term1 + ss_term2
    end #with else it does not work?
    return DoubleSum(sum.innerSum*elem, sum.sum_index, NEI)
end
function *(sum::DoubleSum, elem::SQABasicSymbolic)
    _is_indexed_symbolic(elem) ||
        return DoubleSum(sum.innerSum*elem, sum.sum_index, sum.NEI)
    return *(sum, elem::IndexedObSym)
end

function SymbolicUtils.simplify(a::DoubleSum)
    inner = SymbolicUtils.simplify(a.innerSum)
    if a.sum_index ∉ get_indices(inner)
        return (a.sum_index.range - length(a.NEI)) * inner
    end
    return DoubleSum(inner, a.sum_index, a.NEI; metadata=a.metadata)
end

acts_on(sum::DoubleSum) = acts_on(sum.innerSum)

_aon_vec(x) = (aon=acts_on(x); aon isa AbstractVector ? aon : [aon])
function commutator(a::QSym, b::SingleSum)
    isempty(intersect(_aon_vec(a), _aon_vec(b))) && return 0
    return _commutator(a, b)
end
commutator(a::SingleSum, b::QSym) = -commutator(b, a)
function commutator(a::QSym, b::DoubleSum)
    isempty(intersect(_aon_vec(a), _aon_vec(b))) && return 0
    return _commutator(a, b)
end
commutator(a::DoubleSum, b::QSym) = -commutator(b, a)

*(sum::DoubleSum, x::SQA_SumArg) = DoubleSum(sum.innerSum*x, sum.sum_index, sum.NEI)
*(x::SQA_SumArg, sum::DoubleSum) = DoubleSum(x*sum.innerSum, sum.sum_index, sum.NEI)
*(x::SNuN, sum::DoubleSum) = DoubleSum(x*sum.innerSum, sum.sum_index, sum.NEI)
function *(a::QAdd, b::DoubleSum)
    return DoubleSum(a * b.innerSum, b.sum_index, b.NEI)
end
function *(a::DoubleSum, b::QAdd)
    return DoubleSum(a.innerSum * b, a.sum_index, a.NEI)
end
function *(term::SpecialIndexedTerm, sum::DoubleSum)
    return reorder(term.term * sum, term.indexMapping)
end
function *(sum::DoubleSum, term::SpecialIndexedTerm)
    return reorder(sum * term.term, term.indexMapping)
end
function *(sum1::DoubleSum, sum2::DoubleSum)
    error("Multiplication of DoubleSum with DoubleSum is not defined")
end

SymbolicUtils.iscall(a::DoubleSum) = false
SymbolicUtils.arguments(a::DoubleSum) = SymbolicUtils.arguments(a.innerSum)
checkInnerSums(sum1::DoubleSum, sum2::DoubleSum) = ((sum1.innerSum + sum2.innerSum) == 0)
function reorder(dsum::DoubleSum, indexMapping::Vector{Tuple{Index,Index}})
    DoubleSum(reorder(dsum.innerSum, indexMapping), dsum.sum_index, dsum.NEI)
end

#Base functions
function Base.show(io::IO, elem::DoubleSum)
    write(io, "Σ", "($(elem.sum_index.name)=1:$(elem.sum_index.range))")
    if !(isempty(elem.NEI))
        write(io, "($(elem.sum_index.name)≠")
        for i in 1:length(elem.NEI)
            write(io, "$(elem.NEI[i].name))")
        end
    end
    show(io, elem.innerSum)
end
function Base.isequal(a::DoubleSum, b::DoubleSum)
    isequal(a.innerSum, b.innerSum) &&
        isequal(a.sum_index, b.sum_index) &&
        isequal(a.NEI, b.NEI)
end
function _to_expression(x::DoubleSum)
    :(DoubleSum(
        $(_to_expression(x.innerSum)),
        $(x.sum_index.name),
        $(x.sum_index.range),
        $(writeNEIs(x.NEI)),
    ))
end

function *(sum1::SingleSum, sum2::SingleSum; ind=nothing)
    if sum1.sum_index != sum2.sum_index
        term = sum1.term*sum2.term
        return DoubleSum(
            SingleSum(term, sum1.sum_index, sum1.non_equal_indices),
            sum2.sum_index,
            sum2.non_equal_indices,
        )
    else
        if !(ind isa Index)
            error("Specification of an extra Index is needed!")
        end
        term2 = change_index(sum2.term, sum2.sum_index, ind)
        return DoubleSum(
            SingleSum(sum1.term*term2, sum1.sum_index, sum1.non_equal_indices),
            ind,
            sum1.non_equal_indices,
        )
    end
end
