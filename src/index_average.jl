#Main file for manipulating indexed averages and sums over averages.

const symbolics_terms = Union{SQABasicSymbolic,Number}
"""
    numberedVariable <: CNumber

abstract type for numbered Variables.
"""
abstract type numberedVariable <: CNumber end
"""

    IndexedAverageSum <: CNumber

Defines a symbolic summation over an average, or a multiplication of several averages, using one [`Index`](@ref) entity.

Fields:
======

* term: A multiplication of average terms.
* sum_index: The index, for which the summation will go over.
* non_equal_indices: (optional) A vector of indices, for which the summation-index can not be equal with.

"""
struct IndexedAverageSum <: CNumber
    term::symbolics_terms
    sum_index::Index
    non_equal_indices::Vector{IndexInt}
    metadata
    function IndexedAverageSum(
        term::symbolics_terms, sum_index::Index, non_equal_indices::Vector, metadata
    )
        neis_sym = ""
        if !(isempty(non_equal_indices))
            neis_sym = string("(", neis_sym)
            neis_sym = string(neis_sym, "$(sum_index.name)≠")
            neis_sym = string(neis_sym, writeNEIs(non_equal_indices))
            neis_sym = string(neis_sym, ")")
        end
        _metadata = new(term, sum_index, non_equal_indices, metadata)
        sym = SymbolicUtils.Sym{SQA_VARTYPE}(
            Symbol("∑($(sum_index.name)=1:$(sum_index.range))$(neis_sym)$(term)");
            type=IndexedAverageSum,
        )
        sym = SymbolicUtils.setmetadata(sym, typeof(_metadata), _metadata)
        sym = SymbolicUtils.setmetadata(sym, typeof(metadata), metadata)
        return sym
    end
end
Base.zero(::Type{IndexedAverageSum}) = 0
function IndexedAverageSum(
    term::symbolics_terms, sum_index::Index, non_equal_indices::Vector; metadata=NO_METADATA
)
    if sum_index ∉ get_indices(term)
        return (sum_index.range - length(non_equal_indices)) * term
    end
    prefact = 1.0
    if iscall(term)
        op = operation(term)
        args = map(unwrap_const, arguments(term))
        # move numbers outside of sum
        if op === *
            if args[1] isa Number
                prefact = args[1]
                args_nc = args[2:end]
                if length(args_nc) == 1
                    term = args_nc[1]
                else
                    term = *(args_nc...)
                end
            end
        end
        if op === +
            return sum(
                IndexedAverageSum(arg, sum_index, non_equal_indices; metadata=metadata) for
                arg in map(unwrap_const, arguments(term))
            )
        end
    end
    return prefact*IndexedAverageSum(term, sum_index, non_equal_indices, metadata)
end
IndexedAverageSum(x, args...; kwargs...) = average(SingleSum(x, args...; kwargs...))
IndexedAverageSum(x::Number) = x

"""

    IndexedAverageDoubleSum <: CNumber

Defines a symbolic summation over an [`IndexedAverageSum`](@ref), using a [`Index`](@ref) entity. This schematically represent a double-sum over a multiplication of averages.

Fields:
======

* innerSum: An [`IndexedAverageSum`](@ref) entity.
* sum_index: The index, for which the (outer) summation will go over.
* non_equal_indices: (optional) A vector of indices, for which the (outer) summation-index can not be equal with.

"""
struct IndexedAverageDoubleSum <: CNumber
    innerSum::SQABasicSymbolic
    sum_index::Index
    non_equal_indices::Vector{IndexInt}
    function IndexedAverageDoubleSum(
        term::SQABasicSymbolic, sum_index::Index, non_equal_indices
    )
        if sum_index ∉ get_indices(term)
            return (sum_index.range - length(non_equal_indices)) * term
        end
        if SymbolicUtils.iscall(term)
            op = operation(term)
            if op === +
                return sum(
                    IndexedAverageDoubleSum(arg, sum_index, non_equal_indices) for
                    arg in map(unwrap_const, arguments(term))
                )
            elseif op === *
                args = map(unwrap_const, arguments(term))
                if !isempty(args) && args[1] isa Number
                    return args[1] *
                           IndexedAverageDoubleSum(args[2], sum_index, non_equal_indices)
                end
            end
        end
        if !is_symtype(term, IndexedAverageSum)
            return IndexedAverageSum(term, sum_index, non_equal_indices)
        end
        _metadata = new(term, sum_index, non_equal_indices)
        neis_sym = ""
        if !(isempty(non_equal_indices))
            neis_sym = string("(", neis_sym)
            neis_sym = string(neis_sym, "$(sum_index.name)≠")
            neis_sym = string(neis_sym, writeNEIs(non_equal_indices))
            neis_sym = string(neis_sym, ")")
        end
        term_name = hasproperty(term, :name) ? String(term.name) : string(term)
        sym = SymbolicUtils.Sym{SQA_VARTYPE}(
            Symbol("∑($(sum_index.name):=1:$(sum_index.range))$(neis_sym)$(term_name)");
            type=IndexedAverageDoubleSum,
        )
        sym = SymbolicUtils.setmetadata(sym, typeof(_metadata), _metadata)
        return sym
    end
end
Base.zero(::Type{IndexedAverageDoubleSum}) = 0
function IndexedAverageDoubleSum(term::symbolics_terms, sum_index::Index, non_equal_indices)
    if iscall(term)
        op = operation(term)
        args = copy(map(unwrap_const, arguments(term)))
        param = 1.0
        if op === *
            if args[1] isa Number #put numbers out in front
                param = args[1]
                deleteat!(args, 1)
            end
            if length(args) == 1
                arg1 = unwrap_const(args[1])
                if arg1 isa SQABasicSymbolic && is_symtype(arg1, IndexedAverageSum)
                    return param*IndexedAverageDoubleSum(arg1, sum_index, non_equal_indices)
                end
            end
        elseif op === +
            return sum(
                IndexedAverageDoubleSum(arg, sum_index, non_equal_indices) for
                arg in map(unwrap_const, arguments(term))
            )
        end
    end
    return IndexedAverageSum(term, sum_index, non_equal_indices)
end
IndexedAverageDoubleSum(x, y, z) = IndexedAverageSum(x, y, z)

#Variables
"""

    SingleNumberedVariable <: numberedVariable

Defines a variable, associated with a Number. Used by [`value_map`](@ref)

Fields:
======

* name: The name of the variable.
* numb: An Integer Number.

"""
struct SingleNumberedVariable <: numberedVariable
    name::Symbol
    numb::Int64
    function SingleNumberedVariable(name, numb)
        sym_name = Symbol("$(name)_$(numb)")
        return Parameter(sym_name)
    end
end
"""

    DoubleNumberedVariable <: numberedVariable

Defines a variable, associated with two Numbers. Used by [`value_map`](@ref)

Fields:
======

* name: The name of the variable.
* numb1: An Integer Number.
* numb2: Another Integer Number.

"""
struct DoubleNumberedVariable <: numberedVariable
    name::Symbol
    numb1::IndexInt
    numb2::IndexInt
    function DoubleNumberedVariable(name, numb1, numb2; identical::Bool=true)
        if !(identical) && (numb1 == numb2)
            return 0
        end
        if typeof(numb1) == typeof(numb2) && numb1 isa Int64
            sym_name = Symbol("$(name)_{$(numb1),$(numb2)}")
            return Parameter(sym_name)
        else
            metadata = new(name, numb1, numb2)
            sym = SymbolicUtils.Sym{SQA_VARTYPE}(
                Symbol("$(name)_{$(numb1),$(numb2)}"); type=DoubleNumberedVariable
            )
            sym = SymbolicUtils.setmetadata(sym, typeof(metadata), metadata)
            return sym
        end
    end
end
struct SpecialIndexedAverage <: CNumber #An average-Term with special condition, for example l ≠ k; needed for correct calculus of Double indexed Sums
    term::symbolics_terms
    indexMapping::Vector{Tuple{IndexInt,IndexInt}}
    function SpecialIndexedAverage(term::Average, indexMapping)
        if !is_symtype(term, AvgSym)
            return invoke(
                SpecialIndexedAverage,
                Tuple{symbolics_terms,typeof(indexMapping)},
                term,
                indexMapping,
            )
        end
        if isempty(indexMapping)
            return term
        end
        if SymbolicUtils._iszero(unwrap_const(arguments(term)[1]))
            return 0
        end
        metadata = new(term, indexMapping)
        neis = writeIndexNEIs(indexMapping)
        sym = SymbolicUtils.Sym{SQA_VARTYPE}(
            Symbol("$(neis)$(term)"); type=SpecialIndexedAverage
        )
        sym = SymbolicUtils.setmetadata(sym, typeof(metadata), metadata)
        return sym
    end
end
Base.zero(::Type{SpecialIndexedAverage}) = 0
function SpecialIndexedAverage(term::symbolics_terms, indexMapping)
    if iscall(term)
        op = operation(term)
        args = copy(map(unwrap_const, arguments(term)))
        if op === *
            prefac = 1
            if args[1] isa Number
                prefac = args[1]
                deleteat!(args, 1)
            end
            if length(args) == 1
                return prefac * SpecialIndexedAverage(args[1], indexMapping)
            end
            specInds = [SpecialIndexedAverage(arg, indexMapping) for arg in args]
            return prefac * prod(specInds)
        elseif op === +
            return sum(SpecialIndexedAverage(arg, indexMapping) for arg in args)
        end
    end
    return term
end
SpecialIndexedAverage(x, args...) = x

const AvgSums = Union{
    SQABasicSymbolic,IndexedAverageSum,IndexedAverageDoubleSum,SpecialIndexedTerm
}
const AvgS = Union{IndexedAverageSum,IndexedAverageDoubleSum,SpecialIndexedTerm}

average(indOp::IndexedOperator) = SymbolicUtils._iszero(indOp) ? 0 : _average(indOp)
average(x::SpecialIndexedTerm) = SpecialIndexedAverage(average(x.term), x.indexMapping)
function average(indSum::SingleSum; kwargs...)
    IndexedAverageSum(average(indSum.term), indSum.sum_index, indSum.non_equal_indices)
end
function average(indDSum::DoubleSum)
    avg_inner = average(indDSum.innerSum)
    if avg_inner isa Number
        return (indDSum.sum_index.range - length(indDSum.NEI)) * avg_inner
    end
    return IndexedAverageDoubleSum(avg_inner, indDSum.sum_index, indDSum.NEI)
end

function undo_average(a::IndexedAverageSum)
    SingleSum(undo_average(a.term), a.sum_index, a.non_equal_indices)
end
function undo_average(a::IndexedAverageDoubleSum)
    DoubleSum(undo_average(a.innerSum), a.sum_index, a.non_equal_indices)
end
undo_average(a::SpecialIndexedAverage) = reorder(undo_average(a.term), a.indexMapping)
function undo_average(a::SQABasicSymbolic)
    if is_symtype(a, AvgSym)
        return unwrap_const(sqa_arguments(a)[1])
    elseif is_symtype(a, IndexedAverageSum) &&
        SymbolicUtils.hasmetadata(a, IndexedAverageSum)
        meta = TermInterface.metadata(a)[IndexedAverageSum]
        return undo_average(meta)
    elseif is_symtype(a, IndexedAverageDoubleSum) &&
        SymbolicUtils.hasmetadata(a, IndexedAverageDoubleSum)
        meta = TermInterface.metadata(a)[IndexedAverageDoubleSum]
        return undo_average(meta)
    elseif is_symtype(a, SpecialIndexedAverage) &&
        SymbolicUtils.hasmetadata(a, SpecialIndexedAverage)
        meta = TermInterface.metadata(a)[SpecialIndexedAverage]
        return undo_average(meta)
    elseif SymbolicUtils.iscall(a)
        f = SymbolicUtils.operation(a)
        args = map(arg -> undo_average(unwrap_const(arg)), sqa_arguments(a))
        return f(args...)
    end
    return a
end

# variable
IndexedVariable(x, numb::Int64) = SingleNumberedVariable(x, numb)
function IndexedVariable(x, num1::Int64, num2::Int64; kwargs...)
    DoubleNumberedVariable(x, num1, num2; kwargs...)
end
function IndexedVariable(name::Symbol, ind1::Index, ind2::Index; kwargs...)
    DoubleIndexedVariable(name, ind1, ind2; kwargs...)
end
function DoubleIndexedVariable(x, num1::Int64, num2::Int64; kwargs...)
    DoubleNumberedVariable(x, num1, num2; kwargs...)
end
# get_indices for SQABasicSymbolic is defined in index_utils.jl
#Symbolics functions
SymbolicUtils._iszero(a::IndexedAverageSum) = SymbolicUtils._iszero(a.term)
SymbolicUtils._isone(a::IndexedAverageSum) = SymbolicUtils._isone(a.term)

SymbolicUtils.iscall(a::IndexedAverageSum) = false
SymbolicUtils.iscall(a::IndexedAverageDoubleSum) = false

function writeIndexNEIs(neis::Vector{Tuple{IndexInt,IndexInt}})
    syms = ""
    syms = join([syms, "("])
    for i in 1:length(neis)
        if first(neis[i]) isa Index
            syms = join([syms, first(neis[i]).name])
        else
            syms = join([syms, first(neis[i])])
        end
        syms = join([syms, "≠"])
        if last(neis[i]) isa Index
            syms = join([syms, last(neis[i]).name])
        else
            syms = join([syms, last(neis[i])])
        end
        if i != length(neis)
            syms = join([syms, ","])
        end
    end
    syms = join([syms, ")"])
    return syms
end
function writeIndexNEIs(neis::Vector{Tuple{Index,Index}})
    writeIndexNEIs(convert(Vector{Tuple{IndexInt,IndexInt}}, neis))
end
function writeNeqs(vec::Vector{Tuple{Index,Int64}})
    syms = ""
    for i in 1:length(vec)
        syms = join([syms, "("])
        syms = join([syms, first(vec[i]).name])
        syms = join([syms, "≠"])
        syms = join([syms, last(vec[i])])
        syms = join([syms, ")"])
    end
    return syms
end # TODO: only used in tests, remove?

#Base functions
function Base.hash(a::IndexedAverageSum, h::UInt)
    return hash(
        IndexedAverageSum, hash(a.term, hash(a.sum_index, hash(a.non_equal_indices, h)))
    )
end

Base.isless(a::IndexedAverageSum, b::IndexedAverageSum) = a.sum_index < b.sum_index
function indexed_isequal(a::SQABasicSymbolic, b::SQABasicSymbolic)
    if SymbolicUtils.iscall(a) && SymbolicUtils.iscall(b)
        SymbolicUtils.operation(a) === SymbolicUtils.operation(b) || return false
        args_a = SymbolicUtils.arguments(a)
        args_b = SymbolicUtils.arguments(b)
        length(args_a) == length(args_b) || return false
        return all(
            indexed_isequal(unwrap_const(arg_a), unwrap_const(arg_b)) for
            (arg_a, arg_b) in zip(args_a, args_b)
        )
    end
    if is_symtype(a, IndexedAverageSum) && is_symtype(b, IndexedAverageSum)
        if SymbolicUtils.hasmetadata(a, IndexedAverageSum) &&
            SymbolicUtils.hasmetadata(b, IndexedAverageSum)
            a_meta = TermInterface.metadata(a)[IndexedAverageSum]
            b_meta = TermInterface.metadata(b)[IndexedAverageSum]
            return isequal(a_meta, b_meta)
        end
    elseif is_symtype(a, SpecialIndexedAverage) && is_symtype(b, SpecialIndexedAverage)
        if SymbolicUtils.hasmetadata(a, SpecialIndexedAverage) &&
            SymbolicUtils.hasmetadata(b, SpecialIndexedAverage)
            a_meta = TermInterface.metadata(a)[SpecialIndexedAverage]
            b_meta = TermInterface.metadata(b)[SpecialIndexedAverage]
            return isequal(a_meta.term, b_meta.term) &&
                   isequal(a_meta.indexMapping, b_meta.indexMapping)
        end
    end
    return isequal(a, b)
end
function Base.isequal(a::IndexedAverageSum, b::IndexedAverageSum)
    isequal(a.sum_index, b.sum_index) || return false
    isequal(a.term, b.term) || return false
    isequal(a.non_equal_indices, b.non_equal_indices) || return false
    return true
end
SymbolicUtils.arguments(op::IndexedAverageSum) = arguments(op.term)
SymbolicUtils.arguments(op::IndexedAverageDoubleSum) = op.innerSum
SymbolicUtils.arguments(op::SpecialIndexedAverage) = arguments(op.term)
function sqa_arguments(op::SQABasicSymbolic)
    if is_symtype(op, IndexedAverageSum) && SymbolicUtils.hasmetadata(op, IndexedAverageSum)
        return arguments(TermInterface.metadata(op)[IndexedAverageSum])
    elseif is_symtype(op, IndexedAverageDoubleSum) &&
        SymbolicUtils.hasmetadata(op, IndexedAverageDoubleSum)
        return TermInterface.metadata(op)[IndexedAverageDoubleSum].innerSum
    elseif is_symtype(op, SpecialIndexedAverage) &&
        SymbolicUtils.hasmetadata(op, SpecialIndexedAverage)
        return arguments(TermInterface.metadata(op)[SpecialIndexedAverage])
    end
    if SymbolicUtils.iscall(op)
        return SymbolicUtils.arguments(op)
    end
    return []
end

function sqa_simplify(sym::SQABasicSymbolic)
    if is_symtype(sym, SpecialIndexedAverage) &&
        SymbolicUtils.hasmetadata(sym, SpecialIndexedAverage)
        meta = TermInterface.metadata(sym)[SpecialIndexedAverage]
        return SpecialIndexedAverage(SymbolicUtils.simplify(meta.term), meta.indexMapping)
    end
    return SymbolicUtils.simplify(sym)
end

function Base.show(io::IO, indSum::IndexedAverageSum)
    write(io, "Σ", "($(indSum.sum_index.name)", "=1:$(indSum.sum_index.range))")
    if !(isempty(indSum.non_equal_indices))
        write(io, "($(indSum.sum_index.name) ≠ ")
        for i in 1:length(indSum.non_equal_indices)
            write(io, "$(indSum.non_equal_indices[i].name)")
            if i == length(indSum.non_equal_indices)
                write(io, ")")
            else
                write(io, ",")
            end
        end
    end
    Base.show(io, indSum.term)
end
function Base.show(io::IO, indSum::IndexedAverageDoubleSum)
    write(io, "Σ", "($(indSum.sum_index.name)", "=1:$(indSum.sum_index.range))")
    if !(isempty(indSum.non_equal_indices))
        write(io, "($(indSum.sum_index.name) ≠ ")
        for i in 1:length(indSum.non_equal_indices)
            write(io, "$(indSum.non_equal_indices[i].name)")
            if i == length(indSum.non_equal_indices)
                write(io, ")")
            else
                write(io, ",")
            end
        end
    end
    Base.show(io, indSum.innerSum)
end
function Base.show(io::IO, numbOp::NumberedOperator)
    Base.show(io, numbOp.op)
    Base.show(io, numbOp.numb)
end

function _to_expression(x::NumberedOperator)
    x.op isa Transition &&
        return :(NumberedOperator($(x.op.name), $(x.numb), $(x.op.i), $(x.op.j)))
    x.op isa Destroy && return :(NumberedDestroy($(x.op.name), $(x.numb)))
    x.op isa Create && return :(dagger(NumberedDestroy($(x.op.name), $(x.numb))))
end
function _to_expression(x::SQABasicSymbolic)
    if is_symtype(x, IndexedAverageSum)
        meta = TermInterface.metadata(x)[IndexedAverageSum]
        return :(IndexedAverageSum(
            $(_to_expression(meta.term)),
            $(meta.sum_index.name),
            $(meta.sum_index.range),
            $(writeNEIs(meta.non_equal_indices)),
        ))
    elseif is_symtype(x, SpecialIndexedAverage)
        meta = TermInterface.metadata(x)[SpecialIndexedAverage]
        return _to_expression(meta.term)
    elseif is_symtype(x, IndexedAverageDoubleSum)
        meta = TermInterface.metadata(x)[IndexedAverageDoubleSum]
        return :(IndexedAverageDoubleSum(
            $(_to_expression(meta.innerSum)),
            $(meta.sum_index.name),
            $(meta.sum_index.range),
            $(writeNEIs(meta.non_equal_indices)),
        ))
    elseif is_symtype(x, IndexedVariable) && SymbolicUtils.hasmetadata(x, IndexedVariable)
        meta = TermInterface.metadata(x)[IndexedVariable]
        return _to_expression(meta)
    elseif is_symtype(x, DoubleIndexedVariable) &&
        SymbolicUtils.hasmetadata(x, DoubleIndexedVariable)
        meta = TermInterface.metadata(x)[DoubleIndexedVariable]
        return _to_expression(meta)
    end
    return x
end

# ∨ https://github.com/qojulia/SecondQuantizedAlgebra.jl/issues/28
function SymbolicUtils.Add(::Type{SpecialIndexedAverage}, coeff, dict; kw...)
    SymbolicUtils.Add(SQA_VARTYPE, coeff, dict; type=CNumber, kw...)
end
function SymbolicUtils.Mul(::Type{SpecialIndexedAverage}, coeff, dict; kw...)
    SymbolicUtils.Mul(SQA_VARTYPE, coeff, dict; type=CNumber, kw...)
end

function SymbolicUtils.Add(::Type{IndexedAverageDoubleSum}, coeff, dict; kw...)
    SymbolicUtils.Add(SQA_VARTYPE, coeff, dict; type=CNumber, kw...)
end
function SymbolicUtils.Mul(::Type{IndexedAverageDoubleSum}, coeff, dict; kw...)
    SymbolicUtils.Mul(SQA_VARTYPE, coeff, dict; type=CNumber, kw...)
end

function SymbolicUtils.Add(::Type{IndexedAverageSum}, coeff, dict; kw...)
    SymbolicUtils.Add(SQA_VARTYPE, coeff, dict; type=CNumber, kw...)
end
function SymbolicUtils.Mul(::Type{IndexedAverageSum}, coeff, dict; kw...)
    SymbolicUtils.Mul(SQA_VARTYPE, coeff, dict; type=CNumber, kw...)
end

function SymbolicUtils.Add(::Type{IndexedVariable}, coeff, dict; kw...)
    SymbolicUtils.Add(SQA_VARTYPE, coeff, dict; type=CNumber, kw...)
end
function SymbolicUtils.Mul(::Type{IndexedVariable}, coeff, dict; kw...)
    SymbolicUtils.Mul(SQA_VARTYPE, coeff, dict; type=CNumber, kw...)
end

function SymbolicUtils.Add(::Type{DoubleIndexedVariable}, coeff, dict; kw...)
    SymbolicUtils.Add(SQA_VARTYPE, coeff, dict; type=CNumber, kw...)
end
function SymbolicUtils.Mul(::Type{DoubleIndexedVariable}, coeff, dict; kw...)
    SymbolicUtils.Mul(SQA_VARTYPE, coeff, dict; type=CNumber, kw...)
end

#function that creates an array consisting of all possible number values for each index given
#ind_vec should be sorted beforehand
function create_index_arrays(ind_vec, ranges)
    if length(ind_vec) == 1
        return ranges[1]
    end
    @assert length(ind_vec) == length(ranges)
    array = unique(collect(Iterators.product(ranges...)))
    length(ind_vec) == 1 && return array
    length(unique(get_spec_hilb.(ind_vec))) == length(ind_vec) && return array #every index has its own specHilb
    for vec in get_not_allowed(ind_vec)
        array = array[Bool[all_different(array[i], vec) for i in 1:length(array)]]
    end
    return collect(array)
end
all_different(x, vec) = length(unique(getindex(x, vec))) == length(getindex(x, vec))
function get_not_allowed(ind_vec)
    spec_hilbs = get_spec_hilb.(ind_vec)
    not_allowed = []
    for ind in ind_vec
        indices = findall(x -> isequal(x, ind.aon), spec_hilbs)
        length(indices) == 1 && continue
        if indices ∉ not_allowed
            push!(not_allowed, indices)
        end
    end
    return not_allowed
end
get_spec_hilb(ind::Index) = ind.aon

#this is the new method, insert values directly into the average before calculating anything, simplifies evaluation afterwards extremely
#function for inserting index, k -> 1,2,...,N
"""
    insert_index(term,ind::Index,value::Int)

Function, that inserts an integer value for a index in a specified term.
This function creates Numbered- Variables/Operators/Sums upon calls.

Examples
========

    insert_index(σⱼ²¹,j,1) = σ₁²¹

"""
function insert_index(term::SQABasicSymbolic, ind::Index, value::Int64)
    if is_average(term)
        f = operation(term)
        if f == conj
            return conj(insert_index(unwrap_const(arguments(term)[1]), ind, value))
        end
        return average(inorder!(insert_index(unwrap_const(arguments(term)[1]), ind, value)))
    elseif is_symtype(term, IndexedAverageSum)
        meta = TermInterface.metadata(term)[IndexedAverageSum]
        if ind == meta.sum_index
            error("cannot exchange summation index with number!")
        end
        if ind in meta.non_equal_indices
            newNEI = filter(x->!isequal(x, ind), meta.non_equal_indices)
            push!(newNEI, value)
            return IndexedAverageSum(
                insert_index(meta.term, ind, value), meta.sum_index, newNEI
            )
        else
            return IndexedAverageSum(
                insert_index(meta.term, ind, value), meta.sum_index, meta.non_equal_indices
            )
        end
    elseif is_symtype(term, IndexedAverageDoubleSum)
        meta = TermInterface.metadata(term)[IndexedAverageDoubleSum]
        inner = insert_index(meta.innerSum, ind, value)
        return IndexedAverageDoubleSum(inner, meta.sum_index, meta.non_equal_indices)
    elseif is_symtype(term, SpecialIndexedAverage)
        meta = TermInterface.metadata(term)[SpecialIndexedAverage]
        newterm = insert_index(meta.term, ind, value)
        newMapping = Tuple{IndexInt,IndexInt}[]
        for tuple in meta.indexMapping
            if first(tuple) == ind
                if last(tuple) == value
                    return 0
                end
                push!(newMapping, (value, last(tuple)))
            elseif last(tuple) == ind
                if first(tuple) == value
                    return 0
                end
                push!(newMapping, (first(tuple), value))
            else
                push!(newMapping, tuple)
            end
        end
        filter!(x -> !(first(x) isa Int64 && last(x) isa Int64), newMapping)
        return SpecialIndexedAverage(newterm, newMapping)
    elseif is_symtype(term, DoubleIndexedVariable)
        data = TermInterface.metadata(term)[DoubleIndexedVariable]
        if data.ind1 == ind && data.ind2 == ind
            return DoubleNumberedVariable(data.name, value, value)
        elseif data.ind1 == ind
            return DoubleNumberedVariable(data.name, value, data.ind2)
        elseif data.ind2 == ind
            return DoubleNumberedVariable(data.name, data.ind1, value)
        end
        return term
    elseif is_symtype(term, DoubleNumberedVariable)
        if iscall(term)
            op = operation(term)
            if op === *
                return prod(
                    insert_index(unwrap_const(arg), ind, value) for arg in arguments(term)
                )
            elseif op === +
                return sum(
                    insert_index(unwrap_const(arg), ind, value) for arg in arguments(term)
                )
            elseif op === ^
                return insert_index(
                    unwrap_const(arguments(term)[1]), ind, value
                )^(unwrap_const(arguments(term)[2]))
            end
        end
        data = TermInterface.metadata(term)[DoubleNumberedVariable]
        if data.numb1 isa Index && data.numb1 == ind
            return DoubleNumberedVariable(data.name, value, data.numb2)
        elseif data.numb2 isa Index && data.numb2 == ind
            return DoubleNumberedVariable(data.name, data.numb1, value)
        end
        return term
    elseif is_symtype(term, IndexedVariable)
        meta = TermInterface.metadata(term)[IndexedVariable]
        return meta.ind == ind ? SingleNumberedVariable(meta.name, value) : term
    elseif is_symtype(term, CNumber)
        if iscall(term)
            op = operation(term)
            if op === *
                return prod(
                    insert_index(unwrap_const(arg), ind, value) for arg in arguments(term)
                )
            elseif op === +
                return sum(
                    insert_index(unwrap_const(arg), ind, value) for arg in arguments(term)
                )
            elseif op === ^
                return insert_index(
                    unwrap_const(arguments(term)[1]), ind, value
                )^insert_index(unwrap_const(arguments(term)[2]), ind, value)
                # issue QuantumCumulants 198
            elseif op === /
                return insert_index(
                    unwrap_const(arguments(term)[1]), ind, value
                )/insert_index(unwrap_const(arguments(term)[2]), ind, value)
            elseif length(arguments(term)) == 1 # exp, sin, cos, ln, ... #TODO: write tests
                return op(insert_index(unwrap_const(arguments(term)[1]), ind, value))
            end
        end
        return term
    end
    return term
end
function insert_index(qmul::QMul, ind::Index, value::Int64)
    qmul.arg_c*prod(insert_index(arg, ind, value) for arg in qmul.args_nc)
end
function insert_index(eq::Symbolics.Equation, ind::Index, value::Int64)
    Symbolics.Equation(insert_index(eq.lhs, ind, value), insert_index(eq.rhs, ind, value))
end
function insert_index(term::IndexedOperator, ind::Index, value::Int64)
    term.ind == ind ? NumberedOperator(term.op, value) : term
end
insert_index(x, args...) = x
