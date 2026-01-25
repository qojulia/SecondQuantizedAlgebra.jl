"""
    QNumber

Abstract type representing any expression involving operators.
"""
abstract type QNumber end

"""
    QSym <: QNumber

Abstract type representing fundamental operator types.
"""
abstract type QSym <: QNumber end

# Generic hash fallback for interface -- this will be slow
function Base.hash(op::T, h::UInt) where {T<:QSym}
    n = fieldcount(T)
    if n == 3
        # These three fields need to be defined for any QSym
        return hash(T, hash(op.hilbert, hash(op.name, hash(op.aon, h))))
    else
        # If there are more we'll need to iterate through
        h_ = copy(h)
        for k in n:-1:4
            if fieldname(typeof(op), k) !== :metadata
                h_ = hash(getfield(op, k), h_)
            end
        end
        return hash(T, hash(op.hilbert, hash(op.name, hash(op.aon, h_))))
    end
end

"""
    QTerm <: QNumber

Abstract type representing noncommutative expressions.
"""
abstract type QTerm <: QNumber end

Base.isless(a::QSym, b::QSym) = a.name < b.name

## Interface for SymbolicUtils

TermInterface.head(::QNumber) = :call
SymbolicUtils.iscall(::QSym) = false
SymbolicUtils.iscall(::QTerm) = true
SymbolicUtils.iscall(::Type{T}) where {T<:QTerm} = true
TermInterface.metadata(x::QNumber) = x.metadata

# Symbolic type promotion
SymbolicUtils.promote_symtype(f, Ts::Type{<:QNumber}...) = promote_type(Ts...)
SymbolicUtils.promote_symtype(f, T::Type{<:QNumber}, Ts...) = T
SymbolicUtils.promote_symtype(f, T::Type{<:QNumber}, S::Type{<:Number}) = T
SymbolicUtils.promote_symtype(f, T::Type{<:Number}, S::Type{<:QNumber}) = S
function SymbolicUtils.promote_symtype(f, T::Type{<:QNumber}, S::Type{<:QNumber})
    promote_type(T, S)
end

SymbolicUtils.symtype(x::T) where {T<:QNumber} = T

# Standard simplify and expand functions
function _has_index_sums(x)
    if x isa SingleSum || x isa DoubleSum || x isa SpecialIndexedTerm
        return true
    end
    if x isa QTerm && SymbolicUtils.iscall(x)
        return any(_has_index_sums, SymbolicUtils.arguments(x))
    end
    if x isa SQABasicSymbolic && SymbolicUtils.iscall(x)
        return any(_has_index_sums, SymbolicUtils.arguments(x))
    end
    return false
end

function _simplify_qnumber_noavg(x::QNumber; kwargs...)
    if x isa SingleSum || x isa DoubleSum
        return SymbolicUtils.simplify(x)
    elseif x isa QMul || x isa QAdd
        f = SymbolicUtils.operation(x)
        args = map(
            arg -> SymbolicUtils.simplify(arg; kwargs...), SymbolicUtils.arguments(x)
        )
        return f(args...)
    end
    return x
end

function SymbolicUtils.simplify(x::QNumber; kwargs...)
    if _has_index_sums(x)
        return _simplify_qnumber_noavg(x; kwargs...)
    end
    avg = average(x)
    avg_ = SymbolicUtils.simplify(avg; kwargs...)
    return undo_average(avg_)
end

function Symbolics.expand(x::QNumber; kwargs...)
    expansion = average(x)
    expansion_ = SymbolicUtils.expand(expansion; kwargs...)
    return undo_average(expansion_)
end

function Symbolics.substitute(x::QNumber, dict; kwargs...)
    if x isa QMul
        arg_c = Symbolics.substitute(x.arg_c, dict; kwargs...)
        args_sub = map(arg -> Symbolics.substitute(arg, dict; kwargs...), x.args_nc)
        args_c = filter(arg -> !(arg isa QNumber), args_sub)
        args_nc = filter(arg -> arg isa QNumber, args_sub)
        if any(arg -> isequal(arg, 0) || SymbolicUtils._iszero(arg), args_c)
            return 0
        end
        arg_c = if isempty(args_c)
            arg_c
        elseif length(args_c) == 1
            arg_c * args_c[1]
        else
            arg_c * *(args_c...)
        end
        if isequal(arg_c, 0) || SymbolicUtils._iszero(arg_c)
            return 0
        end
        isempty(args_nc) && return arg_c
        return QMul(arg_c, args_nc)
    elseif x isa QAdd
        args = map(arg -> Symbolics.substitute(arg, dict; kwargs...), x.arguments)
        if all(arg -> !(arg isa QNumber), args)
            return +(args...)
        end
        return QAdd(args)
    elseif x isa SingleSum
        term = Symbolics.substitute(x.term, dict; kwargs...)
        sum_index = Symbolics.substitute(x.sum_index, dict; kwargs...)
        neis = map(nei -> Symbolics.substitute(nei, dict; kwargs...), x.non_equal_indices)
        return SingleSum(term, sum_index, neis, x.metadata)
    elseif x isa DoubleSum
        inner = Symbolics.substitute(x.innerSum, dict; kwargs...)
        sum_index = Symbolics.substitute(x.sum_index, dict; kwargs...)
        neis = map(nei -> Symbolics.substitute(nei, dict; kwargs...), x.NEI)
        return DoubleSum(inner, sum_index, neis; metadata=x.metadata)
    end
    return x
end

function Symbolics.substitute(x::QSym, dict; kwargs...)
    haskey(dict, x) && return dict[x]
    return x
end

## End of interface

## Methods
import Base: *, +, -
const SNuN = Union{<:SymbolicUtils.BasicSymbolic{SQA_VARTYPE},<:Number}

Base.:~(a::QNumber, b::QNumber) = Symbolics.Equation(a, b)

## Multiplication
"""
    QMul <: QTerm

Represent a multiplication involving [`QSym`](@ref) types.

Fields:
======

* arg_c: The commutative prefactor.
* args_nc: A vector containing all [`QSym`](@ref) types.
"""
struct QMul{M} <: QTerm
    arg_c
    args_nc::Vector{Any}
    metadata::M
    function QMul{M}(arg_c, args_nc, metadata) where {M}
        if SymbolicUtils._isone(arg_c) && length(args_nc)==1
            return args_nc[1]
        elseif any(arg -> isequal(arg, 0) || SymbolicUtils._iszero(arg), args_nc) ||
            isequal(arg_c, 0) ||
            SymbolicUtils._iszero(arg_c)
            return 0
        else
            return new(arg_c, args_nc, metadata)
        end
    end
end
QMul(arg_c, args_nc; metadata::M=NO_METADATA) where {M} = QMul{M}(arg_c, args_nc, metadata)
Base.hash(q::QMul, h::UInt) = hash(QMul, hash(q.arg_c, hashvec(q.args_nc, h)))

SymbolicUtils.operation(::QMul) = (*)
SymbolicUtils.arguments(a::QMul) = vcat(a.arg_c, a.args_nc)

function SymbolicUtils.simplify(a::QMul; kwargs...)
    arg_c = SymbolicUtils.simplify(a.arg_c; kwargs...)
    args_nc = map(arg -> SymbolicUtils.simplify(arg; kwargs...), a.args_nc)
    return QMul(arg_c, args_nc; metadata=a.metadata)
end

function TermInterface.maketerm(::Type{<:QMul}, ::typeof(*), args, metadata)
    args_ = map(unwrap_const, args)
    args_c = filter(x->!(x isa QNumber), args_)
    args_nc = filter(x->x isa QNumber, args_)
    arg_c = *(args_c...)
    isempty(args_nc) && return arg_c
    return QMul(arg_c, args_nc; metadata)
end

TermInterface.metadata(a::QMul) = a.metadata

function Base.adjoint(q::QMul)
    args_nc = map(_adjoint, q.args_nc)
    reverse!(args_nc)
    sort!(args_nc; by=acts_on)
    return QMul(_conj(q.arg_c), args_nc; q.metadata)
end

function _coeff_equal(a, b)
    isequal(a, b) && return true
    diff = SymbolicUtils.simplify(SymbolicUtils.expand(a - b))
    return SymbolicUtils._iszero(diff) || isequal(diff, 0)
end

function Base.isequal(a::QMul, b::QMul)
    _coeff_equal(a.arg_c, b.arg_c) || return false
    length(a.args_nc)==length(b.args_nc) || return false
    for (arg_a, arg_b) in zip(a.args_nc, b.args_nc)
        isequal(arg_a, arg_b) || return false
    end
    return true
end

function *(a::QSym, b::QSym)
    check_hilbert(a, b)
    args = [a, b]
    sort!(args; by=acts_on)
    QMul(1, args)
end

function *(a::QSym, b::SNuN)
    SymbolicUtils._iszero(b) && return b
    SymbolicUtils._isone(b) && return a
    return QMul(b, [a])
end
*(b::SNuN, a::QNumber) = a*b

function *(a::QMul, b::SNuN)
    SymbolicUtils._iszero(b) && return b
    SymbolicUtils._isone(b) && return a
    arg_c = a.arg_c * b
    return QMul(arg_c, a.args_nc)
end

function *(a::QSym, b::QMul)
    check_hilbert(a, b)
    args_nc = vcat(a, b.args_nc)
    sort!(args_nc; by=acts_on)
    return merge_commutators(b.arg_c, args_nc)
end
function *(a::QMul, b::QSym)
    check_hilbert(a, b)
    args_nc = vcat(a.args_nc, b)
    sort!(args_nc; by=acts_on)
    return merge_commutators(a.arg_c, args_nc)
end

function *(a::QMul, b::QMul)
    check_hilbert(a, b)
    args_nc = vcat(a.args_nc, b.args_nc)
    sort!(args_nc; by=acts_on)
    arg_c = a.arg_c*b.arg_c
    return merge_commutators(arg_c, args_nc)
end

Base.:/(a::QNumber, b::SNuN) = (1/b) * a

function merge_commutators(arg_c, args_nc)
    #Added extra checks for 0 here
    if isequal(arg_c, 0) || 0 in args_nc
        return 0
    end
    i = 1
    was_merged = false
    while i<length(args_nc)
        if _ismergeable(args_nc[i], args_nc[i + 1])
            args_nc[i] = *(args_nc[i], args_nc[i + 1])
            iszero(args_nc[i]) && return 0
            deleteat!(args_nc, i+1)
            was_merged = true
        end
        i += 1
    end
    if was_merged
        return *(arg_c, args_nc...)
    else
        return QMul(arg_c, args_nc)
    end
end

function _ismergeable(a, b)
    isequal(acts_on(a), acts_on(b)) && ismergeable(a, b) && isequal(hilbert(a), hilbert(b))
end
ismergeable(a, b) = false

## Powers
function Base.:^(a::QNumber, n::Integer)
    iszero(n) && return 1
    isone(n) && return a
    return *((a for i in 1:n)...)
end

## Addition
"""
    QAdd <: QTerm

Represent an addition involving [`QNumber`](@ref) and other types.
"""
struct QAdd <: QTerm
    arguments::Vector{Any}
end

function Base.hash(q::T, h::UInt) where {T<:QAdd}
    hash(T, hashvec(sort(copy(q.arguments); by=string), h))
end
function Base.isequal(a::QAdd, b::QAdd)
    length(a.arguments)==length(b.arguments) || return false
    args_a = sort(copy(a.arguments); by=string)
    args_b = sort(copy(b.arguments); by=string)
    for (arg_a, arg_b) in zip(args_a, args_b)
        isequal(arg_a, arg_b) || return false
    end
    return true
end

SymbolicUtils.operation(::QAdd) = (+)
SymbolicUtils.arguments(a::QAdd) = a.arguments
function TermInterface.maketerm(::Type{<:QAdd}, ::typeof(+), args, metadata)
    QAdd(map(unwrap_const, args))
end
TermInterface.metadata(::QAdd) = NO_METADATA

Base.adjoint(q::QAdd) = QAdd(map(_adjoint, q.arguments))

function _qadd_key_tuple(x)
    if x isa Tuple
        x
    elseif x isa AbstractVector
        Tuple(x)
    elseif x === nothing
        ()
    else
        (x,)
    end
end

function SymbolicUtils.simplify(a::QAdd; kwargs...)
    args = map(arg -> SymbolicUtils.simplify(arg; kwargs...), a.arguments)
    flat = Any[]
    for arg in args
        if arg isa QAdd
            append!(flat, arg.arguments)
        else
            push!(flat, arg)
        end
    end
    dict = Dict{Any,Any}()
    for term in flat
        if term isa Number
            key = (:number,)
            coeff = term
        elseif term isa SingleSum && term.term isa QMul
            neis = term.non_equal_indices
            neis_tuple = _qadd_key_tuple(neis)
            args_nc = term.term.args_nc
            args_nc_tuple = _qadd_key_tuple(args_nc)
            key = (:singlesum, term.sum_index, neis_tuple, args_nc_tuple)
            coeff = term.term.arg_c
        elseif term isa QMul
            if length(term.args_nc) == 1
                key = (:term, term.args_nc[1])
                coeff = term.arg_c
            else
                key = (:qmul, _qadd_key_tuple(term.args_nc))
                coeff = term.arg_c
            end
        else
            key = (:term, term)
            coeff = 1
        end
        dict[key] = get(dict, key, 0) + coeff
    end
    new_args = Any[]
    for (key, coeff) in dict
        coeff = SymbolicUtils.simplify(SymbolicUtils.expand(coeff; kwargs...); kwargs...)
        if SymbolicUtils._iszero(coeff) || isequal(coeff, 0)
            continue
        end
        if key[1] == :singlesum
            _, sum_index, neis_tuple, args_nc_tuple = key
            args_nc = collect(args_nc_tuple)
            sum_term = length(args_nc) == 1 ? coeff * args_nc[1] : QMul(coeff, args_nc)
            term = SingleSum(sum_term, sum_index, collect(neis_tuple))
        elseif key[1] == :qmul
            args_nc = collect(key[2])
            term = length(args_nc) == 1 ? coeff * args_nc[1] : QMul(coeff, args_nc)
        elseif key[1] == :number
            term = coeff
        else
            term = isequal(coeff, 1) ? key[2] : coeff * key[2]
        end
        push!(new_args, term)
    end
    sort!(new_args; by=string)
    isempty(new_args) && return 0
    length(new_args) == 1 && return new_args[1]
    return QAdd(new_args)
end

-(a::QNumber) = -1*a
-(a, b::QNumber) = a + (-b)
-(a::QNumber, b) = a + (-b)
-(a::QNumber, b::QNumber) = a + (-b)

function +(a::QNumber, b::SNuN)
    SymbolicUtils._iszero(b) && return a
    return QAdd([a, b])
end
+(a::SNuN, b::QNumber) = +(b, a)
function +(a::QAdd, b::SNuN)
    SymbolicUtils._iszero(b) && return a
    args = vcat(a.arguments, b)
    return QAdd(args)
end

function +(a::QNumber, b::QNumber)
    check_hilbert(a, b)
    args = [a, b]
    return QAdd(args)
end

function +(a::QAdd, b::QNumber)
    check_hilbert(a, b)
    args = vcat(a.arguments, b)
    return QAdd(args)
end
function +(b::QNumber, a::QAdd)
    check_hilbert(a, b)
    args = vcat(a.arguments, b)
    return QAdd(args)
end
function +(a::QAdd, b::QAdd)
    check_hilbert(a, b)
    args = vcat(a.arguments, b.arguments)
    return QAdd(args)
end

function *(a::QAdd, b)
    check_hilbert(a, b)
    args = Any[a_ * b for a_ in a.arguments]
    flatten_adds!(args)
    isempty(args) && return 0
    q = QAdd(args)
    return q
end
function *(a::QNumber, b::QAdd)
    check_hilbert(a, b)
    args = Any[a * b_ for b_ in b.arguments]
    flatten_adds!(args)
    isempty(args) && return 0
    q = QAdd(args)
    return q
end

function *(a::QAdd, b::QAdd)
    check_hilbert(a, b)
    args = []
    for a_ in a.arguments, b_ in b.arguments
        push!(args, a_ * b_)
    end
    flatten_adds!(args)
    isempty(args) && return 0
    q = QAdd(args)
    return q
end

function flatten_adds!(args)
    i = 1
    while i <= length(args)
        if args[i] isa QAdd
            append!(args, args[i].arguments)
            deleteat!(args, i)
        elseif SymbolicUtils._iszero(args[i]) || isequal(args[i], 0)
            deleteat!(args, i)
        else
            i += 1
        end
    end
    return args
end

## Hilbert space checks
function check_hilbert(a::QNumber, b::QNumber)
    (hilbert(a) == hilbert(b)) ||
        error("Incompatible Hilbert spaces $(hilbert(a)) and $(hilbert(b))!")
end
check_hilbert(x, y) = nothing

hilbert(a::QSym) = a.hilbert
hilbert(a::QMul) = hilbert(a.args_nc[1])
function hilbert(a::QAdd)
    idx = findfirst(x->x isa QNumber, a.arguments)
    idx === nothing && print(a)
    hilbert(a.arguments[idx])
end

const AonType = Union{<:Int,<:ClusterAon}
"""
    acts_on(op)

Shows on which Hilbert space `op` acts. For [`QSym`](@ref) types, this
returns an Integer, whereas for a `Term` it returns a `Vector{Int}`
whose entries specify all subspaces on which the expression acts.
"""
acts_on(op::QSym) = op.aon
function acts_on(q::QMul)
    aon = AonType[]
    for arg in q.args_nc
        aon_ = acts_on(arg)
        aon_ ∈ aon || push!(aon, aon_)
    end
    return aon
end
function acts_on(q::QAdd)
    aon = AonType[]
    for arg in q.arguments
        append!(aon, acts_on(arg))
    end
    unique!(aon)
    sort!(aon)
    return aon
end
acts_on(x) = Int[]

Base.one(::T) where {T<:QNumber} = one(T)
Base.one(::Type{<:QNumber}) = 1
Base.isone(::QNumber) = false
Base.zero(::T) where {T<:QNumber} = zero(T)
Base.zero(::Type{<:QNumber}) = 0
Base.iszero(::QNumber) = false

"""
    @qnumbers

Convenience macro for the construction of operators.

Examples
========
```
julia> h = FockSpace(:fock)
ℋ(fock)

julia> @qnumbers a::Destroy(h)
(a,)

julia> h = FockSpace(:one) ⊗ FockSpace(:two)
ℋ(one) ⊗ ℋ(two)

julia> @qnumbers b::Destroy(h,2)
(b,)
```
"""
macro qnumbers(qs...)
    ex = Expr(:block)
    qnames = []
    for q in qs
        @assert q isa Expr && q.head==:(::)
        q_ = q.args[1]
        @assert q_ isa Symbol
        push!(qnames, q_)
        f = q.args[2]
        @assert f isa Expr && f.head==:call
        op = _make_operator(q_, f.args...)
        ex_ = Expr(:(=), esc(q_), op)
        push!(ex.args, ex_)
    end
    push!(ex.args, Expr(:tuple, map(esc, qnames)...))
    return ex
end

function _make_operator(name, T, h, args...)
    name_ = Expr(:quote, name)
    d = source_metadata(:qnumbers, name)
    return Expr(:call, T, esc(h), name_, args..., Expr(:kw, :metadata, Expr(:quote, d)))
end

get_operators(q::QSym) = [q]
function get_operators(q::QMul)
    ops = QSym[]
    seen_hashes = UInt[]
    for arg in q.args_nc
        h = hash(arg)
        if !(h ∈ seen_hashes)
            push!(seen_hashes, h)
            push!(ops, arg)
        end
    end
    return ops
end
function get_operators(q::QAdd)
    ops = QSym[]
    for arg in q.arguments
        append!(ops, get_operators(arg))
    end
    unique_ops!(ops)
    return ops
end
