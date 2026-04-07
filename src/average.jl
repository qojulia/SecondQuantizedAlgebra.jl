"""
    AvgSym

Marker type for the `symtype` of averaged operator expressions.
Average nodes are `BasicSymbolic{SymReal}` `Term` nodes with `symtype(x) === AvgSym`.
"""
struct AvgSym <: Number end

"""Metadata key for summation indices on averaged expressions."""
struct SumIndices end

"""Metadata key for non-equal index pairs on averaged expressions."""
struct SumNonEqual end

"""
    AvgFunc

Singleton callable struct used as the `operation` of average `Term` nodes.
Using a custom struct (instead of a SymbolicUtils `Sym`) lets us define
`SymbolicUtils.show_call(::IO, ::AvgFunc, ...)` for `⟨op⟩` display
without type piracy.
"""
struct AvgFunc end

"""Singleton instance of [`AvgFunc`](@ref)."""
const sym_average = AvgFunc()

Base.nameof(::AvgFunc) = :avg
Base.show(io::IO, ::AvgFunc) = print(io, "avg")

function SymbolicUtils.show_call(
        io::IO, ::AvgFunc, x::SymbolicUtils.BasicSymbolic; kw...
    )
    args = SymbolicUtils.arguments(x)
    print(io, "⟨")
    for (i, arg) in enumerate(args)
        i > 1 && print(io, ", ")
        print(io, arg)
    end
    return print(io, "⟩")
end

"""
    is_average(x) -> Bool

Check whether `x` is a symbolic average `Term` node.
"""
is_average(::QField) = false
is_average(::Number) = false
function is_average(x::SymbolicUtils.BasicSymbolic)
    return SymbolicUtils.iscall(x) && SymbolicUtils.symtype(x) === AvgSym
end
is_average(x) = false

function _average(op::QField)
    return SymbolicUtils.Term{SymbolicUtils.SymReal}(sym_average, Any[op]; type = AvgSym)
end

"""
    average(expr)

Compute the symbolic average of a quantum operator expression.
Returns a `BasicSymbolic{SymReal}` scalar with `symtype = AvgSym` that participates
in Symbolics arithmetic.

Linearity: distributes over `QAdd` and pulls c-number prefactors out.
"""
function average end

average(op::QSym) = _average(op)
average(x::Number) = x
average(x::SymbolicUtils.BasicSymbolic) = x

function average(op::QAdd)
    result = Num(0)
    for (ops, c) in op.arguments
        if isempty(ops)
            result += c
        else
            # Wrap the operator product for _average.
            # For single ops, pass the QSym directly; for multi-op, pass a QAdd.
            inner = length(ops) == 1 ? only(ops) : _single_qadd(_CNUM_ONE, ops)
            avg = _average(inner)
            r, i = real(c), imag(c)
            if iszero(i)
                result += r * avg
            elseif iszero(r)
                result += im * i * avg
            else
                result += (r + im * i) * avg
            end
        end
    end
    if !isempty(op.indices)
        result = SymbolicUtils.setmetadata(result, SumIndices, op.indices)
        result = SymbolicUtils.setmetadata(result, SumNonEqual, op.non_equal)
    end
    return result
end

# --- undo_average ---

"""Wrap any value as a QAdd — ensures uniform return type."""
_to_qadd(x::QAdd) = x
_to_qadd(x::QSym) = _single_qadd(_CNUM_ONE, QSym[x])
_to_qadd(x::Number) = _single_qadd(_to_cnum(x), QSym[])
_to_qadd(x::CNum) = _single_qadd(x, QSym[])
_to_qadd(x::Num) = _single_qadd(_to_cnum(x), QSym[])
function _to_qadd(x::SymbolicUtils.BasicSymbolic)
    return _single_qadd(_to_cnum(x), QSym[])
end

"""Restore summation metadata from a SymbolicUtils node onto a QAdd."""
function _restore_sum_metadata(result::QAdd, x::SymbolicUtils.BasicSymbolic)
    if SymbolicUtils.hasmetadata(x, SumIndices)
        indices = SymbolicUtils.getmetadata(x, SumIndices)
        non_equal = SymbolicUtils.getmetadata(x, SumNonEqual)
        return QAdd(result.arguments, indices, non_equal)
    end
    return result
end

"""
    undo_average(expr) -> QAdd

Recursively strip `⟨...⟩` from a symbolic expression, recovering operator expressions.
Summation metadata (indices, non-equal constraints) is restored to the resulting `QAdd`.

Always returns `QAdd` for type stability:
- Scalars become single-term `QAdd`s with an empty operator sequence.
- Single operators become single-term `QAdd`s with unit prefactor.
- Sums of averages are reconstructed via `+` and `*` on `QAdd`s.
"""
function undo_average(x::SymbolicUtils.BasicSymbolic)::QAdd
    if SymbolicUtils.iscall(x)
        f = SymbolicUtils.operation(x)
        if f isa AvgFunc
            arg = SymbolicUtils.arguments(x)[1]
            result = SymbolicUtils.isconst(arg) ? arg.val : arg
            return _restore_sum_metadata(_to_qadd(result), x)
        elseif f === (+) || f === (*)
            args = map(undo_average, SymbolicUtils.arguments(x))
            result = f(args...)
            return _restore_sum_metadata(_to_qadd(result), x)
        else
            return _to_qadd(x)
        end
    else
        return _to_qadd(x)
    end
end

undo_average(x::Number)::QAdd = _single_qadd(_to_cnum(x), QSym[])
undo_average(x::Num)::QAdd = undo_average(SymbolicUtils.unwrap(x))
undo_average(x::QSym)::QAdd = _single_qadd(_CNUM_ONE, QSym[x])
undo_average(x::QAdd)::QAdd = x

function undo_average(eq::Symbolics.Equation)
    lhs = undo_average(eq.lhs)
    rhs = undo_average(eq.rhs)
    return lhs => rhs
end

# --- Metadata helpers ---

"""
    has_sum_metadata(x) -> Bool

Check whether a symbolic expression carries summation index metadata.
"""
has_sum_metadata(::Any) = false
function has_sum_metadata(x::SymbolicUtils.BasicSymbolic)
    return SymbolicUtils.hasmetadata(x, SumIndices)
end

"""
    get_sum_indices(x::BasicSymbolic) -> Vector{Index}

Retrieve summation indices from a symbolic expression's metadata.
"""
function get_sum_indices(x::SymbolicUtils.BasicSymbolic)
    return SymbolicUtils.getmetadata(x, SumIndices)
end

"""
    get_sum_non_equal(x::BasicSymbolic) -> Vector{Tuple{Index,Index}}

Retrieve non-equal index pairs from a symbolic expression's metadata.
"""
function get_sum_non_equal(x::SymbolicUtils.BasicSymbolic)
    return SymbolicUtils.getmetadata(x, SumNonEqual)
end

# --- get_indices for BasicSymbolic ---

function get_indices(x::SymbolicUtils.BasicSymbolic)
    if SymbolicUtils.isconst(x)
        return get_indices(x.val)
    end
    if SymbolicUtils.iscall(x)
        f = SymbolicUtils.operation(x)
        if f isa AvgFunc
            arg = SymbolicUtils.arguments(x)[1]
            inner = SymbolicUtils.isconst(arg) ? arg.val : arg
            return get_indices(inner)
        end
        inds = Index[]
        for arg in SymbolicUtils.arguments(x)
            for idx in get_indices(arg)
                idx ∉ inds && push!(inds, idx)
            end
        end
        return inds
    end
    return Index[]
end

# --- acts_on ---

"""
    acts_on(expr) -> Vector{Int}

Return sorted unique Hilbert space indices that an operator or averaged expression acts on.
Not intended for hot-path use — allocates a fresh `Vector{Int}` on each call.
"""
acts_on(op::QSym) = Int[op.space_index]

function acts_on(op::QAdd)
    aon = Int[]
    for (ops, _) in op.arguments
        for x in ops
            push!(aon, x.space_index)
        end
    end
    unique!(sort!(aon))
    return aon
end

function acts_on(s::SymbolicUtils.BasicSymbolic)
    if SymbolicUtils.isconst(s)
        return acts_on(s.val)
    end
    if SymbolicUtils.iscall(s)
        f = SymbolicUtils.operation(s)
        if f isa AvgFunc
            return acts_on(SymbolicUtils.arguments(s)[1])
        else
            aon = Int[]
            for arg in SymbolicUtils.arguments(s)
                append!(aon, acts_on(arg))
            end
            unique!(sort!(aon))
            return aon
        end
    else
        return Int[]
    end
end

acts_on(::Number) = Int[]
