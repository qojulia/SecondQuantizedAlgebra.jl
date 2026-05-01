"""
    AvgSym

Marker `symtype` for averaged operator expressions in the SymbolicUtils tree.

Average nodes are `BasicSymbolic{SymReal}` `Term` nodes with `symtype(x) === AvgSym`.
Not part of the public API — use [`average`](@ref) and [`is_average`](@ref) instead.
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

Return `true` if `x` is a symbolic average node (created by [`average`](@ref)).
"""
is_average(::QField) = false
is_average(::Number) = false
function is_average(x::SymbolicUtils.BasicSymbolic)
    return SymbolicUtils.iscall(x) && SymbolicUtils.symtype(x) === AvgSym
end
is_average(::Any) = false
is_average(x::Num) = is_average(SymbolicUtils.unwrap(x))

function _average(op::QField)
    return SymbolicUtils.Term{SymbolicUtils.SymReal}(sym_average, Any[op]; type = AvgSym)
end

"""
    average(expr) -> BasicSymbolic | Number

Compute the symbolic average ``\\langle \\mathrm{expr} \\rangle`` of a quantum operator expression.

Returns a `SymbolicUtils.BasicSymbolic` node (not a `Symbolics.Num` wrapper).
This ensures a consistent return type regardless of whether the input has symbolic
prefactors.

Average nodes participate in standard symbolic arithmetic (`+`, `*`, `^`, etc.)
and are displayed as `⟨...⟩`. Scalars pass through unchanged.

**Linearity**: distributes over [`QAdd`](@ref) sums and pulls c-number prefactors
outside the average. Summation metadata (indices, non-equal constraints) is preserved
as SymbolicUtils metadata on the result.

# Examples
```julia
h = FockSpace(:f)
@qnumbers a::Destroy(h)
avg = average(a' * a)     # ⟨a†a⟩
avg + 1                    # ⟨a†a⟩ + 1  (symbolic arithmetic)
```

See also [`undo_average`](@ref), [`is_average`](@ref), [`numeric_average`](@ref).
"""
function average end

average(op::QSym) = _average(op)
average(x::Number) = x
average(x::SymbolicUtils.BasicSymbolic) = x
average(x::Num) = average(SymbolicUtils.unwrap(x))

function average(op::QAdd)
    # Use Num arithmetic internally for convenience, unwrap at the end.
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
    out = SymbolicUtils.unwrap(result)
    if !isempty(op.indices)
        out = SymbolicUtils.setmetadata(out, SumIndices, op.indices)
        out = SymbolicUtils.setmetadata(out, SumNonEqual, op.non_equal)
    end
    return out
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

Recursively strip ``\\langle \\cdots \\rangle`` from a symbolic expression,
recovering the underlying operator expressions.

Always returns [`QAdd`](@ref) for type stability:
- Scalars become single-term `QAdd`s with an empty operator sequence
- Single operators become single-term `QAdd`s with unit prefactor
- Sums/products of averages are reconstructed via `+` and `*` on `QAdd`s

Summation metadata (indices, non-equal constraints) is restored from the
SymbolicUtils metadata onto the resulting `QAdd`.

Also accepts `Symbolics.Equation`, returning a `Pair{QAdd, QAdd}`.

See also [`average`](@ref), [`is_average`](@ref).
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

Return `true` if `x` is a `BasicSymbolic` node carrying summation index metadata
(set by [`average`](@ref) when averaging indexed expressions).

See also [`get_sum_indices`](@ref), [`get_sum_non_equal`](@ref).
"""
has_sum_metadata(::Number) = false
has_sum_metadata(::QField) = false
function has_sum_metadata(x::SymbolicUtils.BasicSymbolic)
    return SymbolicUtils.hasmetadata(x, SumIndices)
end
has_sum_metadata(x::Num) = has_sum_metadata(SymbolicUtils.unwrap(x))

"""
    get_sum_indices(x::BasicSymbolic) -> Vector{Index}

Retrieve summation indices stored as metadata on a symbolic expression `x`.

Only valid when [`has_sum_metadata(x)`](@ref has_sum_metadata) is `true`.

See also [`get_sum_non_equal`](@ref), [`has_sum_metadata`](@ref).
"""
function get_sum_indices(x::SymbolicUtils.BasicSymbolic)
    return SymbolicUtils.getmetadata(x, SumIndices)
end
get_sum_indices(x::Num) = get_sum_indices(SymbolicUtils.unwrap(x))

"""
    get_sum_non_equal(x::BasicSymbolic) -> Vector{Tuple{Index,Index}}

Retrieve pairwise index inequality constraints stored as metadata on `x`.

Only valid when [`has_sum_metadata(x)`](@ref has_sum_metadata) is `true`.

See also [`get_sum_indices`](@ref), [`has_sum_metadata`](@ref).
"""
function get_sum_non_equal(x::SymbolicUtils.BasicSymbolic)
    return SymbolicUtils.getmetadata(x, SumNonEqual)
end
get_sum_non_equal(x::Num) = get_sum_non_equal(SymbolicUtils.unwrap(x))

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

Return the sorted unique `space_index` values that `expr` acts on.

Works on [`QSym`](@ref), [`QAdd`](@ref), averaged `BasicSymbolic` expressions,
and `Number`s (returns `Int[]`).

# Examples
```julia
h = FockSpace(:a) ⊗ NLevelSpace(:b, 2)
@qnumbers a::Destroy(h, 1) σ::Transition(h, 1, 2, 2)
acts_on(a' * a)         # [1]
acts_on(a' * σ)         # [1, 2]
```
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
acts_on(x::Num) = acts_on(SymbolicUtils.unwrap(x))
