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

Return whether `x` is a symbolic average object created by [`average`](@ref).

# Examples

```jldoctest
julia> h = FockSpace(:f);

julia> @qnumbers a::Destroy(h);

julia> is_average(average(a)), is_average(a)
(true, false)
```

See also [`average`](@ref), [`undo_average`](@ref).
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

Build the symbolic average ``\\langle \\mathrm{expr} \\rangle`` of an operator expression.

Distributes over sums and pulls c-number prefactors outside the brackets. Scalars
pass through unchanged. The result participates in standard symbolic arithmetic
(`+`, `*`, `^`) and is displayed as `⟨...⟩`.

# Examples

```jldoctest
julia> h = FockSpace(:f);

julia> @qnumbers a::Destroy(h);

julia> avg = average(a' * a);

julia> is_average(avg)
true

julia> avg + 1
1 + ⟨a' * a⟩
```

See also [`undo_average`](@ref), [`is_average`](@ref), [`numeric_average`](@ref).
"""
function average end

average(op::QSym) = _average(op)
average(x::Number) = x
average(x::SymbolicUtils.BasicSymbolic) = x
average(x::Num) = average(SymbolicUtils.unwrap(x))

function _average_operand(
        ops::Vector{QSym}, ne::Vector{NonEqualPair}, indices::Vector{Index}
    )
    if length(ops) == 1 && isempty(ne) && isempty(indices)
        return only(ops)
    end
    return _single_qadd(_CNUM_ONE, ops, ne)
end

function _average_with_metadata(
        inner::QField, shared_indices::Union{Nothing, Vector{Index}},
        ne::Vector{NonEqualPair}
    )
    avg = _average(inner)
    if shared_indices !== nothing
        avg = SymbolicUtils.setmetadata(avg, SumIndices, shared_indices)
        avg = SymbolicUtils.setmetadata(avg, SumNonEqual, _copy_ne(ne))
    end
    return avg
end

function average(op::QAdd)
    result = Num(0)
    shared_indices = isempty(op.indices) ? nothing : copy(op.indices)
    for (term, c) in op.arguments
        if isempty(term.ops)
            result += c
        else
            inner = _average_operand(term.ops, term.ne, op.indices)
            avg = _average_with_metadata(inner, shared_indices, term.ne)
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
    return SymbolicUtils.unwrap(result)
end

"""Wrap any value as a QAdd, ensuring uniform return type."""
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
        stored_ne = SymbolicUtils.hasmetadata(x, SumNonEqual) ?
            SymbolicUtils.getmetadata(x, SumNonEqual) : _EMPTY_NE
        new_args = QTermDict()
        for (term, c) in result.arguments
            _addto!(new_args, term.ops, c, _merge_ne(term.ne, stored_ne))
        end
        return QAdd(new_args, indices)
    end
    return result
end

"""
    undo_average(expr) -> QAdd

Recursively strip symbolic averages ``\\langle \\cdots \\rangle`` and recover
the underlying operator expression as a [`QAdd`](@ref). Summation metadata
(indices, non-equal constraints) is restored on the result. Also accepts a
`Symbolics.Equation`, returning a `Pair{QAdd, QAdd}`.

# Examples

```jldoctest
julia> h = FockSpace(:f);

julia> @qnumbers a::Destroy(h);

julia> undo_average(average(a' * a)) == a' * a
true
```

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

"""
    has_sum_metadata(x) -> Bool

Return `true` if `x` is a `BasicSymbolic` node carrying summation index metadata
(set by [`average`](@ref) when averaging indexed expressions).

# Examples

```jldoctest
julia> h = FockSpace(:site);

julia> @qnumbers a::Destroy(h);

julia> i = Index(h, :i, 3, h);

julia> x = average(Σ(IndexedOperator(a', i) * IndexedOperator(a, i), i));

julia> SecondQuantizedAlgebra.has_sum_metadata(x)
true
```

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

# Examples

```jldoctest
julia> h = FockSpace(:site);

julia> @qnumbers a::Destroy(h);

julia> i = Index(h, :i, 3, h);

julia> x = average(Σ(IndexedOperator(a', i) * IndexedOperator(a, i), i));

julia> SecondQuantizedAlgebra.get_sum_indices(x) == [i]
true
```

See also [`get_sum_non_equal`](@ref), [`has_sum_metadata`](@ref).
"""
function get_sum_indices(x::SymbolicUtils.BasicSymbolic)
    return SymbolicUtils.getmetadata(x, SumIndices)
end
get_sum_indices(x::Num) = get_sum_indices(SymbolicUtils.unwrap(x))

"""
    get_sum_non_equal(x::BasicSymbolic) -> Vector{Tuple{Index,Index}}

Retrieve the pairwise index inequality constraints stored on an averaged term.
An empty vector means no constraints.

Only valid when [`has_sum_metadata(x)`](@ref has_sum_metadata) is `true`.

# Examples

```jldoctest
julia> h = FockSpace(:site);

julia> @qnumbers a::Destroy(h);

julia> i = Index(h, :i, 3, h);

julia> x = average(Σ(IndexedOperator(a', i) * IndexedOperator(a, i), i));

julia> isempty(SecondQuantizedAlgebra.get_sum_non_equal(x))
true
```

See also [`get_sum_indices`](@ref), [`has_sum_metadata`](@ref).
"""
function get_sum_non_equal(x::SymbolicUtils.BasicSymbolic)
    return SymbolicUtils.getmetadata(x, SumNonEqual)
end
get_sum_non_equal(x::Num) = get_sum_non_equal(SymbolicUtils.unwrap(x))

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

"""
    acts_on(expr) -> Vector{Int}

Return the sorted unique `space_index` values that `expr` acts on.

Works on [`QSym`](@ref), [`QAdd`](@ref), averaged `BasicSymbolic` expressions,
and `Number`s (returns `Int[]`).

# Examples

```jldoctest
julia> h = FockSpace(:a) ⊗ NLevelSpace(:b, 2);

julia> @qnumbers a::Destroy(h, 1) σ::Transition(h, 1, 2, 2);

julia> acts_on(a' * a)
1-element Vector{Int64}:
 1

julia> acts_on(a' * σ)
2-element Vector{Int64}:
 1
 2
```
"""
acts_on(op::QSym) = Int[op.space_index]

function acts_on(op::QAdd)
    aon = Int[]
    for term in keys(op.arguments)
        for x in term.ops
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
