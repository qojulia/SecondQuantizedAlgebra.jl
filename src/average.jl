"""
    AvgSym

Marker `symtype` for averaged operator expressions: `BasicSymbolic{SymReal}`
`Term` nodes with `symtype(x) === AvgSym`. Not public — use [`average`](@ref)
and [`is_average`](@ref).
"""
struct AvgSym <: Number end

"""Metadata key for summation indices on averaged expressions."""
struct SumIndices end

"""Metadata key for non-equal index pairs on averaged expressions."""
struct SumNonEqual end

"""
    AvgFunc

Singleton callable used as the `operation` of average `Term` nodes; defining
`show_call` on it gives the `⟨…⟩` display without type piracy.
"""
struct AvgFunc end
const sym_average = AvgFunc()

Base.nameof(::AvgFunc) = :avg
Base.show(io::IO, ::AvgFunc) = print(io, "avg")

function SymbolicUtils.show_call(io::IO, ::AvgFunc, x::SymbolicUtils.BasicSymbolic; kw...)
    print(io, "⟨")
    for (i, arg) in enumerate(SymbolicUtils.arguments(x))
        i > 1 && print(io, ", ")
        print(io, arg)
    end
    return print(io, "⟩")
end

"""
    is_average(x) -> Bool

Whether `x` is a symbolic average object created by [`average`](@ref).

```jldoctest
julia> h = FockSpace(:f);

julia> @qnumbers a::Destroy(h);

julia> is_average(average(a)), is_average(a)
(true, false)
```
"""
is_average(::Any) = false
is_average(x::SymbolicUtils.BasicSymbolic) = SymbolicUtils.iscall(x) && SymbolicUtils.symtype(x) === AvgSym
is_average(x::Num) = is_average(SymbolicUtils.unwrap(x))

_average(op::QField) = SymbolicUtils.Term{SymbolicUtils.SymReal}(sym_average, QField[op]; type = AvgSym)

"""
    average(expr) -> BasicSymbolic | Number

Build the symbolic average ``\\langle \\mathrm{expr} \\rangle``. Distributes over
sums, pulls c-number prefactors out, leaves scalars unchanged. Displayed as `⟨…⟩`.

```jldoctest
julia> h = FockSpace(:f);

julia> @qnumbers a::Destroy(h);

julia> avg = average(a' * a);

julia> is_average(avg)
true

julia> avg + 1
1 + ⟨a' * a⟩
```

See also [`undo_average`](@ref), [`numeric_average`](@ref).
"""
function average end

average(op::QSym) = _average(op)
average(x::Number) = x
average(x::SymbolicUtils.BasicSymbolic) = x
average(x::Num) = average(SymbolicUtils.unwrap(x))

function average(op::QAdd)
    # Accumulate into a `Num` (or `BasicSymbolic{SymReal}` after the first
    # multiplication-by-`avg`). Never let `result` become a `Complex{Num}`:
    # `SymbolicUtils.unwrap(::Complex{<:Num})` would materialise the imaginary
    # piece as a literal `complex(re, im)` symbolic call, which is opaque to
    # `simplify` / `expand`. We bring in `im` via `Symbolics.IM`
    # (the BasicSymbolic{SymReal} sym for `im`) so the chain stays symbolic.
    result = Num(0)
    shared = isempty(op.indices) ? nothing : copy(op.indices)
    for (term, c) in op.arguments
        r, i = real(c), imag(c)
        if isempty(term.ops)
            iszero(r) || (result += r)
            iszero(i) || (result += i * Symbolics.IM)
            continue
        end
        term_uses_sum = shared !== nothing &&
            any(idx -> _depends_on_index_term(c, term.ops, idx), shared)
        inner = (length(term.ops) == 1 && !term_uses_sum) ?
            only(term.ops) : _single_qadd(_CNUM_ONE, term.ops, term.ne)
        avg = _average(inner)
        if term_uses_sum
            avg = SymbolicUtils.setmetadata(avg, SumIndices, shared)
            avg = SymbolicUtils.setmetadata(avg, SumNonEqual, _copy_ne(term.ne))
        end
        iszero(r) || (result += r * avg)
        iszero(i) || (result += i * Symbolics.IM * avg)
    end
    return SymbolicUtils.unwrap(result)
end

# Uniform-return wrappers (all return QAdd).
_to_qadd(x::QAdd) = x
_to_qadd(x::QSym) = _single_qadd(_CNUM_ONE, QSym[x])
_to_qadd(x::SymbolicUtils.BasicSymbolic) = _single_qadd(_to_cnum(x), QSym[])

# Metadata is stored as `ImmutableDict{DataType, Any}`; the isa-narrow seals
# each result to its concrete type without a return annotation.
function _restore_sum_metadata_indices(x::SymbolicUtils.BasicSymbolic)
    v = SymbolicUtils.getmetadata(x, SumIndices)
    v isa Vector{Index} && return v
    throw(ArgumentError("SumIndices metadata has unexpected type $(typeof(v))"))
end
function _restore_sum_metadata_ne(x::SymbolicUtils.BasicSymbolic)
    SymbolicUtils.hasmetadata(x, SumNonEqual) || return _EMPTY_NE
    v = SymbolicUtils.getmetadata(x, SumNonEqual)
    v isa Vector{NonEqualPair} && return v
    throw(ArgumentError("SumNonEqual metadata has unexpected type $(typeof(v))"))
end

function _restore_sum_metadata(result::QAdd, x::SymbolicUtils.BasicSymbolic)
    SymbolicUtils.hasmetadata(x, SumIndices) || return result
    indices = _restore_sum_metadata_indices(x)
    stored_ne = _restore_sum_metadata_ne(x)
    new_args = QTermDict()
    for (term, c) in result.arguments
        _addto!(new_args, term.ops, c, _merge_ne(term.ne, stored_ne))
    end
    return QAdd(new_args, indices)
end

function _fold_qadds(op::F, args::Vector{QAdd}, empty::QAdd) where {F}
    isempty(args) && return empty
    result = first(args)
    for i in 2:length(args)
        result = op(result, args[i])
    end
    return result
end

"""
    undo_average(expr) -> QAdd

Recursively strip symbolic averages and return the underlying operator
expression. Summation metadata is restored. Also accepts a `Symbolics.Equation`,
returning a `Pair{QAdd, QAdd}`.

```jldoctest
julia> h = FockSpace(:f);

julia> @qnumbers a::Destroy(h);

julia> undo_average(average(a' * a)) == a' * a
true
```
"""
function undo_average(x::SymbolicUtils.BasicSymbolic)
    SymbolicUtils.iscall(x) || return _to_qadd(x)
    f = SymbolicUtils.operation(x)
    if f isa AvgFunc
        arg = SymbolicUtils.arguments(x)[1]
        inner = if SymbolicUtils.isconst(arg) && (arg.val isa QField || arg.val isa Number)
            arg.val
        else
            arg
        end
        return _restore_sum_metadata(_to_qadd(inner), x)
    end
    if f === (+) || f === (*)
        args = QAdd[undo_average(a) for a in SymbolicUtils.arguments(x)]
        folded = f === (+) ?
            _fold_qadds(+, args, _zero_qadd()) :
            _fold_qadds(*, args, _single_qadd(_CNUM_ONE, _EMPTY_OPS))
        return _restore_sum_metadata(folded, x)
    end
    return _to_qadd(x)
end

undo_average(x::Number) = _single_qadd(_to_cnum(x), QSym[])
undo_average(x::Num) = undo_average(SymbolicUtils.unwrap(x))
undo_average(x::QSym) = _single_qadd(_CNUM_ONE, QSym[x])
undo_average(x::QAdd) = x
undo_average(eq::Symbolics.Equation) = undo_average(eq.lhs) => undo_average(eq.rhs)

"""
    has_sum_metadata(x) -> Bool

Whether `x` is a `BasicSymbolic` node carrying summation index metadata set by
[`average`](@ref) on indexed expressions.

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
has_sum_metadata(::Any) = false
has_sum_metadata(x::SymbolicUtils.BasicSymbolic) = SymbolicUtils.hasmetadata(x, SumIndices)
has_sum_metadata(x::Num) = has_sum_metadata(SymbolicUtils.unwrap(x))

"""
    get_sum_indices(x::BasicSymbolic) -> Vector{Index}

Summation indices stored as metadata on `x`. Only valid when
[`has_sum_metadata(x)`](@ref has_sum_metadata) is `true`.

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
get_sum_indices(x::SymbolicUtils.BasicSymbolic) = _restore_sum_metadata_indices(x)
get_sum_indices(x::Num) = get_sum_indices(SymbolicUtils.unwrap(x))

"""
    get_sum_non_equal(x::BasicSymbolic) -> Vector{Tuple{Index, Index}}

Pairwise index-inequality constraints stored on an averaged term. An empty
vector means no constraints. Only valid when
[`has_sum_metadata(x)`](@ref has_sum_metadata) is `true`.

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
get_sum_non_equal(x::SymbolicUtils.BasicSymbolic) = _restore_sum_metadata_ne(x)
get_sum_non_equal(x::Num) = get_sum_non_equal(SymbolicUtils.unwrap(x))

# Seals: recursive calls go through `Any`-typed inputs; isa-narrow restores
# the typed accumulator without a return annotation.
_idx_seal(v) = v isa Vector{Index} ? v : Index[]
_aon_seal(v) = v isa Vector{Int} ? v : Int[]

function get_indices(x::SymbolicUtils.BasicSymbolic)
    SymbolicUtils.isconst(x) && return _idx_seal(get_indices(x.val))
    SymbolicUtils.iscall(x) || return Index[]
    f = SymbolicUtils.operation(x)
    if f isa AvgFunc
        arg = SymbolicUtils.arguments(x)[1]
        inner = SymbolicUtils.isconst(arg) ? arg.val : arg
        return _idx_seal(get_indices(inner))
    end
    inds = Index[]
    for arg in SymbolicUtils.arguments(x), idx in _idx_seal(get_indices(arg))
        idx ∉ inds && push!(inds, idx)
    end
    return inds
end

"""
    acts_on(expr) -> Vector{Int}

Sorted unique `space_index` values that `expr` acts on. Works on [`QSym`](@ref),
[`QAdd`](@ref), averaged `BasicSymbolic` expressions, and `Number`s (`Int[]`).

```jldoctest
julia> h = FockSpace(:a) ⊗ NLevelSpace(:b, 2);

julia> @qnumbers a::Destroy(h, 1) σ::Transition(h, 1, 2, 2);

julia> acts_on(a' * a), acts_on(a' * σ)
([1], [1, 2])
```
"""
function acts_on end

acts_on(op::QSym) = Int[op.space_index]
acts_on(::Number) = Int[]
acts_on(x::Num) = acts_on(SymbolicUtils.unwrap(x))

function acts_on(t::QTerm)
    aon = Int[x.space_index for x in t.ops]
    unique!(sort!(aon))
    return aon
end

function acts_on(op::QAdd)
    aon = Int[]
    for term in keys(op.arguments)
        append!(aon, acts_on(term))
    end
    unique!(sort!(aon))
    return aon
end

function acts_on(s::SymbolicUtils.BasicSymbolic)
    SymbolicUtils.isconst(s) && return _aon_seal(acts_on(s.val))
    SymbolicUtils.iscall(s) || return Int[]
    f = SymbolicUtils.operation(s)
    f isa AvgFunc && return _aon_seal(acts_on(SymbolicUtils.arguments(s)[1]))
    aon = Int[]
    for arg in SymbolicUtils.arguments(s)
        append!(aon, _aon_seal(acts_on(arg)))
    end
    unique!(sort!(aon))
    return aon
end
