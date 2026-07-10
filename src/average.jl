"""
    AverageOperator

Metadata key carrying the operator of a lifted (time-dependent) average variable.
The lifted form is a `Number`-symtype `var(iv)` produced by [`make_time_dependent`](@ref); the
operator moves to metadata because the single call argument is the independent
variable. Not public, use [`average`](@ref), [`make_time_dependent`](@ref) and [`is_average`](@ref).
"""
struct AverageOperator end

"""
    AvgFunc

Singleton callable used as the `operation` of average `Term` nodes; defining
`show_call` on it gives the `⟨…⟩` display without type piracy.
"""
struct AvgFunc end
const sym_average = AvgFunc()

Base.nameof(::AvgFunc) = :avg
Base.show(io::IO, ::AvgFunc) = print(io, "avg")
Base.isless(::AvgFunc, ::AvgFunc) = false
SymbolicUtils.promote_symtype(::AvgFunc, Ts...) = Number

function SymbolicUtils.show_call(io::IO, ::AvgFunc, x::SymbolicUtils.BasicSymbolic; kw...)
    print(io, "⟨")
    for (i, arg) in enumerate(SymbolicUtils.arguments(x))
        i > 1 && print(io, ", ")
        print(io, arg)
    end
    return print(io, "⟩")
end

"""
    SumScope

Opaque carrier for an indexed-sum's summation scope (`indices::Vector{Index}`
and `ne::Vector{NonEqualPair}`). Stored as a single scalar argument of a
[`SumFunc`](@ref) node: SymbolicUtils array-ifies a bare `Vector` argument
(losing the `Index` objects), but a scalar struct is kept opaque and
recoverable, exactly like the `QField` argument of an average node. Custom
`==`/`isequal`/`hash` compare it by value so two scopes that differ wrongly
cancel only when they are genuinely equal.
"""
struct SumScope
    indices::Vector{Index}
    ne::Vector{NonEqualPair}
end
Base.:(==)(a::SumScope, b::SumScope) = a.indices == b.indices && a.ne == b.ne
Base.isequal(a::SumScope, b::SumScope) = isequal(a.indices, b.indices) && isequal(a.ne, b.ne)
Base.hash(s::SumScope, h::UInt) = hash(s.ne, hash(s.indices, hash(:SumScope, h)))

"""
    SumFunc

Singleton callable used as the `operation` of moment-layer indexed-sum `Term`
nodes (cf. [`AvgFunc`](@ref)). A node is `sym_sum(body, scope)`: the summed
`body` is an arbitrary averaged expression and the [`SumScope`](@ref) rides as a
plain `Term` argument, so it participates in `isequal`/`hash` (metadata would be
ignored, letting two differently-scoped sums over the same body wrongly cancel).
Because the operation is neither `+` nor `*`, `Add`/`Mul` keep the node opaque,
so the body and scope survive SymbolicUtils canonicalization.
"""
struct SumFunc end
const sym_sum = SumFunc()

Base.nameof(::SumFunc) = :sum
Base.show(io::IO, ::SumFunc) = print(io, "sum")
Base.isless(::SumFunc, ::SumFunc) = false
SymbolicUtils.promote_symtype(::SumFunc, Ts...) = Number

_indexed_sum(body, indices::Vector{Index}, ne::Vector{NonEqualPair}) =
    SymbolicUtils.Term{SymbolicUtils.SymReal}(sym_sum, Any[body, SumScope(indices, ne)]; type = Number)

_scope_of(s::SumScope) = s
_scope_of(s::SymbolicUtils.BasicSymbolic) = s.val::SumScope
_sum_body(x::SymbolicUtils.BasicSymbolic) = SymbolicUtils.arguments(x)[1]
_sum_scope(x::SymbolicUtils.BasicSymbolic) = _scope_of(SymbolicUtils.arguments(x)[2])
_sum_indices(x::SymbolicUtils.BasicSymbolic) = _sum_scope(x).indices
_sum_ne(x::SymbolicUtils.BasicSymbolic) = _sum_scope(x).ne

"""
    is_indexed_sum(x) -> Bool

Whether `x` is a moment-layer indexed-sum node created by [`average`](@ref) on a
summed [`QAdd`](@ref). The summation scope is recovered with
[`get_sum_indices`](@ref) and [`get_sum_non_equal`](@ref).
"""
is_indexed_sum(::Any) = false
is_indexed_sum(x::SymbolicUtils.BasicSymbolic) =
    SymbolicUtils.iscall(x) && SymbolicUtils.operation(x) isa SumFunc
is_indexed_sum(x::Num) = is_indexed_sum(SymbolicUtils.unwrap(x))

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
function is_average(x::SymbolicUtils.BasicSymbolic)
    (SymbolicUtils.iscall(x) && SymbolicUtils.operation(x) === sym_average) && return true
    return SymbolicUtils.hasmetadata(x, AverageOperator)
end
is_average(x::Num) = is_average(SymbolicUtils.unwrap(x))

_average(op::QField) = SymbolicUtils.Term{SymbolicUtils.SymReal}(sym_average, QField[op]; type = Number)

"""
    make_time_dependent(expr, iv) -> expr

Lift every iv-free average leaf in `expr` into a time-dependent `Number` variable
`name(iv)` carrying the operator in `AverageOperator` metadata (and the
`VariableSource` set by `@variables`, so it reads as a ModelingToolkit unknown).
Non-average structure is rebuilt only where a child changed. The lifted node is a
leaf: the walk does not descend into `iv`. The default `name` is uniqueness-only;
downstream code may rename for display without changing identity.
"""
make_time_dependent(x::Num, iv) = Symbolics.wrap(make_time_dependent(SymbolicUtils.unwrap(x), iv))
function make_time_dependent(x, iv)
    x isa SymbolicUtils.BasicSymbolic || return x
    if SymbolicUtils.iscall(x) && SymbolicUtils.operation(x) === sym_average
        return _lift_average(x, iv)
    end
    if is_indexed_sum(x)
        return _indexed_sum(make_time_dependent(_sum_body(x), iv), _sum_indices(x), _sum_ne(x))
    end
    SymbolicUtils.iscall(x) || return x
    args = SymbolicUtils.arguments(x)
    new_args = Vector{Any}(undef, length(args))
    changed = false
    for k in eachindex(args)
        na = make_time_dependent(args[k], iv)
        new_args[k] = na
        changed |= na !== args[k]
    end
    changed || return x   # nothing rewritten: keep the node (and its identity) intact
    return TermInterface.maketerm(typeof(x), SymbolicUtils.operation(x), new_args, TermInterface.metadata(x))
end

function _lift_average(leaf, iv)
    op = undo_average(leaf)
    nm = _avg_var_name(op)
    ivw = iv isa SymbolicUtils.BasicSymbolic ? Symbolics.wrap(iv) : iv
    v = SymbolicUtils.unwrap(first(@variables ($nm(ivw))::Number))
    return SymbolicUtils.setmetadata(v, AverageOperator, op)
end

# Deterministic, structurally injective name for a lifted variable. Identity rests on this
# name: Symbolics `isequal`/`hash` ignore the `AverageOperator` metadata, so the name must
# be 1:1 with the operator. `hash(op)` is a digest (collisions silently merge two distinct
# unknowns); `string(op)` alone drops the subspace (`a` on different spaces prints alike),
# so the per-operator token prepends `space_index` to the faithful operator rendering.
function _avg_var_name(op::QAdd)
    io = IOBuffer()
    print(io, "_avg")
    for (term, c) in op.arguments
        print(io, '_', c)
        for o in term.ops
            print(io, '_', o.space_index, o)
        end
        for (p, q) in term.ne
            print(io, "_ne", index_name(p), index_name(q))
        end
    end
    return Symbol(String(take!(io)))
end

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
    BS = SymbolicUtils.BasicSymbolic{SymbolicUtils.SymReal}
    shared = isempty(op.indices) ? Index[] : op.indices
    GroupKey = Tuple{Vector{Index}, Vector{NonEqualPair}}
    flat_bodies = BS[]                   # ungrouped contributions
    group_keys = GroupKey[]
    group_bodies = Vector{BS}[]          # one buffer per summation scope
    group_slot = Dict{GroupKey, Int}()   # scope -> index into group_bodies (O(1) lookup)
    for (term, c) in op.arguments
        r, i = _realimag(c)
        used = Index[idx for idx in shared if _depends_on_index_term(c, term.ops, idx)]
        contrib = Num(0)
        if isempty(term.ops)
            iszero(r) || (contrib += r)
            iszero(i) || (contrib += i * Symbolics.IM)
        else
            inner = length(term.ops) == 1 ? only(term.ops) :
                _single_qadd(_CNUM_ONE, term.ops, term.ne)
            avg = _average(inner)
            iszero(r) || (contrib += r * avg)
            iszero(i) || (contrib += i * Symbolics.IM * avg)
        end
        body = SymbolicUtils.unwrap(contrib)
        if isempty(used)
            push!(flat_bodies, body)
        else
            key = (used, term.ne)
            slot = get(group_slot, key, 0)
            if slot == 0
                push!(group_keys, (used, _copy_ne(term.ne)))
                push!(group_bodies, BS[body])
                group_slot[key] = length(group_bodies)
            else
                push!(group_bodies[slot], body)
            end
        end
    end
    for (gk, bodies) in zip(group_keys, group_bodies)
        push!(flat_bodies, _indexed_sum(add_worker(SymbolicUtils.SymReal, bodies), gk[1], gk[2]))
    end
    isempty(flat_bodies) && return SymbolicUtils.unwrap(Num(0))
    return add_worker(SymbolicUtils.SymReal, flat_bodies)::BS
end

# Uniform-return wrappers (all return QAdd).
_to_qadd(x::QAdd) = x
_to_qadd(x::QSym) = _single_qadd(_CNUM_ONE, Op[x])
_to_qadd(x::SymbolicUtils.BasicSymbolic) = _single_qadd(_to_cnum(x), Op[])

function _rebuild_indexed_sum(inner::QAdd, indices::Vector{Index}, ne::Vector{NonEqualPair})
    isempty(ne) && return QAdd(inner.arguments, indices)
    new_args = QTermDict()
    for (term, c) in inner.arguments
        _addto!(new_args, term.ops, c, _merge_ne(term.ne, ne))
    end
    return QAdd(new_args, indices)
end

# The `+` fold is O(n²) by repeated `Base.:+` (each copies the growing dict); the
# in-place accumulator behind `sum` makes it O(n) and is byte-identical. The `*`
# fold stays as repeated `*` (the product path has no in-place win; see devdocs).
_fold_qadds(::typeof(+), args::Vector{QAdd}, empty::QAdd) = isempty(args) ? empty : sum(args)
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
    if SymbolicUtils.hasmetadata(x, AverageOperator)
        return _to_qadd(SymbolicUtils.getmetadata(x, AverageOperator))
    end
    if is_indexed_sum(x)
        return _rebuild_indexed_sum(undo_average(_sum_body(x)), _sum_indices(x), _sum_ne(x))
    end
    SymbolicUtils.iscall(x) || return _to_qadd(x)
    f = SymbolicUtils.operation(x)
    if f isa AvgFunc
        arg = SymbolicUtils.arguments(x)[1]
        inner = if SymbolicUtils.isconst(arg) && (arg.val isa QField || arg.val isa Number)
            arg.val
        else
            arg
        end
        return _to_qadd(inner)
    end
    if f === (+) || f === (*)
        args = QAdd[undo_average(a) for a in SymbolicUtils.arguments(x)]
        return f === (+) ?
            _fold_qadds(+, args, _zero_qadd()) :
            _fold_qadds(*, args, _single_qadd(_CNUM_ONE, _EMPTY_OPS))
    end
    return _to_qadd(x)
end

undo_average(x::Number) = _single_qadd(_to_cnum(x), Op[])
undo_average(x::Num) = undo_average(SymbolicUtils.unwrap(x))
undo_average(x::QSym) = _single_qadd(_CNUM_ONE, Op[x])
undo_average(x::QAdd) = x
undo_average(eq::Symbolics.Equation) = undo_average(eq.lhs) => undo_average(eq.rhs)

"""
    has_sum_metadata(x) -> Bool

Whether `x` is a moment-layer indexed-sum node (see [`is_indexed_sum`](@ref))
created by [`average`](@ref) on a summed expression. Retained as the public
predicate name; equivalent to `is_indexed_sum(x)`.

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
has_sum_metadata(x::SymbolicUtils.BasicSymbolic) = is_indexed_sum(x)
has_sum_metadata(x::Num) = has_sum_metadata(SymbolicUtils.unwrap(x))

"""
    get_sum_indices(x::BasicSymbolic) -> Vector{Index}

Summation indices carried by the indexed-sum node `x`. Only valid when
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
get_sum_indices(x::SymbolicUtils.BasicSymbolic) = _sum_indices(x)
get_sum_indices(x::Num) = get_sum_indices(SymbolicUtils.unwrap(x))

"""
    get_sum_non_equal(x::BasicSymbolic) -> Vector{Tuple{Index, Index}}

Pairwise index-inequality constraints carried by the indexed-sum node `x`. An
empty vector means no constraints. Only valid when
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
get_sum_non_equal(x::SymbolicUtils.BasicSymbolic) = _sum_ne(x)
get_sum_non_equal(x::Num) = get_sum_non_equal(SymbolicUtils.unwrap(x))

# Seals: recursive calls go through `Any`-typed inputs; isa-narrow restores
# the typed accumulator without a return annotation.
_idx_seal(v) = v isa Vector{Index} ? v : Index[]
_aon_seal(v) = v isa Vector{Int} ? v : Int[]

function get_indices(x::SymbolicUtils.BasicSymbolic)
    SymbolicUtils.isconst(x) && return _idx_seal(get_indices(x.val))
    SymbolicUtils.hasmetadata(x, AverageOperator) &&
        return _idx_seal(get_indices(SymbolicUtils.getmetadata(x, AverageOperator)))
    is_indexed_sum(x) && return _idx_seal(get_indices(_sum_body(x)))
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
    SymbolicUtils.hasmetadata(s, AverageOperator) &&
        return _aon_seal(acts_on(SymbolicUtils.getmetadata(s, AverageOperator)))
    is_indexed_sum(s) && return _aon_seal(acts_on(_sum_body(s)))
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
