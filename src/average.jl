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
SymbolicUtils.promote_symtype(::SumFunc, Ts...) = Number

indexed_sum(@nospecialize(body), indices::Vector{Index}, ne::Vector{NonEqualPair}) =
    SymbolicUtils.Term{SymbolicUtils.SymReal}(sym_sum, Any[body, SumScope(indices, ne)]; type = sum_symtype(body))

scope_of(s::SumScope) = s
function scope_of(s::SymbolicUtils.BasicSymbolic)
    scope = s.val
    scope isa SumScope && return scope
    throw(ArgumentError("expected an indexed-sum scope, got $s"))
end
sum_body(x::SymbolicUtils.BasicSymbolic) = SymbolicUtils.arguments(x)[1]
sum_scope(x::SymbolicUtils.BasicSymbolic) = scope_of(SymbolicUtils.arguments(x)[2])
sum_indices(x::SymbolicUtils.BasicSymbolic) = sum_scope(x).indices
sum_ne(x::SymbolicUtils.BasicSymbolic) = sum_scope(x).ne

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

avg_symtype(op::QField) = isequal(adjoint(op), op) ? Real : Number

function sum_symtype(@nospecialize(body))
    body isa RawExpression && return SymbolicUtils.symtype(body) <: Real ? Real : Number
    return body isa Real ? Real : Number
end

function avg_symtype_from_arg(@nospecialize(arg))
    op = arg isa RawExpression && SymbolicUtils.isconst(arg) ? arg.val : arg
    op isa Op && return avg_symtype(op)
    op isa QAdd && return avg_symtype(op)
    return Number
end

average_leaf(op::QField) = SymbolicUtils.Term{SymbolicUtils.SymReal}(sym_average, QField[op]; type = avg_symtype(op))

function TermInterface.maketerm(
        ::Type{SymbolicUtils.BasicSymbolic{SymbolicUtils.SymReal}}, f::AvgFunc, args, metadata;
        type = avg_symtype_from_arg(args[1])
    )
    return SymbolicUtils.Term{SymbolicUtils.SymReal}(f, collect(args); type = type, metadata = metadata)
end
function TermInterface.maketerm(
        ::Type{SymbolicUtils.BasicSymbolic{SymbolicUtils.SymReal}}, f::SumFunc, args, metadata;
        type = sum_symtype(args[1])
    )
    return SymbolicUtils.Term{SymbolicUtils.SymReal}(f, collect(args); type = type, metadata = metadata)
end

"""
    make_time_dependent(expr, iv) -> expr

Lift every iv-free average leaf in `expr` into a time-dependent variable `name(iv)`.

The variable is typed `Real` when the operator is Hermitian, else `Number`. It carries
the operator in `AverageOperator` metadata and the `VariableSource` set by `@variables`,
so it reads as a ModelingToolkit unknown. Non-average structure is rebuilt only where a child changed. The lifted node is a
leaf: the walk does not descend into `iv`. The default `name` is uniqueness-only;
downstream code may rename for display without changing identity.
"""
make_time_dependent(x::Num, iv) = Symbolics.wrap(make_time_dependent(SymbolicUtils.unwrap(x), iv))
make_time_dependent(x, iv) = x
function make_time_dependent(x::RawExpression, iv)
    SymbolicUtils.iscall(x) || return x
    f = SymbolicUtils.operation(x)
    f === sym_average && return lift_average(x, iv)
    if f === sym_sum
        return indexed_sum(lift_child(sum_body(x), iv), sum_indices(x), sum_ne(x))
    end
    args = SymbolicUtils.arguments(x)
    new_args = Vector{Any}(undef, length(args))
    changed = false
    for k in eachindex(args)
        na = lift_child(args[k], iv)
        new_args[k] = na
        changed |= na !== args[k]
    end
    changed || return x
    return Base.invoke(
        TermInterface.maketerm, Tuple{Type{RawExpression}, Any, Any, Any},
        RawExpression, f, new_args, TermInterface.metadata(x),
    )
end

lift_child(@nospecialize(a), iv) = a isa RawExpression ? make_time_dependent(a, iv) : a

function lift_average(leaf::RawExpression, iv)
    op = undo_average(leaf)
    nm = avg_var_name(op)
    raw_iv = SymbolicUtils.unwrap(iv)
    fn = if avg_symtype(op) === Real
        SymbolicUtils.Sym{SymbolicUtils.SymReal}(nm; type = SymbolicUtils.FnType{Tuple, Real, Nothing}, shape = UnitRange{Int}[])
    else
        SymbolicUtils.Sym{SymbolicUtils.SymReal}(nm; type = SymbolicUtils.FnType{Tuple, Number, Nothing}, shape = UnitRange{Int}[])
    end
    v = SymbolicUtils.setmetadata(fn(raw_iv), Symbolics.VariableSource, (:variables, nm))
    return SymbolicUtils.setmetadata(v, AverageOperator, op)
end

function avg_var_name(op::QAdd)
    io = IOBuffer()
    print(io, "avg")
    for (term, c) in op.arguments
        print(io, '_', c)
        for o in term.ops
            print(io, '_', o.space_index, o)
        end
        for (p, q) in term.ne
            print(io, "ne", index_name(p), index_name(q))
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

average(op::QSym) = average_leaf(op)
average(x::Number) = x
average(x::SymbolicUtils.BasicSymbolic) = x
average(x::Num) = average(SymbolicUtils.unwrap(x))
average(::QExpr) = throw(
    ArgumentError(
        "averages of formal operator functions are not defined; " *
            "call `taylor(expr, 0:n)` first and average the resulting `QAdd`",
    ),
)

function average(op::QAdd)
    BS = SymbolicUtils.BasicSymbolic{SymbolicUtils.SymReal}
    shared = isempty(op.indices) ? Index[] : op.indices
    GroupKey = Tuple{Vector{Index}, Vector{NonEqualPair}}
    flat_bodies = BS[]
    group_keys = GroupKey[]
    group_bodies = Vector{BS}[]
    group_slot = Dict{GroupKey, Int}()
    for (term, c) in op.arguments
        used = Index[idx for idx in shared if depends_on_scope_index(c, term.ops, term.ne, idx)]
        body = average_term(c, term)
        if isempty(used)
            push!(flat_bodies, body)
        else
            key = (used, term.ne)
            slot = get(group_slot, key, 0)
            if slot == 0
                push!(group_keys, (used, copy_ne(term.ne)))
                push!(group_bodies, BS[body])
                group_slot[key] = length(group_bodies)
            else
                push!(group_bodies[slot], body)
            end
        end
    end
    for (gk, bodies) in zip(group_keys, group_bodies)
        push!(flat_bodies, indexed_sum(add_worker(SymbolicUtils.SymReal, bodies), gk[1], gk[2]))
    end
    isempty(flat_bodies) && return SymbolicUtils.unwrap(Num(0))
    return add_worker(SymbolicUtils.SymReal, flat_bodies)
end

function average_term(c::CNum, term::QTerm)
    contrib = Num(0)
    if c.tail isa RawSymbolicCoeff && !is_call_of(c.tail.expr, complex)
        raw = Num(c.tail.expr)
        iszero(raw) && return SymbolicUtils.unwrap(contrib)
        contrib += isempty(term.ops) ? raw : raw * average_term_ops(term)
        return SymbolicUtils.unwrap(contrib)
    end
    r, i = realimag(c)
    if isempty(term.ops)
        iszero(r) || (contrib += r)
        iszero(i) || (contrib += i * Symbolics.IM)
    else
        avg = average_term_ops(term)
        iszero(r) || (contrib += r * avg)
        iszero(i) || (contrib += i * Symbolics.IM * avg)
    end
    return SymbolicUtils.unwrap(contrib)
end

function average_term_ops(term::QTerm)
    inner = length(term.ops) == 1 ? only(term.ops) : single_qadd(CNUM_ONE, term.ops, term.ne)
    return average_leaf(inner)
end

to_qadd(x::QAdd) = x
to_qadd(x::QSym) = single_qadd(CNUM_ONE, Op[x])
to_qadd(x::SymbolicUtils.BasicSymbolic) = single_qadd(to_cnum(x), Op[])

function rebuild_indexed_sum(inner::QAdd, indices::Vector{Index}, ne::Vector{NonEqualPair})
    isempty(ne) && return QAdd(inner.arguments, indices)
    new_args = QTermDict()
    for (term, c) in inner.arguments
        addto!(new_args, term.ops, c, merge_ne(term.ne, ne))
    end
    return QAdd(new_args, indices)
end

fold_qadds(::typeof(+), args::Vector{QAdd}, empty::QAdd) = isempty(args) ? empty : sum(args)
function fold_qadds(op::F, args::Vector{QAdd}, empty::QAdd) where {F}
    isempty(args) && return empty
    result = first(args)
    for i in 2:length(args)
        result = op(result, args[i])
    end
    return result
end

function contains_average(x)
    x isa Num && return contains_average(SymbolicUtils.unwrap(x))
    x isa SymbolicUtils.BasicSymbolic || return false
    (is_average(x) || is_indexed_sum(x)) && return true
    SymbolicUtils.iscall(x) || return false
    return any(contains_average, SymbolicUtils.arguments(x))
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
        return to_qadd(SymbolicUtils.getmetadata(x, AverageOperator))
    end
    if is_indexed_sum(x)
        return rebuild_indexed_sum(undo_average(sum_body(x)), sum_indices(x), sum_ne(x))
    end
    SymbolicUtils.iscall(x) || return to_qadd(x)
    f = SymbolicUtils.operation(x)
    if f isa AvgFunc
        arg = SymbolicUtils.arguments(x)[1]
        inner = if SymbolicUtils.isconst(arg) && (arg.val isa QField || arg.val isa Number)
            arg.val
        else
            arg
        end
        return to_qadd(inner)
    end
    if f === (/) && length(SymbolicUtils.arguments(x)) == 2
        numerator, denominator = SymbolicUtils.arguments(x)
        if !contains_average(denominator)
            return undo_average(numerator) * inv(to_cnum(denominator))
        end
    end
    if f === (+) || f === (*)
        args = QAdd[undo_average(a) for a in SymbolicUtils.arguments(x)]
        return f === (+) ?
            fold_qadds(+, args, zero_qadd()) :
            fold_qadds(*, args, single_qadd(CNUM_ONE, EMPTY_OPS))
    end
    return to_qadd(x)
end

undo_average(x::Number) = single_qadd(to_cnum(x), Op[])
undo_average(x::Num) = undo_average(SymbolicUtils.unwrap(x))
undo_average(x::QSym) = single_qadd(CNUM_ONE, Op[x])
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
get_sum_indices(x::SymbolicUtils.BasicSymbolic) = sum_indices(x)
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
get_sum_non_equal(x::SymbolicUtils.BasicSymbolic) = sum_ne(x)
get_sum_non_equal(x::Num) = get_sum_non_equal(SymbolicUtils.unwrap(x))

"""
    get_sum_body(x::BasicSymbolic)

Summed body carried by the indexed-sum node `x`: the averaged expression that the
node sums over its [`get_sum_indices`](@ref) (an index-dependent coefficient stays
*inside* this body). Only valid when [`has_sum_metadata(x)`](@ref has_sum_metadata)
is `true`. The inverse of [`indexed_sum`](@ref).

```jldoctest
julia> h = FockSpace(:site);

julia> @qnumbers a::Destroy(h);

julia> i = Index(h, :i, 3, h);

julia> x = average(Σ(IndexedOperator(a', i) * IndexedOperator(a, i), i));

julia> isequal(SecondQuantizedAlgebra.get_sum_body(x), average(IndexedOperator(a', i) * IndexedOperator(a, i)))
true
```

See also [`indexed_sum`](@ref), [`get_sum_indices`](@ref), [`get_sum_non_equal`](@ref).
"""
get_sum_body(x::SymbolicUtils.BasicSymbolic) = sum_body(x)
get_sum_body(x::Num) = get_sum_body(SymbolicUtils.unwrap(x))

"""
    indexed_sum(body, indices::Vector{Index}; non_equal = NonEqualPair[])

Construct a moment-layer indexed-sum node `∑_indices body` (see
[`is_indexed_sum`](@ref)), the same node [`average`](@ref) builds for a summed
[`QAdd`](@ref). `body` is an averaged scalar expression (a `Num` is unwrapped);
`non_equal` supplies pairwise index-inequality constraints (each `(α, β)` means
`α ≠ β`). The public constructor counterpart to the accessors
[`get_sum_body`](@ref), [`get_sum_indices`](@ref), and [`get_sum_non_equal`](@ref).

Use it to re-wrap a body in a sum node when a downstream rewrite (for example
cumulant expansion factorizing a summed moment into a product form) must preserve
its summation scope.

```jldoctest
julia> h = FockSpace(:site);

julia> @qnumbers a::Destroy(h);

julia> i = Index(h, :i, 3, h);

julia> x = average(Σ(IndexedOperator(a', i) * IndexedOperator(a, i), i));

julia> body = SecondQuantizedAlgebra.get_sum_body(x);

julia> isequal(SecondQuantizedAlgebra.indexed_sum(body, [i]), x)
true
```

See also [`get_sum_body`](@ref), [`is_indexed_sum`](@ref).
"""
function indexed_sum(body, indices::Vector{Index}; non_equal::Vector{NonEqualPair} = NonEqualPair[])
    return indexed_sum(SymbolicUtils.unwrap(body), indices, non_equal)
end

function symbolic_indices(@nospecialize(x))
    x isa RawExpression && return get_indices(x)
    x isa Op && return get_indices(x)
    x isa QAdd && return get_indices(x)
    x isa QExpr && return get_indices(x)
    return Index[]
end

function symbolic_acts_on(@nospecialize(x))
    x isa RawExpression && return acts_on(x)
    x isa Op && return acts_on(x)
    x isa QAdd && return acts_on(x)
    x isa QExpr && return acts_on(x)
    return Int[]
end

function get_indices(x::SymbolicUtils.BasicSymbolic)
    SymbolicUtils.isconst(x) && return symbolic_indices(x.val)
    SymbolicUtils.hasmetadata(x, AverageOperator) &&
        return symbolic_indices(SymbolicUtils.getmetadata(x, AverageOperator))
    is_indexed_sum(x) && return symbolic_indices(sum_body(x))
    SymbolicUtils.iscall(x) || return Index[]
    f = SymbolicUtils.operation(x)
    if f isa AvgFunc
        arg = SymbolicUtils.arguments(x)[1]
        inner = SymbolicUtils.isconst(arg) ? arg.val : arg
        return symbolic_indices(inner)
    end
    inds = Index[]
    for arg in SymbolicUtils.arguments(x), idx in symbolic_indices(arg)
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

"""
    acts_on(idx::Index) -> Vector{Int}

The subspace `idx` ranges over, as a one-element `Vector{Int}` holding its
`space_index`. Returns a vector (not a bare `Int`) so the return type matches the
operator method [`acts_on(::QSym)`](@ref). Pairs with [`operator_index`](@ref)
and [`set_acts_on`](@ref).
"""
acts_on(idx::Index) = Int[idx.space_index]

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
    SymbolicUtils.isconst(s) && return symbolic_acts_on(s.val)
    SymbolicUtils.hasmetadata(s, AverageOperator) &&
        return symbolic_acts_on(SymbolicUtils.getmetadata(s, AverageOperator))
    is_indexed_sum(s) && return symbolic_acts_on(sum_body(s))
    SymbolicUtils.iscall(s) || return Int[]
    f = SymbolicUtils.operation(s)
    f isa AvgFunc && return symbolic_acts_on(SymbolicUtils.arguments(s)[1])
    aon = Int[]
    for arg in SymbolicUtils.arguments(s)
        append!(aon, symbolic_acts_on(arg))
    end
    unique!(sort!(aon))
    return aon
end
