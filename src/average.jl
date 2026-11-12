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

SymbolicUtils.:<ₑ(a::Op, b::Op) = canonical_isless(a, b)
SymbolicUtils.:<ₑ(a::QAdd, b::QAdd) = canonical_isless(a, b)

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
`==`/`isequal`/`hash` compare it by value, so two sums cancel only when their
scopes are equal.
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

Whether `x` is a sum of averages over an index, as made by [`average`](@ref) of a [`Σ`](@ref).

The indices and constraints of the sum are read with [`get_sum_indices`](@ref) and
[`get_sum_non_equal`](@ref), and its summand with [`get_sum_body`](@ref).

# Examples
```jldoctest
julia> using SecondQuantizedAlgebra: is_indexed_sum

julia> h = NLevelSpace(:atom, 2); i = Index(h, :i, 10, h);

julia> x = average(Σ(IndexedOperator(Transition(h, :σ, 2, 2), i), i))
Σ(i=1:10) ⟨σ_i₂₂⟩

julia> is_indexed_sum(x)
true
```

See also [`indexed_sum`](@ref), [`average`](@ref).
"""
is_indexed_sum(::Any) = false
is_indexed_sum(x::SymbolicUtils.BasicSymbolic) =
    SymbolicUtils.iscall(x) && SymbolicUtils.operation(x) isa SumFunc
is_indexed_sum(x::Num) = is_indexed_sum(SymbolicUtils.unwrap(x))

"""
    is_average(x) -> Bool

Whether `x` is the average of an operator expression, as made by [`average`](@ref).

# Examples
```jldoctest
julia> a = Destroy(FockSpace(:cavity), :a);

julia> is_average(average(a)), is_average(a)
(true, false)
```

See also [`average`](@ref), [`is_indexed_sum`](@ref).
"""
is_average(::Any) = false
function is_average(x::SymbolicUtils.BasicSymbolic)
    (SymbolicUtils.iscall(x) && SymbolicUtils.operation(x) === sym_average) && return true
    return SymbolicUtils.hasmetadata(x, AverageOperator)
end
is_average(x::Num) = is_average(SymbolicUtils.unwrap(x))

avg_symtype(op::QField) = isequal(adjoint(op), op) ? Real : Number

sum_symtype(body::RawExpression) = SymbolicUtils.symtype(body) <: Real ? Real : Number
sum_symtype(::Real) = Real
sum_symtype(@nospecialize(body)) = Number

avg_symtype_from_arg(arg::RawExpression) =
    SymbolicUtils.isconst(arg) ? avg_symtype_from_arg(arg.val) : Number
avg_symtype_from_arg(op::Union{Op, QAdd}) = avg_symtype(op)
avg_symtype_from_arg(@nospecialize(arg)) = Number

average_leaf(op::QField) = SymbolicUtils.Term{SymbolicUtils.SymReal}(sym_average, QField[op]; type = avg_symtype(op))

"""
    averaged_operator(x)

The operator of an average: the argument of `⟨X⟩`, or the operator carried by its lifted
form `⟨X⟩(t)`.
"""
function averaged_operator(x::RawExpression)
    SymbolicUtils.hasmetadata(x, AverageOperator) && return SymbolicUtils.getmetadata(x, AverageOperator)
    arg = SymbolicUtils.arguments(x)[1]
    return SymbolicUtils.isconst(arg) ? arg.val : arg
end

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
    make_time_dependent(expr, t)

Replace every average in `expr` by an unknown function of time `t`.

Each average ``\\langle X \\rangle`` becomes ``\\langle X \\rangle(t)``, which can be
differentiated with respect to `t` and serves as an unknown of a differential equation, for
instance in ModelingToolkit.jl. The unknown is real when ``X`` is Hermitian and complex
otherwise. [`undo_average`](@ref) recovers the operators.

# Examples
```jldoctest
julia> a = Destroy(FockSpace(:cavity), :a);

julia> @variables t;

julia> make_time_dependent(average(a' * a) + average(a), t)
⟨a⟩(t) + ⟨a' * a⟩(t)
```

See also [`average`](@ref), [`undo_average`](@ref).
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
    fn = SymbolicUtils.Sym{SymbolicUtils.SymReal}(
        nm; type = SymbolicUtils.FnType{Tuple, avg_symtype(op), Nothing}, shape = UnitRange{Int}[],
    )
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
    average(expr)

Expectation value of an operator expression, as a symbolic scalar.

The average is linear: it distributes over sums and takes coefficients out, so the result
is a sum of averages of operator products, printed as `⟨…⟩`. It takes part in ordinary
symbolic arithmetic with Symbolics.jl. The average of a Hermitian operator is real. A sum
over an index becomes a sum of averages over the same index.

# Examples
```jldoctest
julia> a = Destroy(FockSpace(:cavity), :a);

julia> @variables ω;

julia> average(ω * a' * a + a)
⟨a⟩ + ⟨a' * a⟩*ω
```

See also [`undo_average`](@ref), [`numeric_average`](@ref), [`make_time_dependent`](@ref).
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

average(::AnyMap) = throw(
    ArgumentError(
        "a superoperator has no average; average the Heisenberg drift `transpose(L)(A)` instead",
    ),
)

function average(op::QAdd)
    BS = SymbolicUtils.BasicSymbolic{SymbolicUtils.SymReal}
    GroupKey = Tuple{Vector{Index}, Vector{NonEqualPair}}
    flat = BS[]
    groups = Pair{GroupKey, Vector{BS}}[]
    slots = Dict{GroupKey, Int}()
    for (term, c) in op.arguments
        body = average_term(c, term)
        used = bound_indices(term, c, op.indices)
        isempty(used) && (push!(flat, body); continue)
        slot = get!(slots, (used, term.ne)) do
            push!(groups, (used, copy_ne(term.ne)) => BS[])
            length(groups)
        end
        push!(last(groups[slot]), body)
    end
    for ((indices, ne), bodies) in groups
        push!(flat, indexed_sum(add_worker(SymbolicUtils.SymReal, bodies), indices, ne))
    end
    isempty(flat) && return SymbolicUtils.unwrap(Num(0))
    return add_worker(SymbolicUtils.SymReal, flat)
end

function average_term(c::CNum, term::QTerm)
    avg = isempty(term.ops) ? 1 : average_term_ops(term)
    if c.tail isa RawSymbolicCoeff && !is_call_of(c.tail.expr, complex)
        return SymbolicUtils.unwrap(Num(c.tail.expr) * avg)
    end
    tail = c.tail
    tail isa Poly && has_complex_atom(tail) && return SymbolicUtils.unwrap(Num(poly_to_raw(tail)) * avg)
    r, i = realimag(c)
    contrib = Num(0)
    iszero(r) || (contrib += r * avg)
    iszero(i) || (contrib += i * Symbolics.IM * avg)
    return SymbolicUtils.unwrap(contrib)
end

function average_term_ops(term::QTerm)
    inner = length(term.ops) == 1 ? only(term.ops) : single_qadd(CNUM_ONE, term.ops, term.ne)
    return average_leaf(inner)
end

function rebuild_indexed_sum(inner::QAdd, indices::Vector{Index}, ne::Vector{NonEqualPair})
    isempty(ne) && return QAdd(inner.arguments, indices)
    new_args = QTermDict()
    for (term, c) in inner.arguments
        addto!(new_args, term.ops, c, merge_ne(term.ne, ne))
    end
    return QAdd(new_args, indices)
end

contains_average(x) = has_subterm(y -> is_average(y) || is_indexed_sum(y), SymbolicUtils.unwrap(x))

"""
    undo_average(expr)
    undo_average(lhs ~ rhs)

Operator expression whose average is `expr`, the inverse of [`average`](@ref).

Sums over indices are restored, and so are averages turned into functions of time by
[`make_time_dependent`](@ref). An equation gives the pair `lhs => rhs` of operator
expressions.

# Examples
```jldoctest
julia> a = Destroy(FockSpace(:cavity), :a);

julia> @variables ω;

julia> undo_average(average(ω * a' * a + a))
a + ω * a' * a
```

See also [`average`](@ref).
"""
function undo_average(x::SymbolicUtils.BasicSymbolic)
    is_average(x) && return to_qadd(averaged_operator(x))
    if is_indexed_sum(x)
        return rebuild_indexed_sum(undo_average(sum_body(x)), sum_indices(x), sum_ne(x))
    end
    SymbolicUtils.iscall(x) || return to_qadd(x)
    f = SymbolicUtils.operation(x)
    args = SymbolicUtils.arguments(x)
    if f === (/) && length(args) == 2 && !contains_average(args[2])
        return undo_average(args[1]) * inv(to_cnum(args[2]))
    end
    f === (+) && return sum(QAdd[undo_average(a) for a in args])
    f === (*) && return prod(QAdd[undo_average(a) for a in args])
    return to_qadd(x)
end

undo_average(x::Number) = single_qadd(to_cnum(x), Op[])
undo_average(x::Num) = undo_average(SymbolicUtils.unwrap(x))
undo_average(x::QSym) = to_qadd(x)
undo_average(x::QAdd) = x
undo_average(eq::Symbolics.Equation) = undo_average(eq.lhs) => undo_average(eq.rhs)

"""
    has_sum_metadata(x) -> Bool

Same as [`is_indexed_sum`](@ref), kept under its earlier name. The function is not exported.
"""
has_sum_metadata(::Any) = false
has_sum_metadata(x::SymbolicUtils.BasicSymbolic) = is_indexed_sum(x)
has_sum_metadata(x::Num) = has_sum_metadata(SymbolicUtils.unwrap(x))

"""
    get_sum_indices(x) -> Vector{Index}

Indices summed over in the sum of averages `x`.

The function is not exported.

# Examples
```jldoctest
julia> using SecondQuantizedAlgebra: get_sum_indices

julia> h = NLevelSpace(:atom, 2); i = Index(h, :i, 10, h);

julia> get_sum_indices(average(Σ(IndexedOperator(Transition(h, :σ, 2, 2), i), i)))
1-element Vector{Index}:
 i
```

See also [`is_indexed_sum`](@ref), [`get_sum_non_equal`](@ref), [`get_sum_body`](@ref).
"""
get_sum_indices(x::SymbolicUtils.BasicSymbolic) = sum_indices(x)
get_sum_indices(x::Num) = get_sum_indices(SymbolicUtils.unwrap(x))

"""
    get_sum_non_equal(x) -> Vector{Tuple{Index, Index}}

Pairs of indices required to differ in the sum of averages `x`.

The function is not exported.

# Examples
```jldoctest
julia> using SecondQuantizedAlgebra: get_sum_non_equal

julia> h = NLevelSpace(:atom, 2); i = Index(h, :i, 10, h); j = Index(h, :j, 10, h);

julia> get_sum_non_equal(average(Σ(IndexedOperator(Transition(h, :σ, 2, 2), i), i, [j])))
1-element Vector{Tuple{Index, Index}}:
 (i, j)
```

See also [`is_indexed_sum`](@ref), [`get_sum_indices`](@ref).
"""
get_sum_non_equal(x::SymbolicUtils.BasicSymbolic) = sum_ne(x)
get_sum_non_equal(x::Num) = get_sum_non_equal(SymbolicUtils.unwrap(x))

"""
    get_sum_body(x)

Summand of the sum of averages `x`, with its index still free.

The function is not exported.

# Examples
```jldoctest
julia> using SecondQuantizedAlgebra: get_sum_body

julia> h = NLevelSpace(:atom, 2); i = Index(h, :i, 10, h);

julia> get_sum_body(average(Σ(IndexedOperator(Transition(h, :σ, 2, 2), i), i)))
⟨σ_i₂₂⟩
```

See also [`indexed_sum`](@ref), [`get_sum_indices`](@ref).
"""
get_sum_body(x::SymbolicUtils.BasicSymbolic) = sum_body(x)
get_sum_body(x::Num) = get_sum_body(SymbolicUtils.unwrap(x))

"""
    indexed_sum(body, indices::Vector{Index}; non_equal = Tuple{Index, Index}[])

Sum of the averaged expression `body` over `indices`, the inverse of [`get_sum_body`](@ref).

`non_equal` lists pairs of indices that must differ. This rebuilds a sum after its summand
has been rewritten, for instance when a cumulant expansion factorizes the averages in it.
The function is not exported.

# Examples
```jldoctest
julia> using SecondQuantizedAlgebra: indexed_sum, get_sum_body

julia> h = NLevelSpace(:atom, 2); i = Index(h, :i, 10, h);

julia> x = average(Σ(IndexedOperator(Transition(h, :σ, 2, 2), i), i));

julia> isequal(indexed_sum(get_sum_body(x), [i]), x)
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

"""
    foreach_operand(f, x::RawExpression)

Call `f` on every operand of an expression of averages: each operator or expression held by
a constant, and the operator of each average.
"""
function foreach_operand(f::F, x::RawExpression) where {F}
    if SymbolicUtils.isconst(x)
        visit_operand(f, x.val)
    elseif is_average(x)
        visit_operand(f, averaged_operator(x))
    elseif is_indexed_sum(x)
        foreach_operand(f, sum_body(x))
    elseif SymbolicUtils.iscall(x)
        for arg in SymbolicUtils.arguments(x)
            foreach_operand(f, arg)
        end
    end
    return nothing
end

function visit_operand(f::F, @nospecialize(v)) where {F}
    v isa Op && return f(v)
    v isa QAdd && return f(v)
    v isa QExpr && return f(v)
    v isa RawExpression && return f(v)
    return nothing
end

function get_indices(x::SymbolicUtils.BasicSymbolic)
    inds = Index[]
    foreach_operand(v -> merge_indices!(inds, symbolic_indices(v)), x)
    return inds
end

function acts_on(s::SymbolicUtils.BasicSymbolic)
    aon = Int[]
    foreach_operand(v -> append!(aon, symbolic_acts_on(v)), s)
    return unique!(sort!(aon))
end

"""
    qadjoint(x)
    qconj(x)

Hermitian conjugate of an operator expression or a coefficient.

For an operator expression this is `x'`. For a coefficient the conjugate is taken through
the expression tree, so that it reaches the individual parameters instead of wrapping the
whole expression in `conj`. An average becomes the average of the adjoint, as with
[`inner_adjoint`](@ref).

# Examples
```jldoctest
julia> a = Destroy(FockSpace(:cavity), :a);

julia> @variables g::Number;

julia> qadjoint(2im * a), qadjoint(g)
(-2im * a', conj(g))
```

See also [`dagger`](@ref), [`inner_adjoint`](@ref).
"""
function qadjoint(v::SymbolicUtils.BasicSymbolic)
    v === Symbolics.IM && return -Symbolics.IM
    is_average(v) && return adjoint_average(v)
    if is_indexed_sum(v)
        return indexed_sum(qadjoint(sum_body(v)), sum_indices(v), sum_ne(v))
    end
    SymbolicUtils.iscall(v) || return conj(v)
    f = SymbolicUtils.operation(v)
    args = SymbolicUtils.arguments(v)
    f === conj && return args[1]
    f === expim && return expim_symbolic(-only(args))
    return TermInterface.maketerm(typeof(v), f, map(qadjoint, args), TermInterface.metadata(v))
end
qadjoint(x::Number) = adjoint(x)
qadjoint(x::Num) = qadjoint(SymbolicUtils.unwrap(x))
qadjoint(x::QField) = adjoint(x)

const qconj = qadjoint

"""
    dagger(x::QField)

Hermitian conjugate of an operator expression, the same as `x'`.

Coefficients are conjugated along with the operators. The name matches `dagger`
from QuantumOptics.jl, so code can use it for symbolic and numeric operators
alike.

# Examples

```jldoctest
julia> h = FockSpace(:cavity) ⊗ NLevelSpace(:atom, 2);

julia> @qnumbers a::Destroy(h, 1) σ::Transition(h, 1, 2, 2);

julia> dagger(2im * a' * σ)
-2im * a * σ₂₁
```

See also [`qadjoint`](@ref).
"""
QuantumInterface.dagger(x::QField) = qadjoint(x)

"""
    inner_adjoint(x)

Complex conjugate of an expression of averages, written again as averages.

The conjugate of ``\\langle X \\rangle`` is written as ``\\langle X^\\dagger \\rangle``
rather than as `conj(⟨X⟩)`, so that equations for averages stay expressed in averages of
operators. Coefficients are conjugated as by [`qadjoint`](@ref).

# Examples
```jldoctest
julia> a = Destroy(FockSpace(:cavity), :a);

julia> @variables g::Number;

julia> inner_adjoint(g * average(a))
⟨a'⟩*conj(g)
```

See also [`qadjoint`](@ref), [`average`](@ref).
"""
inner_adjoint(x::Union{RawExpression, Num, Number, QField}) = qadjoint(x)

function adjoint_average(v::RawExpression)
    op = averaged_operator(v)
    SymbolicUtils.hasmetadata(v, AverageOperator) &&
        return lift_average(average(adjoint(op::QAdd)), SymbolicUtils.arguments(v)[1])
    return adjoint_average_leaf(op)
end

function adjoint_average_leaf(@nospecialize(inner))
    inner isa Op && return average_leaf(adjoint(inner))
    inner isa QAdd && return average_leaf(adjoint(inner))
    throw(ArgumentError("cannot take the adjoint of the averaged argument `$inner`"))
end
