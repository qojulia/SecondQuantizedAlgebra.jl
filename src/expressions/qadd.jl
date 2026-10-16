"""
    QAdd <: QField

Type of operator expressions that are polynomials in the operators, kept in canonical form.

Sums, products and commutators of operators return a `QAdd`, with every product ordered
and reduced and like terms collected, so two equal polynomials compare equal. A single
operator is an [`Op`](@ref) instead, and an expression that contains a formal function such
as `cos` is a [`QExpr`](@ref). `length` gives the number of terms, and
[`sorted_arguments`](@ref) splits an expression into them. The type is public but not
exported.

# Examples
```jldoctest
julia> using SecondQuantizedAlgebra: QAdd

julia> a = Destroy(FockSpace(:cavity), :a);

julia> x = a * a'
1 + a' * a

julia> x isa QAdd, length(x)
(true, 2)
```

See also [`QExpr`](@ref), [`normal_order`](@ref), [`sorted_arguments`](@ref).
"""
struct QAdd <: QField
    arguments::QTermDict
    indices::Vector{Index}
    function QAdd(arguments::QTermDict, indices::Vector{Index})
        return new(prune_dead_ne(arguments, indices), indices)
    end
end

function prune_dead_ne(args::QTermDict, indices::Vector{Index})
    needs_rebuild = false
    for (term, c) in args
        for p in term.ne
            if !pair_referenced(p, term.ops, c, indices)
                needs_rebuild = true
                break
            end
        end
        needs_rebuild && break
    end
    needs_rebuild || return args
    out = QTermDict()
    sizehint!(out, length(args))
    for (term, c) in args
        kept = NonEqualPair[]
        for p in term.ne
            pair_referenced(p, term.ops, c, indices) && push_ne_unique!(kept, p)
        end
        kept_ne = isempty(kept) ? EMPTY_NE : insertion_sort!(kept, isless)
        addto_key!(out, QTerm(copy(term.ops), kept_ne), c)
    end
    return out
end

function drop_unused_indices(d::QTermDict, indices::Vector{Index})
    isempty(indices) && return indices
    used = Index[idx for idx in indices if any_depends_on_index(d, idx)]
    length(used) == length(indices) && return indices
    return used
end

@inline function pair_referenced(
        p::NonEqualPair, ops::Vector{Op}, c::CNum, scope::Vector{Index},
    )
    α, β = p
    α in scope && return true
    β in scope && return true
    depends_on_index_term(c, ops, α) && return true
    depends_on_index_term(c, ops, β) && return true
    return false
end

"""
    constraint_pairs(expr) -> Vector{Tuple{Index, Index}}

Pairs of indices that `expr` requires to be distinct, collected over all its terms.

The function is not exported.

# Examples
```jldoctest
julia> using SecondQuantizedAlgebra: constraint_pairs

julia> h = NLevelSpace(:atom, 2); i = Index(h, :i, 10, h); j = Index(h, :j, 10, h);

julia> σ(k, l, m) = IndexedOperator(Transition(h, :σ, k, l), m);

julia> constraint_pairs(Σ(σ(2, 1, i) * σ(1, 2, j), i))
1-element Vector{Tuple{Index, Index}}:
 (i, j)
```

See also [`Σ`](@ref), [`assume_distinct_index`](@ref).
"""
function constraint_pairs(q::QAdd)
    return constraint_pairs(q.arguments)
end

function single_qadd(
        c::CNum, ops::Vector{Op}, ne::Vector{NonEqualPair} = EMPTY_NE,
        indices::Vector{Index} = EMPTY_INDICES,
    )
    iszero(c) && return QAdd(QTermDict(), EMPTY_INDICES)
    d = QTermDict()
    addto!(d, ops, c, ne)
    return QAdd(d, indices)
end

to_qadd(x::QAdd) = x
to_qadd(x::QSym) = single_qadd(CNUM_ONE, Op[x])
to_qadd(x::SymbolicUtils.BasicSymbolic) = single_qadd(to_cnum(x), Op[])

Base.length(a::QAdd) = length(a.arguments)
Base.iszero(a::QAdd) = isempty(a.arguments)

"""
    Base.one(q::QField) -> QAdd

Multiplicative identity as a unit `QAdd`.
"""
Base.one(::Type{<:QField}) = single_qadd(CNUM_ONE, Op[])
Base.one(q::QField) = one(typeof(q))

"""
    Base.isone(q::QField) -> Bool

True iff `q` is structurally the multiplicative identity.
"""
Base.isone(::QSym) = false
function Base.isone(q::QAdd)
    length(q.arguments) == 1 || return false
    (term, c) = first(q.arguments)
    isempty(term.ops) || return false
    is_native(c) && return isone(native_scalar(c))
    return isone(c)
end

function same_index_set(a::Vector{Index}, b::Vector{Index})
    length(a) == length(b) || return false
    for idx in a
        count(x -> isequal(x, idx), a) == count(x -> isequal(x, idx), b) || return false
    end
    return true
end

function Base.isequal(a::QAdd, b::QAdd)
    isequal(a.arguments, b.arguments) || return false
    return same_index_set(a.indices, b.indices)
end
Base.:(==)(a::QAdd, b::QAdd) = isequal(a, b)
Base.:(==)(a::QAdd, b::Op) = isequal(a, +b)
Base.:(==)(a::Op, b::QAdd) = isequal(+a, b)
function Base.hash(q::QAdd, h::UInt)
    ih = foldl((acc, idx) -> acc + hash(idx), q.indices; init = zero(UInt))
    return hash(:QAdd, hash(q.arguments, hash(ih, hash(length(q.indices), h))))
end

function Base.adjoint(q::QAdd)
    out = QTermDict()
    for (t, c) in q
        rev = Op[adjoint(o) for o in Iterators.reverse(t.ops)]
        canonicalize!(out, rev, conj(c), t.ne)
    end
    return QAdd(out, copy(q.indices))
end

Base.iterate(q::QAdd) = iterate(q.arguments)
Base.iterate(q::QAdd, state) = iterate(q.arguments, state)
Base.eltype(::Type{QAdd}) = Pair{QTerm, CNum}

ne_sort_key(ne::Vector{NonEqualPair}) = copy(ne)

function ne_less(a::Vector{NonEqualPair}, b::Vector{NonEqualPair})
    la, lb = length(a), length(b)
    @inbounds for i in 1:min(la, lb)
        a[i] != b[i] && return isless(a[i], b[i])
    end
    return la < lb
end

function sorted_args_less(x::Pair{QTerm, CNum}, y::Pair{QTerm, CNum})
    xo, yo = x.first.ops, y.first.ops
    lx, ly = length(xo), length(yo)
    lx != ly && return lx < ly
    @inbounds for i in 1:lx
        kx, ky = full_op_key(xo[i]), full_op_key(yo[i])
        kx != ky && return isless(kx, ky)
    end
    return ne_less(x.first.ne, y.first.ne)
end

"""
    sorted_arguments(expr)

Terms of `expr`, each as an expression of its own, in a fixed order.

Each term can then be read with [`get_prefactor`](@ref) and
[`get_operators`](@ref). The function is not exported.

# Examples

```jldoctest
julia> using SecondQuantizedAlgebra: sorted_arguments

julia> a = Destroy(FockSpace(:cavity), :a);

julia> @variables ω g;

julia> foreach(println, sorted_arguments(ω * a' * a + g * (a + a')))
g * a
g * a'
ω * a' * a
```

See also [`get_prefactor`](@ref), [`get_operators`](@ref).
"""
function sorted_arguments(q::QAdd)
    isempty(q.arguments) && return QAdd[]
    pairs = collect(q.arguments)
    insertion_sort!(pairs, sorted_args_less)
    return QAdd[single_qadd(c, term.ops, term.ne, q.indices) for (term, c) in pairs]
end

full_op_key(op::QSym) = (
    sort_key(op)..., Int(op.kind), name_rank(op.name_id),
    Int(op.l1), Int(op.l2), Int(op.g), Int(op.nlev),
)

"""
    term_order_key(t::QTerm) -> Tuple

Total ordering key for an operator product, built from `order_key`: length, then
operator-by-operator, then the canonical-sorted non-equal constraints.
"""
term_order_key(t::QTerm) = (length(t.ops), map(order_key, t.ops), ne_sort_key(t.ne))

sort_key(op::QSym) = (op.space_index, name_rank(op.index.name_id))

"""
    qadd_order_key(q::QAdd) -> Vector

Total, reproducible ordering key for a sum: its term/coefficient pairs in sorted order, so
two `QAdd`s compare with `<` on their keys and tie exactly when they are `isequal`. The
coefficient contributes its printed form (a reproducible tiebreak, not a numeric order).
"""
function qadd_order_key(q::QAdd)
    pairs = [(term_order_key(t), coeff_key(c)) for (t, c) in q.arguments]
    insertion_sort!(pairs, isless)
    return pairs
end

Base.isless(a::QAdd, b::QAdd) = isless(qadd_order_key(a), qadd_order_key(b))
Base.isless(a::Type{<:QField}, b::Type{<:QField}) = isless(nameof(a), nameof(b))

function coeff_key(c::CNum)
    re, im = realimag(c)
    return (string(re), string(im))
end

"""
    Base.getindex(q::QAdd, key::AbstractVector{<:QSym}) -> Complex{Num}

Look up the prefactor for a given operator sequence. Returns zero if absent.
Throws if more than one constrained term shares the same operator sequence.
"""
function Base.getindex(q::QAdd, key::AbstractVector{<:QSym})
    found = CNUM_ZERO
    matches = 0
    for (term, c) in q.arguments
        term.ops == key || continue
        matches += 1
        matches > 1 && throw(
            ArgumentError(
                "operator sequence has multiple constrained terms; iterate `q.arguments` " *
                    "or `sorted_arguments(q)` to inspect them explicitly"
            )
        )
        found = c
    end
    return to_num(found)
end

"""
    Base.haskey(q::QAdd, key::AbstractVector{<:QSym}) -> Bool

Return `true` if some stored term has operator sequence `key`, ignoring
constraint scope. Pair with [`Base.getindex`](@ref) for the unique-prefactor
lookup; iterate `q.arguments` for full scope-aware access.
"""
function Base.haskey(q::QAdd, key::AbstractVector{<:QSym})
    for term in keys(q.arguments)
        term.ops == key && return true
    end
    return false
end

"""
    get_prefactor(expr)

Coefficient of a single-term expression `expr`, as a `Complex{Num}`.

An expression with more than one term throws an `ArgumentError`.

# Examples

```jldoctest
julia> a = Destroy(FockSpace(:cavity), :a);

julia> @variables g;

julia> get_prefactor(2g * a' * a)
2g
```

See also [`get_operators`](@ref), [`get_variables`](@ref).
"""
function get_prefactor(s::QAdd)
    length(s.arguments) == 1 || throw(
        ArgumentError("get_prefactor requires a single-term expression, got $(length(s.arguments)) terms")
    )
    return to_num(first(values(s.arguments)))
end

"""
    get_operators(expr)

Distinct operators occurring in `expr`, in canonical order.

An operator and its adjoint count as two entries, and repeated factors appear
once.

# Examples

```jldoctest
julia> a = Destroy(FockSpace(:cavity), :a);

julia> @variables g;

julia> get_operators(g * a' * a + a)
2-element Vector{Op}:
 a
 a'
```

See also [`get_variables`](@ref), [`unique_up_to_adjoint`](@ref).
"""
function get_operators(s::QAdd)
    ops = Op[]
    for term in keys(s.arguments)
        append!(ops, term.ops)
    end
    unique!(ops)
    sort!(ops; by = order_key)
    return ops
end

"""
    get_sum_indices(q::QAdd) -> Vector{Index}
    get_sum_indices(q::QSym) -> Vector{Index}

The summation indices bound by `q`, in the order of its sum scope. A `QSym` and an unsummed
`QAdd` bind none. The result is a copy, so mutating it does not change `q`.

```jldoctest
julia> h = FockSpace(:site);

julia> @qnumbers a::Destroy(h);

julia> i = Index(h, :i, 3, h);

julia> SecondQuantizedAlgebra.get_sum_indices(Σ(IndexedOperator(a, i), i)) == [i]
true
```

See also [`get_indices`](@ref), [`constraint_pairs`](@ref).
"""
get_sum_indices(q::QAdd) = copy(q.indices)
get_sum_indices(::QSym) = Index[]

"""
    get_variables(expr)

Symbolic parameters occurring in the coefficients of `expr`.

The result can be passed to Symbolics functions such as `build_function`. An
optional second argument restricts it to a given list of variables.

# Examples

```jldoctest
julia> a = Destroy(FockSpace(:cavity), :a);

julia> @variables ω g;

julia> H = ω * a' * a + g * (a + a');

julia> length(get_variables(H))
2
```

See also [`get_prefactor`](@ref), [`get_operators`](@ref), [`substitute`](@ref).
"""
function Symbolics.get_variables(s::QAdd; kwargs...)
    vars = SymbolicUtils.BasicSymbolic[]
    for c in values(s.arguments)
        append!(vars, Symbolics.get_variables(real(c); kwargs...))
        append!(vars, Symbolics.get_variables(imag(c); kwargs...))
    end
    unique!(vars)
    sort!(vars; by = string)
    return vars
end

function Symbolics.get_variables(
        s::QAdd,
        varlist;
        is_atomic = SymbolicUtils.default_is_atomic,
        kwargs...,
    )
    vars = SymbolicUtils.BasicSymbolic[]
    for c in values(s.arguments)
        append!(
            vars,
            Symbolics.get_variables(
                real(c), varlist; is_atomic = is_atomic, kwargs...
            ),
        )
        append!(
            vars,
            Symbolics.get_variables(
                imag(c), varlist; is_atomic = is_atomic, kwargs...
            ),
        )
    end
    unique!(vars)
    sort!(vars; by = string)
    return vars
end

function Symbolics.get_variables!(buffer, s::QAdd; kwargs...)
    for c in values(s.arguments)
        Symbolics.get_variables!(buffer, real(c); kwargs...)
        Symbolics.get_variables!(buffer, imag(c); kwargs...)
    end
    return buffer
end

function Symbolics.get_variables!(
        buffer,
        s::QAdd,
        varlist;
        is_atomic = SymbolicUtils.default_is_atomic,
        kwargs...,
    )
    for c in values(s.arguments)
        Symbolics.get_variables!(
            buffer, real(c), varlist; is_atomic = is_atomic, kwargs...
        )
        Symbolics.get_variables!(
            buffer, imag(c), varlist; is_atomic = is_atomic, kwargs...
        )
    end
    return buffer
end

"""
    coefficient(expr, term) -> Coeff

The coefficient of the operator product `term` in the polynomial `expr`, or zero when `expr`
has no such term. `term` is an operator, a product such as `a' * a`, or `1` for the scalar
part. It must canonicalize to a single product with unit coefficient, so `2a`, `a + a'` and
`a * a'` (which is `a'a + 1`) are refused. Index inequalities are part of a term: a product
that `expr` stores under other constraints, or under several, throws instead of reading as
zero; build `term` with [`assume_distinct_index`](@ref) to match a constrained product.
Summed expressions are refused.

Unlike [`coefficient_terms`](@ref), which splits a scalar into its parameter monomials,
`coefficient` reads the scalar that multiplies an operator product.

```jldoctest
julia> using SecondQuantizedAlgebra: coefficient

julia> @qnumbers a::Destroy(FockSpace(:c));

julia> H = 3 * a' * a + 2 * a' + 2 * a + 5;

julia> coefficient(H, a' * a), coefficient(H, a'), coefficient(H, 1), coefficient(H, a' * a')
(3, 2, 5, 0)
```
"""
function coefficient(expr::Union{QField, Coefficient}, term::Union{QField, Coefficient})
    q = coefficient_operand(expr, "expression")
    isempty(q.indices) || summed_coefficient_error()
    key = monomial_key(coefficient_operand(term, "term"))
    found = CNUM_ZERO
    matched = false
    for (t, c) in q.arguments
        isequal(t.ops, key.ops) || continue
        (matched || !isequal(t.ne, key.ne)) && ambiguous_coefficient_error()
        found = c
        matched = true
    end
    return found
end

coefficient_operand(x::QAdd, ::String) = x
coefficient_operand(x::QSym, ::String) = to_qadd(x)
coefficient_operand(x::Coefficient, ::String) = single_qadd(to_cnum(x), Op[])
coefficient_operand(::QField, what::String) = throw(
    ArgumentError(
        "`coefficient` needs a polynomial $what; lower a formal expression with `taylor` first"
    )
)

function monomial_key(q::QAdd)
    isempty(q.indices) || summed_coefficient_error()
    length(q.arguments) == 1 || throw(
        ArgumentError(
            "`coefficient` needs a single operator product with unit coefficient, but the " *
                "term has $(length(q.arguments)) canonical terms"
        )
    )
    t, c = first(q.arguments)
    isone(c) || throw(
        ArgumentError(
            "`coefficient` needs a single operator product with unit coefficient, but the " *
                "term has coefficient $(to_num(c))"
        )
    )
    return t
end

@noinline summed_coefficient_error() = throw(
    ArgumentError("`coefficient` does not read summed expressions; iterate the terms of the sum instead")
)

@noinline ambiguous_coefficient_error() = throw(
    ArgumentError(
        "the operator product is stored under other or several index constraints; build the " *
            "term with `assume_distinct_index` to match, or iterate the expression's terms"
    )
)
