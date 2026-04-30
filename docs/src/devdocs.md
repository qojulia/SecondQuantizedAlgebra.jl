# Developer Documentation

This page explains the internal architecture and design rationale of **SecondQuantizedAlgebra.jl**. It is intended for contributors and anyone who wants to understand *why* things are implemented the way they are.

## Type hierarchy

```
QField (abstract)
├── QSym (abstract) — atomic operators (leaves)
│   ├── Destroy / Create          (FockSpace)
│   ├── Transition                (NLevelSpace)
│   ├── Pauli                     (PauliSpace)
│   ├── Spin                      (SpinSpace)
│   └── Position / Momentum       (PhaseSpace)
└── QTerm (abstract) — compound expressions
    └── QAdd                      (the only concrete subtype)
```

**Why `QAdd` is the only compound type.** Earlier versions of the package had both `QMul` (products) and `QAdd` (sums). This created a two-level expression tree where dispatch needed to handle `QSym`, `QMul`, and `QAdd` at every level, and the return type of `*` was unpredictable (`QSym`, `QMul`, or `QAdd` depending on simplification). The redesign collapses this into a single `QAdd` that maps operator sequences (`Vector{QSym}`) to prefactors (`CNum`). Every multiplication immediately produces a `QAdd`, giving a uniform return type. This type stability is critical for performance — the Julia compiler can infer return types through chains of arithmetic, avoiding dynamic dispatch and heap-allocated boxes at every intermediate step. The dictionary representation also makes like-term collection automatic — adding two expressions with a common operator sequence just sums their prefactors at the same dictionary key, with no explicit simplification pass needed.


## Operator struct layout

Every `QSym` subtype shares a common field layout:

```julia
struct Destroy <: QSym
    name::Symbol        # display name
    space_index::Int    # which subspace in a ProductSpace (always 1 for simple spaces)
    index::Index        # symbolic summation index (NO_INDEX when absent)
end
```

**Why `space_index` instead of storing the Hilbert space.** Operators don't hold a reference to their Hilbert space. The space is only used at construction time for validation (checking bounds, matching types). At runtime, only the integer `space_index` matters — it determines which operators commute (different `space_index` → commute trivially) and where to embed in a composite basis during numeric conversion. This keeps operators lightweight and avoids type instability from heterogeneous space references.

**Why `index::Index` on every operator.** Symbolic summation indices (for expressions like ``\sum_i a_i^\dagger a_i``) live directly on the operator rather than in a wrapper type. This avoids a separate `IndexedOperator` struct in the type hierarchy and keeps dispatch simple. The sentinel `NO_INDEX` (with `space_index = 0`) indicates no index is present, checked via `has_index(idx) = idx.space_index != 0`. The `_same_site` check requires both `space_index` and `index` to match — operators with the same `space_index` but different `index` represent distinct sites in an indexed sum and don't interact via commutation rules.

**Transition-specific fields.** `Transition` adds `i::Int` and `j::Int` for the bra/ket levels. `Pauli` and `Spin` add `axis::Int` (1=x, 2=y, 3=z).

**`ladder` function.** Returns 0 for creation-like operators (`Create`, `Transition`, `Pauli`, `Spin`, `Position`) and 1 for annihilation-like operators (`Destroy`, `Momentum`). Used internally for canonical ordering.


## Hilbert spaces

Hilbert spaces are used at construction time only. The `ProductSpace{T}` uses a concrete `Tuple` type parameter:

```julia
struct ProductSpace{T <: Tuple{Vararg{HilbertSpace}}} <: HilbertSpace
    spaces::T
end
```

This gives full type stability — `ProductSpace{Tuple{FockSpace, NLevelSpace}}` is a concrete type. The `⊗` operator flattens nested `ProductSpace`s during construction.

**No `ClusterSpace`.** Earlier versions of SecondQuantizedAlgebra.jl had a `ClusterSpace{H,T}` wrapper, inherited from QuantumCumulants.jl, that encoded "this subspace has ``N`` identical copies, track correlations up to order ``k``." It was paired with a `copy_index::Int` field on every `QSym` so that `cluster_expand(op, N)` could materialize ``N`` distinct copies (`a_1, a_2, …, a_N`). QuantumCumulants.jl has since deprecated the cluster workflow (compare its [old cluster-based](https://qojulia.github.io/QuantumCumulants.jl/stable/examples/superradiant-laser/) and [new indexed](https://qojulia.github.io/QuantumCumulants.jl/stable/examples/superradiant_laser_indexed/) superradiant-laser tutorials), and we removed the corresponding scaffolding here because the symbolic [`Index`](@ref) machinery already covers the same use case at lower cost: an indexed operator `IndexedOperator(op, i)` with `i::Index` of range ``N`` represents the same family of copies, the algebra reasons over it via [`Σ`](@ref) and diagonal splitting, and ``N`` stays symbolic through equation derivation. The cumulant truncation order is properly a property of the moment expansion, so it now lives in the downstream solver call (`meanfield(…; order=k)` in QuantumCumulants.jl) rather than being baked into the Hilbert space. The integer-instantiation helpers (`insert_index`, `IndexedOperator(op, ::Int)`, `to_numeric(op, b, ranges)`) were removed alongside `copy_index` since they only made sense when concrete copies needed distinct numeric basis slots.


## The `CNum` prefactor type

```julia
const CNum = Complex{Num}
```

All prefactors are promoted to `CNum` — integers, floats, symbolic variables, complex numbers all get converted via `_to_cnum`. This ensures a single concrete prefactor type throughout the entire algebra, which is essential for type stability and therefore performance. The Julia compiler sees `Dict{Vector{QSym}, CNum}` as a fully concrete type, so dictionary operations, arithmetic, and iteration never trigger dynamic dispatch.

`Num` is the symbolic number type from Symbolics.jl, so `Complex{Num}` can represent both symbolic prefactors (`ω`, `g + im*κ`) and plain numeric ones (`3.0 + 0.0im`).

**Fast-path arithmetic.** Splitting into real and imaginary `Num` parts enables fast-path arithmetic. Most physics prefactors have zero imaginary part, so `_mul_cnum` and `_add_cnum` check `_iszero_num(c.im)` and skip half the arithmetic:

```julia
# Common case: both purely real → 1 multiply instead of 4
if ai_zero && bi_zero
    return Complex(ar * br, _NUM_ZERO)
end
```

When both operands are plain numbers (not symbolic), `_const_val` extracts them and does native Julia arithmetic, bypassing Symbolics entirely. Cached constants (`_CNUM_ZERO`, `_CNUM_ONE`, `_CNUM_IM`, etc.) avoid repeated allocations.


## QAdd internals

```julia
const QTermDict = Dict{Vector{QSym}, CNum}

struct QAdd <: QTerm
    arguments::QTermDict
    indices::Vector{Index}
    non_equal::Vector{Tuple{Index, Index}}
end
```

**Dictionary keys are `Vector{QSym}`.** Each key is an ordered sequence of operators representing a single product term. The ordering is canonical (determined by `_site_sort!` and the ordering convention), so structurally equal products always have the same key. This makes like-term collection automatic: adding `c₁ · a†a` and `c₂ · a†a` just adds `c₁ + c₂` to the same dictionary entry.

**`_addto!` helper.** The core insertion function that handles like-term collection and zero elimination:

```julia
function _addto!(d::QTermDict, key::Vector{QSym}, val::CNum)
    existing = get(d, key, nothing)
    if existing === nothing
        _iszero_cnum(val) || (d[key] = val)
    else
        new_val = _add_cnum(existing, val)
        _iszero_cnum(new_val) ? delete!(d, key) : (d[key] = new_val)
    end
end
```

Zero terms are never stored. This keeps the dictionary clean without needing explicit cleanup passes.

**Summation metadata.** `indices` and `non_equal` track symbolic sums. A `QAdd` with `indices = [i]` represents ``\sum_i (\text{terms})``. The `non_equal` field stores pairwise inequality constraints like `(i, j)` meaning ``i \neq j``. This metadata replaces the old `SingleSum`, `DoubleSum`, etc. types — there are no dedicated sum types anymore, just annotated `QAdd`s.


## Operator sorting

```julia
_sort_key(op::QSym) = (op.space_index, op.index.name)

_site_sort!(v::Vector{QSym}) = sort!(v; by = _sort_key, alg = Base.MergeSort)
```

**Why `MergeSort`?** Julia's default sort (`QuickSort`) is not stable — elements with equal keys can be reordered. Here, operators on the same site have identical `_sort_key` values, but their relative order encodes non-commutative multiplication (e.g. `a†·a ≠ a·a†`). A stable sort preserves this left-to-right order for equal keys while still grouping operators by site. `MergeSort` is the standard stable sort algorithm in Julia's `Base`.

**Why sort at all during multiplication?** Even with `LazyOrder`, operators from different sites are grouped together. This is correct because operators on different sites always commute, and grouping them makes the ordering worklist algorithm more efficient (it only needs to scan adjacent same-site pairs).


## Ordering: the worklist algorithm

The core of the algebra lives in `ordering.jl`. When two operators are multiplied, `_apply_ordering` processes the product through a worklist:

```
1. Start: worklist = [(prefactor, [op₁, op₂, ...])]
2. Pop a term, scan adjacent pairs left-to-right
3. If a rule fires (algebraic reduction or ordering swap):
   - Push the resulting term(s) back onto the worklist
   - Stop scanning this term (start over with next pop)
4. If no rule fires: term is fully ordered, push to "done"
5. Repeat until worklist is empty
```

**Two categories of rules:**

1. **Algebraic reductions** (ordering-independent, applied by `NormalOrder` at construction and by `simplify()`):
   - Transition composition: `|i⟩⟨j| · |k⟩⟨l| → δⱼₖ |i⟩⟨l|`
   - Pauli product: `σⱼ · σₖ = δⱼₖ I + iεⱼₖₗ σₗ`

2. **Ordering swaps** (only applied by `NormalOrder`, or explicitly by `normal_order()`):
   - Fock: `a · a† → a† · a + 1` (two terms: swapped + identity)
   - Spin: `[Sⱼ, Sₖ] = iεⱼₖₗ Sₗ` when axis out of order
   - PhaseSpace: `p · x → x · p - i`

Each swap creates two terms (the swapped pair + the commutator term), so the worklist can grow. But like-term collection in `_addto!` keeps things manageable.

**`_same_site` guard.** All rules only fire when `a.space_index == b.space_index && a.index == b.index`. Operators on different sites trivially commute.

**`LazyOrder` bypass.** When the ordering is `LazyOrder`, `_apply_ordering` returns the input unchanged — no worklist, no scanning, and no algebraic reductions. All rules (both reductions and ordering swaps) are deferred until an explicit `simplify()` or `normal_order()` call:

```julia
function _apply_ordering(arg_c::CNum, ops::Vector{QSym}, ::LazyOrder)
    return _OTerm[(arg_c, ops)]
end
```


## Simplification vs. normal ordering

**`simplify(expr)`** (in `simplify.jl`) applies only the ordering-independent algebraic reductions (Transition composition, Pauli products). It never applies commutation-based swaps. It also simplifies symbolic prefactors via `Symbolics.simplify` and collects like terms.

**`normal_order(expr)`** (in `normal_order.jl`) applies the full set of rules — algebraic reductions AND commutation swaps — by calling `_apply_ordering(c, ops, NormalOrder())` regardless of the global ordering setting.

**`expand(expr)`** expands symbolic prefactors (e.g. `(a+b)² → a² + 2ab + b²`) without touching operator structure.


## Diagonal splitting

When a symbolic sum ``\sum_i f(i)`` is multiplied by an operator with a free index `j` on the same space, the case `i = j` must be handled separately (the "diagonal" contribution). This is done by `_accumulate_with_diag!` in `qadd.jl`, which is invoked from each `*(QAdd, ·)` overload:

1. For each `(ops_a, ops_b)` term pair, vcat the operands (preserving physical multiplication order).
2. Site-sort + apply ordering on the unsorted vcat to produce the off-diagonal contribution.
3. For each summation index `i` the term depends on, find every distinct free index `j` (same space, not already constrained `i ≠ j`) and substitute `i → j` on the **unsorted** vcat — *before* `_site_sort!` would canonicalize same-space ops by index name. This preserves the physical operator order so that same-site composition rules (e.g. ``\sigma_{αβ}\sigma_{βγ} = \sigma_{αγ}``) fire on the correct sequence.
4. Site-sort and apply ordering on the substituted ops to get the diagonal contribution. Record the constraint `(i, j)` in `non_equal`.

Substitution must happen pre-sort: post-sort, the original physical order of same-space ops is lost (sorted alphabetically by index name), and substitution can collapse pairs in the wrong order — e.g. `σⱼ¹²·σᵢ²¹` (with `i ≠ j`) sorts to `σᵢ²¹·σⱼ¹²` and collapses at `i = j` to `σⱼ²²`, but the physically correct collapse is `σⱼ¹¹`.


## Index system

**`Index` struct:**
```julia
struct Index
    name::Symbol       # display name
    range::Num         # upper bound (Int or symbolic N)
    space_index::Int   # which space
    sym::Num           # Symbolics.jl variable for substitution in prefactors
end
```

**Why `sym` field?** When `change_index(expr, i, j)` operates on a symbolic prefactor like `g(i)`, it needs to substitute the Symbolics variable for `i` with the one for `j`. The `sym` field holds this variable. Without it, we'd need a separate mapping from index names to symbolic variables.

**`change_index(expr, from, to)`** performs symbolic substitution, replacing the index `from` with `to` throughout an expression tree (operator indices and the `sym` field of symbolic prefactors). Used for diagonal splitting and renaming sum indices.

**`NotIdentical` metadata.** `DoubleIndexedVariable(:g, i, j; identical=false)` creates `g(i,j)` with metadata `NotIdentical = true`. Two mechanisms enforce `g(i,i) = 0`: (1) at construction time, if `i == j` the function immediately returns `Num(0)`; (2) after a later substitution makes both arguments equal (e.g. `change_index` with `i → j`), `_check_not_identical` detects the equality and returns zero. Together these enforce `g(i,i) = 0` for off-diagonal coupling constants.


## Averaging

`average(expr)` converts operator expressions into symbolic scalars (SymbolicUtils `Term` nodes):

```julia
struct AvgFunc end                        # singleton callable, the "operation"
const sym_average = AvgFunc()

_average(op) = SymbolicUtils.Term{SymbolicUtils.SymReal}(sym_average, [op]; type = AvgSym)
```

**Why a custom `AvgFunc` instead of a `Sym`?** Using a custom struct lets us define `SymbolicUtils.show_call` for `⟨...⟩` display without type piracy on SymbolicUtils symbols.

**`AvgSym <: Number`** is the `symtype` marker. It ensures averaged expressions participate in Symbolics arithmetic (`+`, `*`, `simplify`) while remaining distinguishable from plain numbers.

### Why `average` returns `BasicSymbolic`, not `Num`

The Symbolics.jl ecosystem has two layers for symbolic expressions:
- `SymbolicUtils.BasicSymbolic{T}` — the raw symbolic tree
- `Symbolics.Num` — a `<: Real` wrapper that lets symbolic expressions participate in Julia's number promotion system

The standard Symbolics convention is for public APIs to wrap results in `Num` (e.g. `@variables x` produces `Num`, `Symbolics.derivative` returns `Num`). We deliberately break this convention for `average` because **`Num` wrapping is not type-stable across the `average(::QAdd)` code path**.

The problem: `average(::QAdd)` pulls c-number prefactors out of the average. When a prefactor involves a symbolic variable (e.g. `average(g * a)` where `g` is `@variables g::Complex`), the internal arithmetic `Num(0) + real(g) * BasicSymbolic_avg` goes through SymbolicUtils' promotion rules, which produce `BasicSymbolic`, not `Num`. Wrapping the other methods (`average(::QSym)`) in `Num` while `average(::QAdd)` returns `BasicSymbolic` creates an inconsistent return type — the same function returns different wrapper types depending on the input.

Rather than fighting the type system with explicit re-wrapping (which would be fragile and require `Num`-unwrapping dispatch methods on every downstream function), we return `BasicSymbolic` consistently:

- `average(::QSym)` → `BasicSymbolic{SymReal}` (the `Term` node directly)
- `average(::QAdd)` → `BasicSymbolic{SymReal}` (via `SymbolicUtils.unwrap` on the `Num` arithmetic result)
- `average(::Number)` → `Number` (scalars pass through unchanged)
- `average(::BasicSymbolic)` → `BasicSymbolic` (passthrough)

This means downstream functions (`is_average`, `acts_on`, `get_indices`, `numeric_average`, `undo_average`) can dispatch directly on `BasicSymbolic` without unwrapping. For robustness, `Num`-accepting dispatches are provided as defensive fallbacks — they unwrap via `SymbolicUtils.unwrap` and forward to the `BasicSymbolic` method. These catch cases where a user (or Symbolics arithmetic) wraps an average in `Num`.

`BasicSymbolic` values still support arithmetic (`+`, `*`, `^`) via SymbolicUtils, and SymbolicUtils keeps the result as `BasicSymbolic` even for expressions like `x - x = 0` (wrapped as `Const{SymReal}(0)`).

**Summation metadata on averages.** When averaging a `QAdd` with summation indices, the `SumIndices` and `SumNonEqual` metadata are attached to the resulting node via `SymbolicUtils.setmetadata`. The `undo_average` function restores this metadata when recovering operator expressions. This is how symbolic sums survive the `average → manipulate → undo_average` round-trip.

**`undo_average` always returns `QAdd`.** It accepts both `BasicSymbolic` and `Num` inputs (unwrapping `Num` first). Scalars become single-term `QAdd`s with an empty operator sequence, and lone `QSym`s become single-term `QAdd`s with unit prefactor. This uniform return type makes `undo_average` type-stable despite walking SymbolicUtils expression trees (where `operation(x)` returns `Any`).


## Numeric conversion

`to_numeric` maps symbolic operators to QuantumOpticsBase matrices:

- **Simple basis:** Direct dispatch — `Destroy → destroy(b)`, `Transition → transition(b, i, j)`, etc.
- **Composite basis:** Embeds the single-site matrix as a `LazyTensor`:
  ```julia
  op_num = to_numeric(op, b.bases[idx])
  LazyTensor(b, [idx], (op_num,))
  ```

**`_to_number`** extracts plain Julia numbers from `Num`/`CNum` wrappers for numeric evaluation. Falls back to the symbolic value if it can't be unwrapped (for symbolic prefactors that haven't been substituted yet).

**`_lazy_one`** creates the identity operator. For simple bases it returns `one(b)` (dense identity). For composite bases it returns a `LazyTensor` identity rather than materializing the full Kronecker-product identity matrix.


## Conjugation helpers (operators.jl)

Three internal conjugation functions handle `conj`/`adjoint` on symbolic expression trees that may contain averaged operator nodes (`⟨...⟩`):

- **`_conj(v)`**: Recursively conjugates symbolic trees. Leaves get `conj()`, QFields get `adjoint()`.
- **`_inconj(v)`**: Like `_conj` but aware of `AvgFunc` nodes — conjugates the operator *inside* the average via `adjoint`, preserving the `⟨...⟩` wrapper. Also handles nested `conj(avg(...))` by recursing into the argument.
- **`_adjoint(x)`**: Unified entry point dispatching to `adjoint` for operators, `_conj` for symbolic expressions.

These are needed because `conj(⟨a†a⟩)` should give `⟨a†a⟩` (Hermitian operator has real expectation value), not a naively conjugated symbolic tree. Standard `Base.conj` doesn't know how to recurse into SymbolicUtils expression trees. These helpers cannot be wired directly to `Base.conj`/`Base.adjoint` because that would be type piracy — defining methods on `Base` functions for `SymbolicUtils.BasicSymbolic`, a type we don't own. Instead, downstream packages (e.g. QuantumCumulants.jl) call these helpers explicitly when they need to conjugate averaged expressions.


## Printing and LaTeX

**Terminal printing** uses Unicode: `†` for dagger, subscript digits (`₀`-`₉`) for Transition levels, `σx`/`σy`/`σz` for Pauli axes. Summations render as `Σ(i=1:N)`.

**`sorted_arguments`** ensures deterministic output order. The sort key is `(length(ops), full_op_keys...)` where `_full_op_key(op) = (_sort_key(op)..., _type_order(op), op.name)`. This gives: shorter terms first, then by site, then by type (Destroy < Create < Transition < Pauli < Spin < Position < Momentum), then by name.

**LaTeX** uses Latexify.jl's `@latexrecipe` macro. `transition_superscript(::Bool)` toggles the global `transition_idx_script` `Ref` between `:^` and `:_`, controlling whether Transition level indices render as superscripts (`{name}^{{ij}}`) or subscripts (`{name}_{{ij}}`).
