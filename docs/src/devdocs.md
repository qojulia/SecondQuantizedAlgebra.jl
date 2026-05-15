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
└── QAdd                          (sum of QTerm products — the only compound type)
```

`QTerm` is the per-entry storage key (operator product + non-equal constraints) used as the dict key inside `QAdd`. There is no abstract `QTerm` supertype — `QAdd` is a `QField` directly.

**Why `QAdd` is the only compound type.** Earlier versions of the package had both `QMul` (products) and `QAdd` (sums). This created a two-level expression tree where dispatch needed to handle `QSym`, `QMul`, and `QAdd` at every level, and the return type of `*` was unpredictable (`QSym`, `QMul`, or `QAdd` depending on simplification). The redesign collapses this into a single `QAdd` whose internal dictionary is keyed by the full term identity `(ops, ne)` and stores only the prefactor as the value. Every multiplication immediately produces a `QAdd`, giving a uniform return type. This type stability is critical for performance — the Julia compiler can infer return types through chains of arithmetic, avoiding dynamic dispatch and heap-allocated boxes at every intermediate step. The exact-key representation also keeps like-term collection honest: only terms with the same operator string *and* the same scoped constraints are merged.


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
const NonEqualPair = Tuple{Index, Index}

struct QTerm
    ops::Vector{QSym}
    ne::Vector{NonEqualPair}
end

const QTermDict = Dict{QTerm, CNum}

struct QAdd <: QField
    arguments::QTermDict
    indices::Vector{Index}
end
```

**What a dict entry represents.** A single entry `QTerm(ops, ne) => c` in `arguments` represents the term `c · ops[1] · ops[2] · …` valid for any index assignment satisfying the pairwise constraints in `ne` (each `(α, β) ∈ ne` means `α ≠ β`). The `indices` field on `QAdd` carries the outer summation scope: a `QAdd` with `indices = [i, j]` represents ``\sum_i \sum_j \sum_\text{terms} c \cdot \text{ops}`` where each individual term may further constrain `(i, j)` per its own `ne`. This per-term scoping is what lets a single `QAdd` represent expressions like ``\sum_i a_i a_i + \sum_{i \neq j} a_i a_j`` as two dict entries with different `ne` rather than two separate `QAdd` summations.

**Dictionary keys are full constrained terms.** Each `QTerm` is an ordered operator product plus the pairwise inequality constraints that scope that product. The ordering is canonical (determined by `_site_sort!` and the ordering convention), so structurally equal products always have the same `ops`. Like-term collection happens only when both `ops` and `ne` match exactly. This is the key invariant that makes constrained sums correct: `a_i` and `a_i` under `i ≠ j` are distinct stored terms, not one merged term with unioned metadata.

**The schema is visible at the type level.** `QTermDict` is a plain `Dict{QTerm, CNum}` alias — there is no wrapper struct. Iterating a `QAdd` yields `Pair{QTerm, CNum}`, and callers reach `term.ops` / `term.ne` on the key directly. This keeps the storage shape (`(ops, ne) => coeff`, not `ops => (coeff, ne)`) honest at every callsite.

**`_addto!` helper.** The core insertion function that handles like-term collection and zero elimination:

```julia
function _addto!(d::QTermDict, ops::Vector{QSym}, val::CNum, ne = _EMPTY_NE)
    key = QTerm(copy(ops), canonical_ne(ne))
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

**Summation metadata.** `indices` remains the sum-level metadata on `QAdd`: a `QAdd` with `indices = [i]` represents ``\sum_i (\text{terms})``. Pairwise inequality constraints like `(i, j)` meaning ``i \neq j`` live on the individual `QTerm`s, not on `QAdd` globally. This is why display and round-trip logic can represent mixed-scope sums truthfully.


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

**Three categories of rules:**

1. **Algebraic reductions** (ordering-independent, applied by `NormalOrder` at construction and by `simplify()`):
   - Transition composition: `|i⟩⟨j| · |k⟩⟨l| → δⱼₖ |i⟩⟨l|`
   - Pauli product: `σⱼ · σₖ = δⱼₖ I + iεⱼₖₗ σₗ`

2. **Ordering swaps** (only applied by `NormalOrder`, or explicitly by `normal_order()`):
   - Fock: `a · a† → a† · a + 1` (two terms: swapped + identity)
   - Spin: `[Sⱼ, Sₖ] = iεⱼₖₗ Sₗ` when axis out of order
   - PhaseSpace: `p · x → x · p - i`

3. **Completeness expansion** (only under `NormalOrder`, applied post-worklist by `_expand_gs_oterms`, or explicitly by `simplify(expr, h)` / `normal_order(expr, h)` under `LazyOrder`):
   - NLevel: `σᵍᵍ → 1 - Σ_{k≠g} σᵏᵏ` for the ground-state projector of any `NLevelSpace`. See "Simplification vs. normal ordering" below for the design rationale.

Each swap creates two terms (the swapped pair + the commutator term), so the worklist can grow. But like-term collection in `_addto!` keeps things manageable.

**`_same_site` guard.** All rules only fire when `a.space_index == b.space_index && a.index == b.index`. Operators on different sites trivially commute.

**`LazyOrder` bypass.** When the ordering is `LazyOrder`, `_apply_ordering` returns the input unchanged — no worklist, no scanning, and no algebraic reductions. All rules (both reductions and ordering swaps) are deferred until an explicit `simplify()` or `normal_order()` call:

```julia
_apply_ordering(arg_c::CNum, ops::Vector{QSym}, ::LazyOrder) =
    OrderedTerm[OrderedTerm(arg_c, ops)]
```

**Performance note.** When assembling large symbolic expressions where canonicalization only matters at the end (e.g. building a Hamiltonian via many additions and multiplications before computing a commutator or normal-ordering), wrap the construction in `with_ordering(LazyOrder()) do … end`. Operators accumulate without firing reductions or swaps on every product, and a single `normal_order(expr)` (or `simplify(expr, h)`) at the end produces canonical form in one pass. This avoids re-running the ordering worklist on every pairwise multiplication — for products with many same-site interactions, the savings compound.


## Simplification vs. normal ordering

The package exposes three normalization passes, each picking up a different subset of the rule categories above.

`simplify(expr)` (in `simplify.jl`) applies only the ordering-independent algebraic reductions (Transition composition, Pauli products), simplifies symbolic prefactors via `Symbolics.simplify`, and collects like terms. It never performs commutation swaps and never expands ground-state projectors. `normal_order(expr)` (in `normal_order.jl`) applies the full rule set — algebraic reductions plus commutation swaps plus completeness expansion — by calling `_apply_ordering(c, ops, NormalOrder())` regardless of the global ordering setting. `expand(expr)` expands symbolic prefactors only (`(a+b)² → a² + 2ab + b²`), leaving operator structure untouched.

**The `h`-overload and completeness.** `simplify(expr, h)` and `normal_order(expr, h)` are the explicit opt-in for ground-state completeness expansion under `LazyOrder` — they rewrite every `σᵍᵍ` (for any `NLevelSpace` subspace) as `1 - Σ_{k ≠ g} σᵏᵏ`. Under `NormalOrder` the `h` argument is a no-op cleanup pass, since `*` already eagerly canonicalizes every product. Each `Transition` carries its own `ground_state` and `n_levels` fields, so the algebra never consults `h` directly; the argument exists purely as the opt-in marker.

**Why `σᵍᵍ` is not in canonical form.** The canonical basis for `NLevelSpace` is `{σⁱʲ : (i,j) ≠ (g,g)} ∪ {1}` — the ground-state projector is deliberately excluded. The reason is the `QAdd = Dict{QTerm, CNum}` design invariant: physically equal expressions must have equal dict keys. If `σᵍᵍ` could live in canonical form, then `σᵍᵍ + σ²² + σ³³` (3-level, g=1) and `1` would be physically equal but compare unequal as dicts, breaking `isequal`, hash-based dedup, and `_addto!` merging. The eager rewrite preserves the invariant: every product passes through `_expand_gs_oterms` post-worklist (or through `_try_algebraic_reduction!` with `expand_gs = true` during same-site composition), so dict-key equality always implies physical equality. As a side effect, composition results like `σ¹² · σ²¹ → 1 - σ²² - σ³³` come out in the form a physicist would write directly, and downstream code (e.g. QuantumCumulants.jl meanfield) never has to dedupe an algebraically redundant `⟨σᵍᵍ⟩` moment against `1 - Σ ⟨σᵏᵏ⟩`.

User-constructed `σᵍᵍ` is the one exception. `Transition(h, :σ, g, g)` returns a plain `Transition` without expanding — canonicalization only fires when the operator enters a `*`. This keeps direct construction and inspection cheap, and is why the LazyOrder opt-in path is needed at all: under `LazyOrder` even composition-produced `σᵍᵍ` survives until the user invokes the `h`-overload.

**Cost.** Each `σᵍᵍ` reduction emits `n_levels` terms, so `k` reductions in one product cost `n_levels^k`. For typical workloads (2-level atoms, or 3-level with 1–2 atoms per product) this is small. The exponential only matters for high level counts combined with many `σᵍᵍ` factors in one product; if a future workload hits that, the design points either to deferring expansion to a single post-pass (which only helps when `σᵍᵍ`-bearing terms collide before expansion) or to promoting `σᵍᵍ` to first-class canonical with an explicit `apply_completeness` pass (cleaner separation, but breaks the dict-key invariant and forces a meanfield-side audit). Neither change is warranted until the cost actually bites.


## Diagonal splitting

When two `QAdd`s multiply (or a `QAdd` multiplies a `QSym`), the term-by-term product needs to handle the boundary case where two distinct-index operators land on the same site. `_accumulate_with_diag!` in `qadd_arithmetic.jl` does this in one pass; it is invoked from every `*(QAdd, ·)` overload after the per-term `vcat` of operand operators.

**The unified rule.** For each multiplied term, look at every pair of *distinct* free indices `(α, β)` on the same Hilbert subspace that isn't already constrained `α ≠ β`. The substitution `α → β` on the **unsorted** physical-order operator vcat gives the correct boundary value at `α = β`. The same substitution applied to the **post-sort, post-ordering** form may give a different result, because `_site_sort!` canonicalizes same-space ops alphabetically by index name and can put a same-site composition pair in the wrong order. When the two disagree, the boundary case is *not* implied by the off-diagonal term and must be emitted explicitly: add the unsorted-substituted contribution as a separate diagonal term, and re-key the original (off-diagonal) entry under the augmented constraint `ne ∪ {(α, β)}`.

**Two regimes for the disagreement check:**

- **Pairs `(i, j)` where `i` is a summation index of the outer `QAdd`.** Always emit. The sum over `i` runs across `i = j`; once the off-diagonal entry is constrained `i ≠ j`, the boundary value must surface as its own term regardless of pre-/post-sort agreement.
- **Free-index pairs `(α, β)`** (neither is a sum index). Only emit when `_oterms_equivalent(pre_sort, post_sort) == false`. If both substitutions land on the same multiset of ordered terms, the off-diagonal entry already carries the boundary contribution and no constraint is needed.

**Why pre-sort, not post-sort.** Once `_site_sort!` runs, same-space operators are sorted alphabetically by index name; substituting then would feed `_apply_ordering` a pair in the wrong physical order. Concretely, ``\sigma_j^{12} \cdot \sigma_i^{21}`` (with `i ≠ j`) sorts to ``\sigma_i^{21} \cdot \sigma_j^{12}`` and collapses at `i = j` to ``\sigma_j^{22}``, but the physically correct collapse from the original ordering is ``\sigma_j^{11}``. Substituting on the unsorted product preserves the physical order so same-site rules like ``\sigma^{αβ}\sigma^{βγ} = \sigma^{αγ}`` fire on the correct sequence.


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

**Summation metadata on averages.** When averaging an indexed `QAdd`, each averaged term carries its own `SumIndices` and `SumNonEqual` metadata via `SymbolicUtils.setmetadata`. This matches the internal term model: scoped constraints belong to individual terms, not to the whole sum. `undo_average` restores that exact metadata term-by-term, so symbolic sums survive the `average → manipulate → undo_average` round-trip even when the same operator string appears under different constraint scopes.

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
