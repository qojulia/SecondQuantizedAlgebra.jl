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

**Transition-specific fields.** `Transition` adds `i::Int` and `j::Int` for the bra/ket levels along with `ground_state::Int` and `n_levels::Int`; `Pauli` and `Spin` add `axis::Int` (1=x, 2=y, 3=z).

**The five operator hooks.** The whole algebra talks to operator types through five methods that every concrete `QSym` subtype implements. They are the entire interface; everything else in the package builds on top of them.

| Hook | Returns | Meaning |
|---|---|---|
| `_site_compare(a, b, ne)` | `SiteCmp` | three-way site comparison driving the sort |
| `_can_commute(a, b)` | `Bool` | true iff swapping needs no commutator residual (called only on provably-same-site pairs) |
| `_commute_pair(a, b)` | `(swap_b, swap_a, residual_coeff, residual_ops)` | swap and residual for same-site non-commuting pairs |
| `_reduce_pair(a, b)` | `Nothing`, `CNum`, or `(QSym, CNum)` | local algebraic identity (Transition composition, Pauli product, …) |
| `_ground_state_expand(op)` | `Nothing` or `(g, n_levels, site)` | only `Transition` overrides non-trivially |

Cross-type fallbacks in `operators/operators.jl` cover pairs of different concrete types — those always live on distinct sites, always commute, and never reduce — so the fallbacks return the trivial answers and `_commute_pair` errors on a cross-type call by construction. This is what keeps adding a new `QSym` subtype to a self-contained job: define the struct, fill in the five hooks, optionally add `adjoint` and `to_numeric` methods.


## Hilbert spaces

Hilbert spaces are used at construction time only. The `ProductSpace{T}` uses a concrete `Tuple` type parameter:

```julia
struct ProductSpace{T <: Tuple{Vararg{HilbertSpace}}} <: HilbertSpace
    spaces::T
end
```

This gives full type stability — `ProductSpace{Tuple{FockSpace, NLevelSpace}}` is a concrete type. The `⊗` operator flattens nested `ProductSpace`s during construction. Indexed families of identical subsystems are represented through the [`Index`](@ref) and [`Σ`](@ref) machinery rather than as a dedicated Hilbert-space type.


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

**Dictionary keys are full constrained terms.** Each `QTerm` is an ordered operator product plus the pairwise inequality constraints that scope that product. The ordering is canonical (established by `_partial_sort!` and the streaming passes described below), so structurally equal products always have the same `ops`. Like-term collection happens only when both `ops` and `ne` match exactly. This is the key invariant that makes constrained sums correct: `a_i` and `a_i` under `i ≠ j` are distinct stored terms, not one merged term with unioned metadata.

**The schema is visible at the type level.** `QTermDict` is a plain `Dict{QTerm, CNum}` alias — there is no wrapper struct. Iterating a `QAdd` yields `Pair{QTerm, CNum}`, and callers reach `term.ops` / `term.ne` on the key directly. This keeps the storage shape (`(ops, ne) => coeff`, not `ops => (coeff, ne)`) honest at every callsite.

**Zero terms are never stored.** `_addto!` (the single insertion helper used by every code path that adds to a `QTermDict`) builds the key, looks it up, sums the coefficient if it already exists, and deletes the entry whenever the result is zero. The dictionary therefore stays clean without an explicit cleanup pass.

**Summation metadata.** `indices` remains the sum-level metadata on `QAdd`: a `QAdd` with `indices = [i]` represents ``\sum_i (\text{terms})``. Pairwise inequality constraints like `(i, j)` meaning ``i \neq j`` live on the individual `QTerm`s, not on `QAdd` globally. This is why display and round-trip logic can represent mixed-scope sums truthfully.


## Operator sorting

```julia
@enum SiteCmp::UInt8 Less Equal Undetermined Greater

_site_compare(a::QSym, b::QSym, ne::Vector{NonEqualPair})::SiteCmp
_partial_sort!(ops::Vector{QSym}, ne::Vector{NonEqualPair})
```

Sorting is driven by a three-way comparator. Two operators have one of three site relationships: **distinct** (different `space_index`, or different numeric indices on the same space, or a `ne` entry that resolves them), **equal** (same `space_index` and syntactically equal `index`), or **undetermined** (same `space_index`, indices not syntactically equal, no `ne` entry resolves them). `_site_compare` returns `Less` or `Greater` for distinct sites, `Equal` for same-site pairs, and `Undetermined` for the third case.

**Why `Equal` and `Undetermined` are different enum values** even though both mean "do not reorder": distinguishing them at the type level removes an entire class of bugs where an operator-type comparator returns the wrong "neutral" answer. `Equal` signals to the sibling passes that the pair is a candidate for same-site composition or commutation; `Undetermined` signals that the algebra has no information yet and the pair must be left alone until something — a sum-index substitution, an explicit `≠` declaration — turns it into one of the other three.

**Why a stable insertion sort.** `_partial_sort!` is the only function in the package that ever reorders `ops`. It swaps adjacent pairs only when `_site_compare` returns `Greater`; `Equal` and `Undetermined` are left in their incoming physical order (their order encodes non-commutative multiplication that the sibling passes will interpret). Insertion sort is stable by construction and `O(n)` on the near-sorted inputs the algebra produces in practice.

For example, `[a_fock, σ²¹_i, b_fock, σ¹²_j]` with `i, j` symbolic indices on the same atom space and no `ne` constraint resolving them partial-sorts to `[a_fock, b_fock, σ²¹_i, σ¹²_j]`: the two Fock operators move to the front (they are distinct from the σs), the two σs stay in their physical order (undetermined relative to each other).


## The canonicalization pipeline

Algebraic rewrites fall into three categories, and the package treats each one differently. *Unconditionally safe* rewrites — local reductions on operators provably on the same site, such as `|i⟩⟨j|·|k⟩⟨l| = δⱼₖ |i⟩⟨l|` or the Pauli product rule — fire eagerly inside every `*`. *Conditionally safe* rewrites — commutation swaps that reorder a same-site pair like `a·a† → a†·a + 1` — fire when the site relationship is known, and stay dormant when it is undetermined. *Structurally destructive* rewrites — the completeness identity `σᵍᵍ = 1 - Σ_{k≠g} σᵏᵏ`, which is mathematically exact but multiplies one atom into `n_levels` terms — fire only on user request. This split is why the canonical-form rules below distinguish `Equal` from `Undetermined`, and why ground-state expansion has its own public function.

Every operation that perturbs an operator sequence ends in the same primitive, `_canonicalize!`, which runs `_partial_sort!` followed by a streaming pipeline of small passes:

```julia
@inline function _stream!(out, ops, c, ne)
    _reduce_ops(ops, c) do ops1, c1
        _commute_ops(ops1, c1) do ops2, c2
            _reduce_ops(ops2, c2) do ops3, c3
                _canonicalize_to_dict!(out, ops3, c3, ne)
            end
        end
    end
end
```

Each pass has signature `(ops, c, sink::F) where {F}`. The type parameter on the sink forces specialization so the nested `do`-blocks inline into one fused function with zero closure allocation — the streaming pattern reads like a deferred pipeline but at the machine level is one straight-line function.

The order `reduce → commute → reduce` matters. The first reduce folds Transition and Pauli same-site pairs, which never commute but always compose locally; this leaves only Fock, Spin, and PhaseSpace ladder pairs for `_commute_ops` to act on. The trailing reduce catches any same-site composition that surfaces when a commute residual lands next to another operator on the same site — for example, a Spin commutator's contracted-axis residual meeting a same-site neighbor.

The four passes have one job each:

| Pass | Behaviour |
|---|---|
| `_reduce_ops` | folds adjacent provably-same-site pairs via `_reduce_pair`; one output per input |
| `_commute_ops` | applies swaps for adjacent same-site pairs whose `_can_commute` is false; emits both swap branch and residual branch |
| `_expand_gs_ops` | applies `σᵍᵍ → 1 - Σ_{k≠g} σᵏᵏ` to every ground-state projector; opt-in, not part of the default pipeline |
| `_substitute_ops` | walks `ops` applying a substitution dict; forks when a value is a multi-term `QAdd` |

The terminal sink, `_canonicalize_to_dict!`, is the only place a `QTerm` is constructed during a pipeline run. It builds `QTerm(ops, ne)`, looks up in `out`, sums the coefficient, drops zero entries.

The aliasing rules are what justify the streaming style: there is no intermediate term-list materialization. Single-output passes mutate `ops` in place and the sink receives the same `Vector`; forking passes call `copy(ops)` before mutating each branch except possibly one, so no two branches share a mutable operator string; the terminal sink takes ownership of whatever it receives. Every pass documents which rule it follows.

The pipeline establishes the **canonical-form invariant**: `ops` is in partial-canonical order, no adjacent provably-same-site pair has a remaining commutator residual or reduction, and like terms are collected. The invariant deliberately says nothing about `σᵍᵍ` — see the next section.


## Simplification vs. normal ordering

The package exposes four entry points that combine the pipeline above in different ways.

`normal_order(expr)` re-streams each term through `_stream!`. Eager `*` already produces canonical form, so on the output of `*` it is idempotent; it earns its keep as a finalizer for hand-constructed expressions and as the second half of `simplify`. `simplify(expr)` runs `normal_order` and then walks the resulting terms once more, applying `Symbolics.simplify` to each coefficient and dropping summation indices that no surviving term depends on. The expensive per-coefficient simplification deliberately lives in this outer pass rather than inside the streaming pipeline: it runs once per surviving term, not on every dict insertion. `expand(expr)` distributes symbolic prefactors only — `(a + b)² → a² + 2ab + b²` — and leaves the operator structure untouched. `expand_completeness(expr)` applies the ground-state identity, described below.

`commutator(a, b)` is `a*b - b*a` in the general case, with a fast path on `QSym × QSym`: when exactly one direction of the pair is non-canonical, the commutator equals the swap residual returned by `_commute_pair`, so the call short-circuits to a single-term `QAdd` without running the full pipeline twice and a subtraction.

**Why `σᵍᵍ` stays atomic in canonical form.** Ground-state projectors are legitimate atoms in the canonical basis. The completeness identity `σᵍᵍ = 1 - Σ_{k≠g} σᵏᵏ` is mathematically exact, but applying it eagerly multiplies every product containing a `σᵍᵍ` by `n_levels`. A product with `k` ground-state factors balloons to `n_levels^k` terms before the surrounding context has a chance to cancel anything, and like-term collection across operations cannot recover the original compactness once the explosion has happened. Keeping `σᵍᵍ` atomic also serves downstream consumers: in mean-field expansions, `⟨σᵍᵍ⟩` is a single moment that solvers carry through their equations directly, while `1 - Σ ⟨σᵏᵏ⟩` is a sum of `n_levels - 1` moments tied together by an identity constraint that meanfield code would otherwise have to recognize and dedupe.

`expand_completeness(expr)` is the explicit handle for the cases where the expansion is genuinely wanted — converting to a basis in which `σᵍᵍ` is a dependent quantity, or feeding into code that expects the identity already materialized. Internally it walks each term through `_expand_gs_ops` (which forks `n_levels` ways per `σᵍᵍ`) and re-streams the output. User-constructed `σᵍᵍ` is the matching boundary case: `Transition(h, :σ, g, g)` returns a plain `Transition`, and canonicalization only happens when the operator participates in a `*`, so direct construction and inspection stay cheap.


## Diagonal splitting

When a multiplication crosses a summed index, the resulting product can equate two indices that were free in the operands. `_accumulate_with_diag!` handles that boundary case in one pass, invoked from every `*(QAdd, ·)` overload after the per-term `vcat` of operand operators.

For each summation index `sum_idx` that the operand depends on, and each free index `ext_idx` on the same Hilbert subspace that is not already constrained `≠ sum_idx`, the pass emits two contributions. The **off-diagonal** branch canonicalizes the original operator vcat under `ne` augmented with `(sum_idx, ext_idx)`; the new constraint upgrades `_site_compare` for that pair from `Undetermined` to `Less` or `Greater`, which lets `_partial_sort!` place them deterministically. The **diagonal** branch substitutes `sum_idx → ext_idx` in the unsorted operator vcat and in the prefactor, then canonicalizes under `ne` with any constraint involving `sum_idx` dropped. The two branches together cover the partition `sum_idx ≠ ext_idx` and `sum_idx = ext_idx` exactly once each.

The one subtle point is that the diagonal substitution operates on the **unsorted** vcat. `_partial_sort!` can place same-space operators with different symbolic indices in an order that, after the substitution makes those indices equal, would feed `_reduce_ops` a same-site composition pair in the wrong physical order. Concretely, ``\sigma_j^{12} \cdot \sigma_i^{21}`` (with `i ≠ j`) partial-sorts on `name` order, and substituting `i → j` after that sort produces a sequence whose same-site collapse is ``\sigma_j^{22}``. Substituting on the unsorted product gives ``\sigma_j^{12} \cdot \sigma_j^{21} = \sigma_j^{11}``, the physically correct result. Doing the substitution before sorting preserves the physical adjacency that the same-site rules need.

This mechanism only handles the case where one of the two indices is bound by a `Σ`. Two free indices outside any sum land in the `Undetermined` regime and need an explicit user declaration — see the next section.


## Free indices and `assume_distinct_index`

Two operators with different symbolic indices on the same Hilbert subspace, neither bound by a `Σ`, have an undetermined site relationship. The algebra cannot tell whether the user means "these label distinct atomic sites" or "these are two index variables that may or may not coincide", and the conservative reading is the second: the operators stay in their physical order, no same-site collapse fires, and the resulting expression carries the ambiguity faithfully.

`assume_distinct_index(q, [(α, β), …])` resolves the ambiguity in the first direction: it augments every term's `ne` with the supplied pairs, re-canonicalizes so `_partial_sort!` can place the resolved pairs deterministically, and runs `expand_completeness` so any ground-state projectors that emerge from same-site composition under the new constraint are folded. The two-atom inter-atom coherence `σⱼ¹² · σₖ²¹` is canonicalized by `assume_distinct_index(σⱼ¹² · σₖ²¹, [(j, k)])`.


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

Symbolics.jl conventionally wraps public-API outputs in `Symbolics.Num`, a `<: Real` adapter over the raw `SymbolicUtils.BasicSymbolic` tree. `average` deliberately breaks that convention because `Num` wrapping is not type-stable across `average(::QAdd)`: pulling c-number prefactors out of the average routes the internal arithmetic through SymbolicUtils promotion rules, which produce `BasicSymbolic`, and a `QAdd` that contains a symbolic-prefactor term then disagrees in wrapper type with a `QAdd` that does not.

Returning `BasicSymbolic` uniformly across `average(::QSym)`, `average(::QAdd)`, and `average(::BasicSymbolic)` (scalars pass through unchanged) keeps the return type stable, lets every downstream function (`is_average`, `acts_on`, `get_indices`, `numeric_average`, `undo_average`) dispatch on one type without unwrapping, and avoids fragile re-wrapping. `Num`-accepting dispatches exist as defensive fallbacks for callers (or Symbolics arithmetic) that re-wrap.

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


## Hermitian conjugation (operators.jl)

Two exported helpers handle Hermitian conjugation on mixed operator/symbolic expression trees, including averaged operator nodes (`⟨...⟩`):

- **`qadjoint(x)`** (aliased as `qconj` and `dagger`): Hermitian conjugate that distributes through `SymbolicUtils.BasicSymbolic` trees and dispatches to `adjoint` on `QField` and `Number`. Distinct from `Base.conj`, which on a `BasicSymbolic` returns an opaque `conj(...)` wrapper instead of recursing into arguments; the distributed form is needed by downstream hashing and substitution machinery.
- **`inner_adjoint(x)`**: Pushes the adjoint *inside* `AvgFunc` nodes, rewriting `conj(⟨X⟩)` as `⟨X†⟩`. Used when building equations of motion where both sides must share the canonical "average-of-operator" form. Also collapses nested `conj(avg(...))` by recursing into the argument.

These cannot be wired directly to `Base.conj`/`Base.adjoint` on `SymbolicUtils.BasicSymbolic` because that would be type piracy — defining methods on `Base` functions for a type we don't own. Downstream packages (QuantumCumulants.jl, QuantumInputOutput.jl) call these explicitly.


## Printing and LaTeX

**Terminal printing** uses Unicode: `†` for dagger, subscript digits (`₀`-`₉`) for Transition levels, `σx`/`σy`/`σz` for Pauli axes. Summations render as `Σ(i=1:N)`.

**`sorted_arguments`** ensures deterministic output order. The sort key is `(length(ops), full_op_keys...)` where `_full_op_key(op) = (_sort_key(op)..., _type_order(op), op.name)`. This gives: shorter terms first, then by site, then by type (Destroy < Create < Transition < Pauli < Spin < Position < Momentum), then by name.

**LaTeX** uses Latexify.jl's `@latexrecipe` macro. `transition_superscript(::Bool)` toggles the global `transition_idx_script` `Ref` between `:^` and `:_`, controlling whether Transition level indices render as superscripts (`{name}^{{ij}}`) or subscripts (`{name}_{{ij}}`).
