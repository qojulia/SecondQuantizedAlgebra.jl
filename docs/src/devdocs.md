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

A sixth, defaulted hook supports the reduce pass: `_may_reduce(a, b)::Bool` answers whether `_reduce_pair(a, b)` could return anything other than `NoReduction`. The default is `false`; only the same-type pairs that genuinely compose (`Transition×Transition`, `Pauli×Pauli`) override it to `true`. The reduce pass consults this isbits-`Bool` gate before the `_reduce_pair` call. The reason is type stability: `_reduce_pair` returns a `(kind, op, factor)` tuple holding a non-isbits `QSym`, so a dynamic dispatch on an abstract `QSym` pair boxes that tuple on every adjacent pair, whereas the `_may_reduce` `Bool` comes back unboxed. Fock/Spin/PhaseSpace products (which never reduce) therefore skip the boxing entirely. A new reducing operator type overrides both `_reduce_pair` and `_may_reduce`; a new non-reducing type needs neither.

Cross-type fallbacks in `operators/operators.jl` cover pairs of different concrete types (always distinct sites, always commuting, never reducing), so the fallbacks return the trivial answers and `_commute_pair` errors on a cross-type call by construction. This is what keeps adding a new `QSym` subtype to a self-contained job: define the struct, fill in the five hooks, optionally add `adjoint` and `to_numeric` methods.


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

**Opaque vs. split complex parameters.** `_to_cnum` only splits a parameter into `real`/`imag` parts when it already arrives as a `Complex{Num}` (a `::Complex` variable, or an explicit `complex(re, im)` node). Any other symbol, including a `Number`-symtype variable (`@variables η::Number`), is stored *opaquely* in the real slot as `Complex(Num(η), 0)`. This keeps coefficient arithmetic on a single symbol (`η * η → η²`, one multiply) instead of expanding `(a+bi)(c+di)` over two independent unknowns, the cost a `::Complex` parameter pays on every product.

The opaque storage means conjugation cannot be the generic `conj(::Complex{Num})`, which only flips `.im` and is a no-op when `.im == 0`; otherwise the phase of a `Number`-symtype parameter would be silently dropped. `_conj_cnum` is therefore **symtype-aware**: it leaves a real-symtype real part unchanged (its own conjugate) and applies the symbolic `conj` to a non-real one.

```julia
_sym_conj(x::Num) = SymbolicUtils.symtype(x) <: Real ? x : Num(conj(SymbolicUtils.unwrap(x)))
_conj_cnum(c::CNum) = Complex(_sym_conj(c.re), -c.im)
```

`Base.adjoint(::QAdd)` conjugates each term coefficient through `_conj_cnum`, so `η::Number` correctly satisfies `adjoint(η) == conj(η)`, and `(η a)†(η a)` carries `|η|² = conj(η)·η` rather than `η²`.


## QAdd internals

```julia
const NonEqualPair = Tuple{Index, Index}

struct QTerm
    ops::Vector{QSym}
    ne::Vector{NonEqualPair}
    hash::UInt          # cached hash(ops, ne), computed once at construction
end

const QTermDict = Dict{QTerm, CNum}

struct QAdd <: QField
    arguments::QTermDict
    indices::Vector{Index}
end
```

**Why `QTerm` caches its hash.** `QTerm` is a dict key, and every `_addto_key!` both probes (`get`) and writes (`setindex!`/`delete!`) the same key, so the key is hashed at least twice per insertion, plus once more for each entry whenever the dict grows and rehashes. Hashing recurses over the whole `ops` vector, and because `ops` has the abstract eltype `QSym`, each per-operator `hash` is a dynamic dispatch (its `UInt` result is boxed). Caching `hash(ops, ne)` once in the inner constructor turns every later hash of that key into a single `hash(::UInt, h)` and makes dict-growth rehashing free. The cache is sound because `QTerm` is immutable and `ne` is already canonicalized before construction, so structurally equal keys always carry the same cached value; `isequal` short-circuits on the cached hash before comparing `ops`. A trusted three-argument constructor lets `_copy_key` carry the known hash across a verbatim copy without recomputing it.

**What a dict entry represents.** A single entry `QTerm(ops, ne) => c` in `arguments` represents the term `c · ops[1] · ops[2] · …` valid for any index assignment satisfying the pairwise constraints in `ne` (each `(α, β) ∈ ne` means `α ≠ β`). The `indices` field on `QAdd` carries the outer summation scope: a `QAdd` with `indices = [i, j]` represents ``\sum_i \sum_j \sum_\text{terms} c \cdot \text{ops}`` where each individual term may further constrain `(i, j)` per its own `ne`. This per-term scoping is what lets a single `QAdd` represent expressions like ``\sum_i a_i a_i + \sum_{i \neq j} a_i a_j`` as two dict entries with different `ne` rather than two separate `QAdd` summations.

**Dictionary keys are full constrained terms.** Each `QTerm` is an ordered operator product plus the pairwise inequality constraints that scope that product. The ordering is canonical (established by `_partial_sort!` and the streaming passes described below), so structurally equal products always have the same `ops`. Like-term collection happens only when both `ops` and `ne` match exactly. This is the key invariant that makes constrained sums correct: `a_i` and `a_i` under `i ≠ j` are distinct stored terms, not one merged term with unioned metadata.

**The schema is visible at the type level.** `QTermDict` is a plain `Dict{QTerm, CNum}` alias — there is no wrapper struct. Iterating a `QAdd` yields `Pair{QTerm, CNum}`, and callers reach `term.ops` / `term.ne` on the key directly. This keeps the storage shape (`(ops, ne) => coeff`, not `ops => (coeff, ne)`) honest at every callsite.

**Zero terms are never stored.** `_addto!` (the single insertion helper used by every code path that adds to a `QTermDict`) builds the key, looks it up, sums the coefficient if it already exists, and deletes the entry whenever the result is zero. The dictionary therefore stays clean without an explicit cleanup pass.

**Summation metadata.** `indices` remains the sum-level metadata on `QAdd`: a `QAdd` with `indices = [i]` represents ``\sum_i (\text{terms})``. Pairwise inequality constraints like `(i, j)` meaning ``i \neq j`` live on the individual `QTerm`s, not on `QAdd` globally. This is why display and round-trip logic can represent mixed-scope sums truthfully.

**Dead-NE invariant.** Every `QTerm.ne` pair must reference at least one index that the term can observe: either some operator in `term.ops` carries that index, the coefficient `c::CNum` depends on it (e.g. an `IndexedVariable` factor), or it appears in the enclosing `QAdd.indices` (sum scope). Pairs failing all three tests encode nothing observable, would only obstruct dict dedup, and so must not survive storage. The `QAdd` inner constructor enforces this by calling `_prune_dead_ne(arguments, indices)`: it walks every term, drops the dead pairs through `_addto_key!` (which performs the like-term merge if dead-stripped keys collide), and is idempotent on already-clean input (no rebuild, no allocation). Algebraic code that constructs `QAdd`s does not need to clean up NE manually; the constructor does it. The predicate is `_depends_on_index_term`, the same one used everywhere else NE/scope dependence is queried, so the rule stays consistent with `_canonicalize!` and `_emit_scaled_by_scope!`.


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


## Disjoint bound indices in products

Diagonal splitting handles a *bound* index in one factor meeting a *free* index in the other. Two bound indices that share a display name are something else. The expression

```math
\left(\sum_i X_i\right)\left(\sum_i Y_i\right)
```

is ambiguous as written. Two readings are equally available:

1. **Alpha-rename one side.** The two `i` letters denote independently bound variables that happen to print the same; the product is ``\sum_{i,j} X_i Y_j``. This is the convention math papers use, but only after the reader silently introduces a fresh variable the writer did not write.
2. **Share the bound variable.** The two occurrences denote one variable, giving the diagonal ``\sum_i X_i Y_i``, a different operator.

The reader applies alpha-conversion on the fly; the algebra cannot. The user constructed exactly one `Index(:i, …)` object and used it on both sides, so there are no two implicit `i`s in the term graph for the algebra to discover. Reading 2 would silently drop the off-diagonal ``\sum_{i \neq j}`` contributions; reading 1 would require inventing a fresh `Index` the user can't reach (see *Naming policy*).

`*(QAdd, QAdd)` therefore throws `ArgumentError` when `a.indices ∩ b.indices` is non-empty. The caller disambiguates: rename one side via `change_index` with a user-constructed `Index`, or build the two factors with distinct `Index` objects from the start.


## Free indices and `assume_distinct_index`

Two operators with different symbolic indices on the same Hilbert subspace, neither bound by a `Σ`, have an undetermined site relationship. The algebra cannot tell whether the user means "these label distinct atomic sites" or "these are two index variables that may or may not coincide", and the conservative reading is the second: the operators stay in their physical order, no same-site collapse fires, and the resulting expression carries the ambiguity faithfully.

`assume_distinct_index(q, [(α, β), …])` resolves the ambiguity in the first direction: it augments every term's `ne` with the supplied pairs, re-canonicalizes so `_partial_sort!` can place the resolved pairs deterministically, and runs `expand_completeness` so any ground-state projectors that emerge from same-site composition under the new constraint are folded. The two-atom inter-atom coherence `σⱼ¹² · σₖ²¹` is canonicalized by `assume_distinct_index(σⱼ¹² · σₖ²¹, [(j, k)])`.


## Naming policy

Indices are user-owned: the algebra never mints `Index` objects on the user's behalf. Every name appearing in any output traces back to a user `Index(...)` call. The principle is operational, not aesthetic; an algebra-invented `Index` is invisible to the user's vocabulary, breaks pattern-matching on equation outputs, and gives no handle for `evaluate(...; limits = ...)` or initial-condition substitution.

Three consequences beyond the `*(QAdd, QAdd)` throw described above:

- **Diagonal splitting only fires for `(sum_idx, ext_idx)` pairs where `ext_idx` is already free in the operand.** The algebra does not invent a fresh `ext_idx`; without a user-declared free index on the same space, the `(i = ext_idx)` branch has nothing to substitute into.
- **`assume_distinct_index(q, [(α, β), …])` is the user's channel** for resolving `Undetermined` free pairs. The user supplies the inequality; the algebra applies it.
- **No public helper renames bound variables.** Consumers needing alpha-rename use `change_index` with their own freshly-constructed `Index`. QuantumCumulants.jl's `complete!`, for example, mints completion-internal canonical names from the user's existing vocabulary.
- **Alpha-equivalent sums are not auto-collected.** ``\sum_i \sigma_i + \sum_j \sigma_j`` stays two terms rather than folding to ``2 \sum_i \sigma_i``, even though the two summations are mathematically identical. Collection keys on the exact `(ops, ne)` term identity, and `σ_i` and `σ_j` are distinct keys; merging them would require the algebra to alpha-rename one bound index onto the other, which it will not do. The stored form is numerically correct. To collapse it, rename one side with a user-constructed `Index`: `change_index(Σ(σ_j, j), j, i) + Σ(σ_i, i)` gives ``2 \sum_i \sigma_i``. This is distinct from the *product* case above: addition carries no ``\sum_{i,j}``-vs-``\sum_i`` ambiguity, so the non-merge is purely a naming-policy consequence, not a semantic one.


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

**Indexed sums become a dedicated node.** When averaging an indexed `QAdd`, each index-dependent term is wrapped in a moment-layer sum node `sym_sum(body, scope)` (operation `SumFunc`, mirroring the `AvgFunc` average node) rather than stamped with metadata. The `body` is the term's averaged scalar (numeric prefactor and multi-factor c-number coefficient included), so an index-dependent coefficient stays *inside* the sum scope; the previous metadata-on-leaf form could not express ``\sum_i (\text{numeric} \cdot f_i \cdot g_i)`` because `Add`/`Mul` discard metadata on any composite or numerically-scaled node. Terms sharing the same `(indices, ne)` are grouped into one node with a multi-term body, mirroring display grouping (`_group_dep_terms`).

The scope rides as a `SumScope` *argument* of the `Term`, not as metadata, because SymbolicUtils `isequal`/`hash` ignore metadata: two sums over the same body but different scope would otherwise be `isequal` and wrongly cancel in a subtraction. `SumScope` defines value-based `==`/`isequal`/`hash`; SymbolicUtils stores it as an opaque scalar argument, the same way an `AvgFunc` node stores its `QField` (a bare `Vector` argument would instead be array-ified into a symbolic array, losing the `Index` objects). `is_indexed_sum` recognises the node, `get_sum_indices`/`get_sum_non_equal` (and the retained alias `has_sum_metadata`) read the scope, and `undo_average` rebuilds the summed `QAdd` by folding the node's `ne` into each term and attaching its `indices`. Every consumer that walks averaged trees (`get_indices`, `acts_on`, `make_time_dependent`, `inner_adjoint`, `numeric_average`) recurses into `body` and reconstructs the node around the rewritten body, so the scope is preserved without special-casing each operation.

**No `Complex{Num}` intermediates in the result.** `average(::QAdd)` walks `(term, c::CNum)` pairs and folds each one into a single result. It must never let `result` become a `Complex{Num}`, because `SymbolicUtils.unwrap(::Complex{<:Num})` materialises a literal `Term(complex, [re, img])` node whenever both `re` and `img` are symbolic, and that node is opaque to `simplify` / `expand` and also generates a runtime `complex(::Real, ::Complex)` call when consumers like ModelingToolkit codegen the equation (Base defines no such method). Each per-term contribution is unwrapped to `BasicSymbolic{SymReal}` (the single concrete Moshi type backing `BasicSymbolic`) at the boundary, so the accumulator and the grouped sum bodies are concretely typed rather than `Any`: every contribution reduces to `SymReal` because the averaged leaf, `Symbolics.IM`, and the `CNum` real/imaginary parts are all real-symtype. The imaginary unit enters the chain through `Symbolics.IM`, the `BasicSymbolic{SymReal}` symbol that prints as `im`, instead of through `1im::Complex{Bool}`. Concretely, the constant branch contributes `r + i * Symbolics.IM` and the operator branch contributes `r * avg + i * Symbolics.IM * avg`, skipping the `r` or `i` add when zero. This keeps every result a clean polynomial in `Symbolics.IM` that downstream `simplify(...; expand = true)` can reduce to zero on identically-equal differences.

**`undo_average` always returns `QAdd`.** It accepts both `BasicSymbolic` and `Num` inputs (unwrapping `Num` first). Scalars become single-term `QAdd`s with an empty operator sequence, and lone `QSym`s become single-term `QAdd`s with unit prefactor. This uniform return type makes `undo_average` type-stable despite walking SymbolicUtils expression trees (where `operation(x)` returns `Any`).

**Singleton single-op terms average to the bare op.** `average(::QAdd)` wraps a term in `_single_qadd(_CNUM_ONE, ops, ne)` only when cross-op canonicalization needs the constraint: `length(term.ops) > 1`. A single-op term averages to `only(term.ops)` regardless of `term.ne`, because the enclosing scope now rides on the `sym_sum` node rather than on the average leaf. This relies on the dead-NE invariant from *QAdd internals*: any `ne` reaching `average` references some op index or sum-scope index, but for a single-op term there is no second op to constrain, so a leaf-level `ne` would be observationally inert and would only obstruct downstream dict lookups (notably QC's MTK variable map).


## Numeric conversion

`to_numeric` maps symbolic operators to QuantumOpticsBase matrices:

- **Simple basis:** Direct dispatch — `Destroy → destroy(b)`, `Transition → transition(b, i, j)`, etc.
- **Composite basis:** Embeds the single-site matrix as a `LazyTensor`:
  ```julia
  op_num = to_numeric(op, b.bases[idx])
  LazyTensor(b, [idx], (op_num,))
  ```

**`_to_number`** extracts plain Julia numbers from `Num`/`CNum` wrappers for numeric evaluation. Falls back to the symbolic value if it can't be unwrapped (for symbolic prefactors that haven't been substituted yet).

**`_reduce_const`** reduces a fully-substituted coefficient part to a number. `Symbolics.value` handles numeric constants directly, but a part it leaves symbolic (for example `conj` of a complex literal, which SymbolicUtils does not fold) is compiled with `build_function` and evaluated. Because a `Number`-symtype parameter is held opaquely in the real slot, that real part can itself reduce to a `Complex`, so `_to_complex` recombines as `_reduce_const(real(x)) + im * _reduce_const(imag(x))` rather than `Complex(re, im)`.

**`_lazy_one`** creates the identity operator. For simple bases it returns `one(b)` (dense identity). For composite bases it returns a `LazyTensor` identity rather than materializing the full Kronecker-product identity matrix.


## Hermitian conjugation (operators.jl)

Two exported helpers handle Hermitian conjugation on mixed operator/symbolic expression trees, including averaged operator nodes (`⟨...⟩`):

- **`qadjoint(x)`** (aliased as `qconj` and `dagger`): Hermitian conjugate that distributes through `SymbolicUtils.BasicSymbolic` trees and dispatches to `adjoint` on `QField` and `Number`. Distinct from `Base.conj`, which on a `BasicSymbolic` returns an opaque `conj(...)` wrapper instead of recursing into arguments; the distributed form is needed by downstream hashing and substitution machinery.
- **`inner_adjoint(x)`**: Pushes the adjoint *inside* `AvgFunc` nodes, rewriting `conj(⟨X⟩)` as `⟨X†⟩`. Used when building equations of motion where both sides must share the canonical "average-of-operator" form. Also collapses nested `conj(avg(...))` by recursing into the argument.

These cannot be wired directly to `Base.conj`/`Base.adjoint` on `SymbolicUtils.BasicSymbolic` because that would be type piracy — defining methods on `Base` functions for a type we don't own. Downstream packages (QuantumCumulants.jl, QuantumInputOutput.jl) call these explicitly.


## Ordering keys: three distinct orders

SQA carries three orderings that must not be conflated, because each answers a different question.

**`_site_compare` (partial, commutation).** The `SiteCmp` three-way comparator drives `_partial_sort!` during normal ordering. It is deliberately *not* identity-faithful: it returns `Equal` for operators that differ only in `axis` (Pauli/Spin) or levels `i,j` (Transition), because it encodes which adjacent factors may be reordered, not whether two operators are the same. It is also partial (`Undetermined` for free indices with no resolving `ne`).

**`_full_op_key`/`_sort_key` (display sort).** Used by `sorted_arguments` for deterministic *display* order only. Also not identity-faithful: it omits axis and levels, so two distinct operators can share a display key. Fine for printing, wrong for keying.

**`order_key`/`term_order_key`/`qadd_order_key` (total, identity-faithful).** The public ordering used by downstream packages (QuantumCumulants.jl) to pick canonical representatives and compare expressions reproducibly. The contract is `order_key(a) == order_key(b)` iff `isequal(a, b)`: the key carries every identity field, so it is a strict total order that never ties distinct operators. `order_key` has one method per concrete `QSym` subtype, co-located with that type's `isequal`/`hash`, and no generic `QSym` fallback, so a new operator type added without a key is a loud `MethodError` rather than a silent field-dropping collision. It is preferred over hashing for ordering because Julia's `hash` is not stable across versions/platforms (which would make canonical-representative choice, and therefore generated equations, irreproducible) and collides; hashing remains correct for equality-only dedup via `Dict`/`Set`, which use `hash` plus `isequal`.

## Printing and LaTeX

**Terminal printing** uses Unicode: `†` for dagger, subscript digits (`₀`-`₉`) for Transition levels, `σx`/`σy`/`σz` for Pauli axes. Summations render as `Σ(i=1:N)`.

**`sorted_arguments`** ensures deterministic output order. The sort key is `(length(ops), full_op_keys...)` where `_full_op_key(op) = (_sort_key(op)..., _type_order(op), op.name)`. This gives: shorter terms first, then by site, then by type (Destroy < Create < Transition < Pauli < Spin < Position < Momentum), then by name.

**LaTeX** uses Latexify.jl's `@latexrecipe` macro. `transition_superscript(::Bool)` toggles the global `transition_idx_script` `Ref` between `:^` and `:_`, controlling whether Transition level indices render as superscripts (`{name}^{{ij}}`) or subscripts (`{name}_{{ij}}`).
