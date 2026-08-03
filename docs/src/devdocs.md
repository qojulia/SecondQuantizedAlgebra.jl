# Developer Documentation

This page explains the internal architecture and design rationale of **SecondQuantizedAlgebra.jl**. It is intended for contributors and anyone who wants to understand *why* things are implemented the way they are.

## Type hierarchy

```
QField (abstract)
├── QSym (abstract): atomic operators (leaves)
│   └── Op: the single concrete leaf; a `kind::OpKind` tag picks the role
│           Destroy / Create (FockSpace), Transition (NLevelSpace),
│           CollectiveTransition (CollectiveNLevelSpace),
│           Pauli (PauliSpace), Spin (SpinSpace), Position / Momentum (PhaseSpace)
└── QAdd: sum of QTerm products (the only compound type)
```

`QTerm` is the per-entry storage key (operator product + non-equal constraints) used as the dict key inside `QAdd`. There is no abstract `QTerm` supertype; `QAdd` is a `QField` directly.

`QSym` stays abstract with `Op` as its sole subtype so external `::QSym` signatures keep resolving, while the operator *vector* `QTerm.ops` is `Vector{Op}` with a concrete element type. That concrete eltype is the point of the collapse: the per-operator hooks below dispatch and inline statically, where the former subtype hierarchy forced a dynamic dispatch (and a boxed result) on every per-operator call.

**Why `QAdd` is the only compound type.** Earlier versions of the package had both `QMul` (products) and `QAdd` (sums). This created a two-level expression tree where dispatch needed to handle `QSym`, `QMul`, and `QAdd` at every level, and the return type of `*` was unpredictable (`QSym`, `QMul`, or `QAdd` depending on simplification). The redesign collapses this into a single `QAdd` whose internal dictionary is keyed by the full term identity `(ops, ne)` and stores only the prefactor as the value. Every multiplication immediately produces a `QAdd`, giving a uniform return type. This type stability is critical for performance — the Julia compiler can infer return types through chains of arithmetic, avoiding dynamic dispatch and heap-allocated boxes at every intermediate step. The exact-key representation also keeps like-term collection honest: only terms with the same operator string *and* the same scoped constraints are merged.


## Operator struct layout

All operators are one concrete struct with a runtime `kind` tag; the level/axis
fields are shared storage interpreted per `kind`:

```julia
@enum OpKind::UInt8 OP_DESTROY OP_CREATE OP_TRANSITION OP_PAULI OP_SPIN OP_POSITION OP_MOMENTUM OP_COLLECTIVE_TRANSITION

struct Op <: QSym
    kind::OpKind        # which physical role this operator plays
    name_id::Int32      # interned display name (intern.jl); resolve via operator_name(op)
    space_index::Int32  # which subspace in a ProductSpace (always 1 for simple spaces)
    index::Index        # symbolic summation index (NO_INDEX when absent)
    l1::Int32           # Transition/CollectiveTransition i | Pauli/Spin axis
    l2::Int32           # Transition/CollectiveTransition j
    g::Int32            # Transition ground state
    nlev::Int32         # Transition number of levels
end
```

The role names `Destroy`, `Create`, `Transition`, `CollectiveTransition`,
`Pauli`, `Spin`, `Position`, `Momentum` are constructor functions returning an
`Op` with the matching `kind`;
`is_destroy`/`is_create`/… (and `optype`) are the exported predicates that
replace `isa`. Fock/Position/Momentum leave `l1..nlev` zero, and
create-vs-destroy lives in `kind` (the old `ladder` field is gone, though
`ladder(::Op)` is still provided for canonical-ordering callers).

**`Op` is `isbits` (issue #137).** Every field is a bits type, so `Op` is
`isbits` (44 bytes) and `Vector{Op}` (the per-term `ops` storage) is a dense
inline buffer the GC never scans, hashed and compared on pure integers. The win
over a `Symbol`-name layout came from the struct *shrinking*: a 72-byte
`Symbol`-name `Op` and even an `InlineStrings.String7` 72-byte layout measured
within noise, while interning the name to an `Int32` (so the struct drops to
44 bytes) measured about 40% faster on dict-keyed product/sum construction. The
name is therefore stored as an interned `Int32` id into a module-global table
(see `intern.jl`); `operator_name(op)::Symbol` resolves it back for printing. The
role constructors take a `Symbol` name (interned via `_name_id`); internal
rebuilds (adjoint, `IndexedOperator`, reduce/commute residuals,
`expand_completeness`) forward the existing `name_id` by constructing `Op`
directly, so an already-interned id is never re-interned and no raw id is ever
accepted on the public constructor surface.

**Custom `hash`/`isequal` are kept, now trivially cheap.** They hash/compare the
`kind`, `name_id`, `space_index`, packed levels, and the `Index` (whose own
`hash`/`==` compare interned ids and short-circuit on `NO_INDEX` via the shared
`const` sentinel). They are retained for explicit field ordering and the
`Index === ` short-circuit, not to avoid `Num` recursion (no `Num`s remain in
`Index`).

**Canonical ordering uses a name-rank table, not the raw id.** Interned ids are
insertion-order, so ordering operators by `name_id` would make canonical form
depend on declaration order. `intern.jl` maintains `_NAME_RANK` (each id's
lexicographic position, recomputed when a new name is interned, a cold path), and
every ordering site (`_site_compare`, `order_key`, `_full_op_key`, `_sort_key`,
the diagonal-pair filter) orders by `_name_rank(name_id)`. This preserves the
exact alphabetical canonical form and keeps it deterministic across sessions even
though the ids themselves are not.

**Why `space_index` instead of storing the Hilbert space.** Operators don't hold a reference to their Hilbert space. The space is only used at construction time for validation (checking bounds, matching types). At runtime, only the integer `space_index` matters — it determines which operators commute (different `space_index` → commute trivially) and where to embed in a composite basis during numeric conversion. This keeps operators lightweight and avoids type instability from heterogeneous space references.

**Why `index::Index` on every operator.** Symbolic summation indices (for expressions like ``\sum_i a_i^\dagger a_i``) live directly on the operator rather than in a wrapper type. This avoids a separate `IndexedOperator` struct in the type hierarchy and keeps dispatch simple. The sentinel `NO_INDEX` (with `space_index = 0`) indicates no index is present, checked via `has_index(idx) = idx.space_index != 0`. The `_same_site` check requires both `space_index` and `index` to match — operators with the same `space_index` but different `index` represent distinct sites in an indexed sum and don't interact via commutation rules.

**Per-role field meaning.** `Transition` reads the bra/ket levels from `l1`/`l2`
and the ground state and level count from `g`/`nlev`; `Pauli` and `Spin` read the
axis (1=x, 2=y, 3=z) from `l1`. `CollectiveTransition` also reads `l1`/`l2`,
but leaves `g`/`nlev` zero because it has no completeness expansion. Fock and
PhaseSpace operators ignore the packed fields entirely.

**The five operator hooks.** The whole algebra talks to operators through five
methods. After the collapse each is a single `(::Op, ::Op)` method that branches
on `kind` (in `operators/operators.jl`), rather than the former
one-method-per-subtype-pair. They are the entire interface; everything else
builds on them.

| Hook | Returns | Meaning |
|---|---|---|
| `_site_compare(a, b, ne)` | `SiteCmp` | three-way site comparison driving the sort |
| `_can_commute(a, b)` | `Bool` | true iff swapping needs no commutator residual (called only on provably-same-site pairs) |
| `_commute_pair(a, b)` | `(swap_b, swap_a, c1, ops1, c2, ops2)` | swap and up to two residuals for same-site non-commuting pairs |
| `_reduce_pair(a, b)` | `(ReduceKind, Op, CNum)` | local algebraic identity (Transition composition, Pauli product, …) |
| `_ground_state_expand(op)` | `(is_gs, g, n_levels, site)` | only `Transition` (a `σᵍᵍ`) returns non-trivially |

Because `Op` is concrete, these hooks now infer concrete return types: `_commute_pair` and `_reduce_pair` return `Tuple{Op, …}` instead of `Tuple{QSym, …}`, so the former tuple boxing is gone for free.

A sixth, defaulted hook supports the reduce pass: `_may_reduce(a, b)::Bool` answers whether `_reduce_pair(a, b)` could return anything other than `NoReduction`. It is `true` only for `Transition×Transition` and `Pauli×Pauli`. With the concrete `Op` eltype its original boxing-avoidance motivation is moot, but it still cheaply skips the same-site field checks for the common non-reducing pair.

`CollectiveTransition` is the reason `_commute_pair` has two residual slots:
``[S^{ij},S^{kl}] = \delta_{jk}S^{il} - \delta_{li}S^{kj}`` can produce two
terms. Fock, PhaseSpace, and Spin populate only the first slot. The commute pass
forks an independent branch for each nonzero residual, and the `commutator`
leaf fast path sums both slots.

**Site families inside `_site_compare`.** The comparator first orders by `space_index`, so operators are grouped by the subspace they act on before anything else is consulted (operators on different subspaces commute, so this reordering is always safe and gives the more natural per-subspace canonical form). Within a subspace, cross-role comparison is family-scoped: Fock `{Destroy, Create}` and PhaseSpace `{Position, Momentum}` compare within their family (PhaseSpace ignores the operator name, treating x and p as conjugate variables on one site); the other roles are singleton families. A subspace carries a single Hilbert-space type, so same-`space_index` operators always share a family; the `kind` integer fallback preserves the existing values `OP_DESTROY=0 < … < OP_MOMENTUM=6` and appends `OP_COLLECTIVE_TRANSITION=7`. It only distinguishes the pathological case of two operators sharing a `space_index` across unrelated simple spaces.

**Adding an operator role.** Add an `OP_*` enum arm, a constructor function and an `is_*` predicate, then add any non-default branches required in the six hooks plus `adjoint`, `order_key`, the `numeric_operator` arm in each backend extension (`ext/`), and the printing/Latexify methods. The hooks are written so a future open-extension escape hatch (an `OP_CUSTOM` arm carrying a payload) could be slotted in without reshaping them; a `Union{Nothing,QSym}` payload field is deliberately *not* added now because it would fail the `CheckConcreteStructs` (`all_concrete`) gate.


## Hilbert spaces

Hilbert spaces are used at construction time only. The `ProductSpace{T}` uses a concrete `Tuple` type parameter:

```julia
struct ProductSpace{T <: Tuple{Vararg{HilbertSpace}}} <: HilbertSpace
    spaces::T
end
```

This gives full type stability — `ProductSpace{Tuple{FockSpace, NLevelSpace}}` is a concrete type. The `⊗` operator flattens nested `ProductSpace`s during construction. Indexed families of identical subsystems are represented through the [`Index`](@ref) and [`Σ`](@ref) machinery. `CollectiveNLevelSpace` is deliberately separate: it represents the permutation-symmetric collective role and converts numerically through `ManyBodyBasis`.

For two levels, `CollectiveTransition` is the transition-basis form of the
existing collective `Spin` algebra: `S^{12}=S_+`, `S^{21}=S_-`, and
`(S^{11}-S^{22})/2=S_z` in the QuantumOpticsBase convention. They remain
separate roles because they use different basis conventions and numeric paths
(`ManyBodyBasis` versus `SpinBasis`).


## The `CNum` prefactor type

```julia
const CNum = Complex{Num}
```

All prefactors are promoted to `CNum` — integers, floats, symbolic variables, complex numbers all get converted via `_to_cnum`. This ensures a single concrete prefactor type throughout the entire algebra, which is essential for type stability and therefore performance. The Julia compiler sees `Dict{QTerm, CNum}` (with `QTerm.ops::Vector{Op}`) as a fully concrete type, so dictionary operations, arithmetic, and iteration never trigger dynamic dispatch.

`Num` is the symbolic number type from Symbolics.jl, so `Complex{Num}` can represent both symbolic prefactors (`ω`, `g + im*κ`) and plain numeric ones (`3.0 + 0.0im`).

**Fast-path arithmetic.** Splitting into real and imaginary `Num` parts enables fast-path arithmetic. Most physics prefactors have zero imaginary part, so `_mul_cnum` and `_add_cnum` check `_iszero_num(c.im)` and skip half the arithmetic:

```julia
# Common case: both purely real → 1 multiply instead of 4
if ai_zero && bi_zero
    return Complex(ar * br, _NUM_ZERO)
end
```

The `_iszero_num` check skips the imaginary arithmetic, and cached constants (`_CNUM_ZERO`, `_CNUM_ONE`, `_CNUM_IM`, etc.) avoid re-constructing the common literals. There is, however, **no native-number fast path**: the real-part product `ar * br` always goes through `Num` arithmetic, even when both operands wrap plain numbers. Because Symbolics wraps numeric literals as hash-consed `BasicSymbolic` constants, every coefficient multiply/add and every `Num` construction costs on the order of 10²–10³ ns (versus ~1 ns native). On product- and power-heavy workloads this coefficient arithmetic — not the operator algebra — dominates the run time.

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
    ops::Vector{Op}
    ne::Vector{NonEqualPair}
    hash::UInt          # cached hash(ops, ne), computed once at construction
end

const QTermDict = Dict{QTerm, CNum}

struct QAdd <: QField
    arguments::QTermDict
    indices::Vector{Index}
end
```

**Why `QTerm` caches its hash.** `QTerm` is a dict key, and every `_addto_key!` both probes (`get`) and writes (`setindex!`/`delete!`) the same key, so the key is hashed at least twice per insertion, plus once more for each entry whenever the dict grows and rehashes. Hashing recurses over the whole `ops` vector (each per-operator `hash` is now a static call, since `ops::Vector{Op}` has a concrete element type). Caching `hash(ops, ne)` once in the inner constructor still pays off: it turns every later hash of that key into a single `hash(::UInt, h)` and makes dict-growth rehashing free, regardless of how long the product is. The cache is sound because `QTerm` is immutable and `ne` is already canonicalized before construction, so structurally equal keys always carry the same cached value; `isequal` short-circuits on the cached hash before comparing `ops`. A trusted three-argument constructor lets `_copy_key` carry the known hash across a verbatim copy without recomputing it.

**What a dict entry represents.** A single entry `QTerm(ops, ne) => c` in `arguments` represents the term `c · ops[1] · ops[2] · …` valid for any index assignment satisfying the pairwise constraints in `ne` (each `(α, β) ∈ ne` means `α ≠ β`). The `indices` field on `QAdd` carries the outer summation scope: a `QAdd` with `indices = [i, j]` represents ``\sum_i \sum_j \sum_\text{terms} c \cdot \text{ops}`` where each individual term may further constrain `(i, j)` per its own `ne`. This per-term scoping is what lets a single `QAdd` represent expressions like ``\sum_i a_i a_i + \sum_{i \neq j} a_i a_j`` as two dict entries with different `ne` rather than two separate `QAdd` summations.

**Dictionary keys are full constrained terms.** Each `QTerm` is an ordered operator product plus the pairwise inequality constraints that scope that product. The ordering is canonical (established by `_partial_sort!` and the streaming passes described below), so structurally equal products always have the same `ops`. Like-term collection happens only when both `ops` and `ne` match exactly. This is the key invariant that makes constrained sums correct: `a_i` and `a_i` under `i ≠ j` are distinct stored terms, not one merged term with unioned metadata.

**The schema is visible at the type level.** `QTermDict` is a plain `Dict{QTerm, CNum}` alias — there is no wrapper struct. Iterating a `QAdd` yields `Pair{QTerm, CNum}`, and callers reach `term.ops` / `term.ne` on the key directly. This keeps the storage shape (`(ops, ne) => coeff`, not `ops => (coeff, ne)`) honest at every callsite.

**Zero terms are never stored.** `_addto!` (the single insertion helper used by every code path that adds to a `QTermDict`) builds the key, looks it up, sums the coefficient if it already exists, and deletes the entry whenever the result is zero. The dictionary therefore stays clean without an explicit cleanup pass.

**Summation metadata.** `indices` remains the sum-level metadata on `QAdd`: a `QAdd` with `indices = [i]` represents ``\sum_i (\text{terms})``. Pairwise inequality constraints like `(i, j)` meaning ``i \neq j`` live on the individual `QTerm`s, not on `QAdd` globally. This is why display and round-trip logic can represent mixed-scope sums truthfully.

**Dead-NE invariant.** Every `QTerm.ne` pair must reference at least one index that the term can observe: either some operator in `term.ops` carries that index, the coefficient `c::CNum` depends on it (e.g. an `IndexedVariable` factor), or it appears in the enclosing `QAdd.indices` (sum scope). Pairs failing all three tests encode nothing observable, would only obstruct dict dedup, and so must not survive storage. The `QAdd` inner constructor enforces this by calling `_prune_dead_ne(arguments, indices)`: it walks every term, drops the dead pairs through `_addto_key!` (which performs the like-term merge if dead-stripped keys collide), and is idempotent on already-clean input (no rebuild, no allocation). Algebraic code that constructs `QAdd`s does not need to clean up NE manually; the constructor does it. The predicate is `_depends_on_index_term`, the same one used everywhere else NE/scope dependence is queried, so the rule stays consistent with `_canonicalize!` and `_emit_scaled_by_scope!`.


## Additive reductions and the MutableArithmetics interface

Folding many `QAdd`s with `Base.:+` is O(n²): every `+(QAdd, …)` first calls `_copy_args` to clone the whole backing dict before inserting, so `t1 + t2 + … + tn`, `sum`, and `reduce(+, …)` copy the growing accumulator on every step. `Base.:+` keeps that copy on purpose (value semantics: both operands are user-held, so a single op must not mutate them).

`mutable_arithmetics.jl` recovers the linear cost for reduction drivers, which own a transient accumulator. `_QAddBuilder` wraps one `QTermDict` plus an index vector and threads them through the whole fold; `_accumulate!` reuses the exact in-place ops of the matching `+` method (`_addto_key!` with `_copy_key`, plus an in-place `_merge_indices!` index union), and `_build` materializes the result once. Deferring the single `_prune_dead_ne` (constructor invariant) and `_drop_unused_indices` to `_build` is equivalent to running them per `+` step, because a mid-chain index drop only trims the index vector, it never deletes terms. The result is byte-identical to the `+`-chain; the existing `sum(terms) == manual` test in `algebra_test.jl` is the guard.

The ergonomic entry points are `Base.sum`/`reduce(+, …)` overridden only for `AbstractArray{<:QAdd}` (the element is a `QAdd`, so the result is, which is why `AbstractArray{<:QSym}` is deliberately *not* overridden: `acts_on`/`prefactor` map operators to non-`QField`s). The fast path is the bracketed comprehension `sum([ω[k]*a[k]'*a[k] for k in 1:M])`: an array comprehension yields a `Vector{QAdd}` the override accelerates, whereas a bare generator `sum(… for …)` cannot be intercepted by element type and stays on Base's generic `+`-fold. The MA interface (`mutability`/`operate!`/`add_mul`/`promote_operation`) is implemented on `_QAddBuilder` so `MA.@rewrite` and manual `operate!!` loops work too; we do not subtype `MA.AbstractMutable` (it would drag the whole `QField` hierarchy, including the immutable `Op` leaf, under MA and risk `+`/`*` ambiguities).

This is orthogonal to the analytical sum `Σ`/`∑`, which sums symbolically over an index into a single indexed `QAdd` and already accumulates in place internally. `Σ` is untouched; the only integration requirement is that the accumulator merge `.indices`/`ne` exactly as `+(QAdd, QAdd)` does, which holds because it reuses the same helpers.

**Why the product path is left as repeated `*`.** A single `*(QAdd, QAdd)` already accumulates all term-pairs into one shared dict, and `^(QSym, n)` builds the whole operator string and canonicalizes once, so both are already optimal. Chained products, `^(QAdd, n)`, and `reduce(*, …)` *do* re-materialize an intermediate `QAdd` per step, but a prototype double-buffered product builder (deferring the per-step wrapper/prune/absorb) measured ~1.0× (no meaningful time or allocation win, occasionally slower from the upfront copy). The reason is structural: a product's dominant cost is the intrinsic distributive expansion (`_emit_product!` over all `|a|×|b|` term pairs plus canonicalization), which any builder pays identically; only the cheap per-step wrapper is recoverable, and fusing it across factors would complicate the `_absorb_pinned_sums` index-scope logic for three-plus indexed sums. So the MA work is additive-only by design.

## Operator sorting

```julia
@enum SiteCmp::UInt8 Less Equal Undetermined Greater

_site_compare(a::Op, b::Op, ne::Vector{NonEqualPair})::SiteCmp
_partial_sort!(ops::Vector{Op}, ne::Vector{NonEqualPair})
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
| `_commute_ops` | applies swaps for adjacent same-site pairs whose `_can_commute` is false; emits the swap plus up to two residual branches |
| `_expand_gs_ops` | applies `σᵍᵍ → 1 - Σ_{k≠g} σᵏᵏ` to every ground-state projector; opt-in, not part of the default pipeline |
| `_substitute_ops` | walks `ops` applying a substitution dict; forks when a value is a multi-term `QAdd` |

The terminal sink, `_canonicalize_to_dict!`, is the only place a `QTerm` is constructed during a pipeline run. It builds `QTerm(ops, ne)`, looks up in `out`, sums the coefficient, drops zero entries.

The aliasing rules are what justify the streaming style: there is no intermediate term-list materialization. Single-output passes mutate `ops` in place and the sink receives the same `Vector`; forking passes call `copy(ops)` before mutating each branch except possibly one, so no two branches share a mutable operator string; the terminal sink takes ownership of whatever it receives. Every pass documents which rule it follows.

The pipeline establishes the **canonical-form invariant**: `ops` is in partial-canonical order, no adjacent provably-same-site pair has a remaining commutator residual or reduction, and like terms are collected. The invariant deliberately says nothing about `σᵍᵍ` — see the next section.


## Simplification vs. normal ordering

The package exposes four entry points that combine the pipeline above in different ways.

`normal_order(expr)` re-streams each term through `_stream!`. Eager `*` already produces canonical form, so on the output of `*` it is idempotent; it earns its keep as a finalizer for hand-constructed expressions and as the second half of `simplify`. `simplify(expr)` runs `normal_order` and then walks the resulting terms once more, applying `Symbolics.simplify` to each coefficient and dropping summation indices that no surviving term depends on. The expensive per-coefficient simplification deliberately lives in this outer pass rather than inside the streaming pipeline: it runs once per surviving term, not on every dict insertion. `expand(expr)` distributes symbolic prefactors only — `(a + b)² → a² + 2ab + b²` — and leaves the operator structure untouched. `expand_completeness(expr)` applies the ground-state identity, described below.

`commutator(a, b)` is `a*b - b*a` in the general case, with a fast path on `QSym × QSym`: when exactly one direction of the pair is non-canonical, the commutator equals the residuals returned by `_commute_pair`, so the call short-circuits without running the full pipeline twice and a subtraction.

**Why `σᵍᵍ` stays atomic in canonical form.** Ground-state projectors are legitimate atoms in the canonical basis. The completeness identity `σᵍᵍ = 1 - Σ_{k≠g} σᵏᵏ` is mathematically exact, but applying it eagerly multiplies every product containing a `σᵍᵍ` by `n_levels`. A product with `k` ground-state factors balloons to `n_levels^k` terms before the surrounding context has a chance to cancel anything, and like-term collection across operations cannot recover the original compactness once the explosion has happened. Keeping `σᵍᵍ` atomic also serves downstream consumers: in mean-field expansions, `⟨σᵍᵍ⟩` is a single moment that solvers carry through their equations directly, while `1 - Σ ⟨σᵏᵏ⟩` is a sum of `n_levels - 1` moments tied together by an identity constraint that meanfield code would otherwise have to recognize and dedupe.

`expand_completeness(expr)` is the explicit handle for the cases where the expansion is genuinely wanted — converting to a basis in which `σᵍᵍ` is a dependent quantity, or feeding into code that expects the identity already materialized. Internally it walks each term through `_expand_gs_ops` (which forks `n_levels` ways per `σᵍᵍ`) and re-streams the output. User-constructed `σᵍᵍ` is the matching boundary case: `Transition(h, :σ, g, g)` returns a plain `Transition`, and canonicalization only happens when the operator participates in a `*`, so direct construction and inspection stay cheap.


## Unitary transformations

`UnitaryTransform` stores a unitary as its action on the generators of each site it touches, never as an operator expression: `exp(α a† - conj(α) a)` is not a `QAdd`, but `a ↦ a + α` is. `substitute` with `Op` keys is the engine, so the type exists only for what raw substitution cannot do: guarantee the rule set is complete, carry the inverse map, carry the gauge.

**Traps.** Each fails silently, or in a way a single round-trip test cannot catch.

* **Rule sets must be complete.** Eager canonicalization stores `σ¹²·σ²¹` as the leaf `σ¹¹`, so rules covering only the off-diagonals give a wrong answer with no diagnostic. Completeness is a constructor invariant; application rejects any operator on a covered site that has no rule.
* **Never build a conjugate coefficient by negating a transcendental's argument.** `exp(±im·ϕ)` Euler-splits into four distinct atoms, since `cos(-ϕ)` is not folded to `cos(ϕ)`, and the squeeze residuals never close. Negate the imaginary part of the real/imaginary decomposition instead. Same trap for `exp(-r)`; use `1/exp(r)`, the same atom at exponent `-1`.
* **The N-level contraction runs over `W`'s row index**: `U†|i⟩⟨j|U = Σ_kl conj(W[i,k])·W[j,l]·|k⟩⟨l|`, `W[k,l] = ⟨k|U|l⟩`. The transpose disagrees with a matrix oracle on every rule, so this direction is pinned by test.
* **The two rule sets compose in opposite orders**, since `conjugate(A, U₁U₂) = U₂†U₁†AU₁U₂` but `conjugate(A, inv(U₁U₂)) = U₁(U₂AU₂†)U₁†`. Writing both the same way breaks `inv` only for non-commuting factors. The gauge composes as `conjugate(gauge(U₁), U₂) + gauge(U₂)`.
* **Coefficient reduction iterates to a fixpoint over every pair present, and runs ungated when zero-testing.** Keeping a rewrite only when it shortens the term list is right for display but non-monotone: a blocked rewrite reads as a nonzero residual and `is_canonical` false-negatives.
* **A hand-built `Monomial` must keep `syms` sorted by `_fkey`**, which the merge routines assume; an unsorted one silently fails to cancel.
* **Rules key on the exact `Op`, and site identity includes the name.** Two differently named Fock modes on one `FockSpace` are distinct sites, so only one may need rules; an indexed operator never matches a bare rule key and throws. The same name-blindness caps a `PhaseSpace` at one canonical pair, since a second differently-named `Position` there is the first one relabelled, and the coverage error points at `⊗` instead.
* **Each constructor spells its inverse out, and that redundancy is deliberate.** What differs between them is the sign convention, which is the thing worth reading at the call site: the two rotations are antisymmetric in the off-diagonal, the two-mode squeeze symmetric, the quadrature squeeze a diagonal swap rather than a negation. Factoring the shape out would take that as a parameter and hide it, and would cover three of the six constructors. The second spelling is checked, not merely repeated: `canonicality_residuals` round-trips every generator through both directions, so a disagreeing pair fails `is_canonical`.
* **`transform` collides with `QuantumOpticsBase.transform`**, so a script loading both must qualify `SecondQuantizedAlgebra.transform`. Adding methods there instead would put a purely symbolic operation behind a weak dependency, and the two share no contract for dispatch to use.

**Coefficient reduction does not go through the CAS.** Rotations produce `cos² + sin²` and squeezes `cosh² - sinh²`, which fold on the parameter-polynomial tier in a few µs by rewriting `hi^e → (1 + sign·lo²)^(e÷2)·hi^(e mod 2)`; `Symbolics.simplify` needs around 100 µs and above degree 2 with a composite argument does not fold at all. Since `_is_atom` requires a one-argument call's *argument* to be an atom, `cos(θ)` is already a polynomial factor while `cos(ω·t)` sits on the `Complex{Num}` tail and is swapped for transient stand-in symbols, reduced, and swapped back. Widening `_is_atom` instead would move those coefficients onto the polynomial tier, where `_simplify_prefactor` early-returns and everything the CAS folds beyond Pythagoras is silently lost from the public `simplify`.

The transient swap is the expensive path, and stays so: `conjugate` through a rotation with a composite angle costs around 240 µs against 2 µs for a `RotatingFrame`, because the two `substitute` walks in and out dominate. The stand-in symbols themselves are cached at `__init__` rather than re-derived, which took a third off the reduction; caching them at load instead would be wrong, since a symbol minted during precompilation is absent from the runtime intern table and a later rebuild would mint a second one under a different `_fkey`.

**A phase is a single atom, and that is a different case.** `expim(x)` is atomic whatever its argument, because the objection above does not apply: the CAS has never heard of the head, so there is nothing to lose. What it buys is phasor arithmetic on the polynomial tier: multiplication adds exponents, conjugation negates them (a phase is unimodular), and `expim(x)·conj(expim(x))` collapses to `1` with no relation declared. Written as `cos` and `sin` it is two unrelated factors, and every identity between them has to come from outside:

| identity | on the trig tail | as phase exponents |
|---|---|---|
| `cos² + sin² = 1` | a `ParamRelation`, iterated to a fixpoint | exponents `+1, -1` sum to `0` |
| `cos(A+B) = cos A cos B - sin A sin B` | no relation states it | frequencies were never separated |
| `exp(im·x)·exp(-im·x) = 1` | four distinct atoms, never closes | one atom, exponents cancel |

**Adding a coefficient atom is a checklist, and `test/expressions/atom_conformance_test.jl` is that checklist as a runnable test.** Add an entry to its `ATOMS` tuple and make it pass. The contract is spread over four layers, and every miss is silent rather than an error:

| layer | what the atom owes |
|---|---|
| head | `promote_symtype` and `promote_shape`. These take part in hash-consing, so without them every `maketerm` rebuild (each `substitute`, each `Postwalk`) returns a *different* atom whose symtype has decayed to `Real`, and `conj` on that is the identity. |
| construction | a canonicalising producer (argument normalisation, sign orientation, literal folding) that every call site uses, public and internal. Never `_atom_coeff` directly. |
| type | the public constructor returns a `Coeff`, never a `Num`. `Base.conj(::Num) = x` is unconditional and unoverridable, and `Complex * Num` splits into real/imag halves, so a `Num`-typed complex atom cannot be made correct. |
| `cnum.jl` | `_is_atom`, the `_recognize` arm, the `_recognize` `conj` arm, `_conj_poly`, the negative-exponent lowering in `_term_to_num`, `_unary_at_zero` if it is reachable at all. |
| adjoints | `qadjoint` and `inner_adjoint` rebuild through `maketerm` with *conjugated arguments*, which is wrong for any atom whose conjugate is not `f(conj(x))`. |
| bounds | `_sup_norm`, and `_phase_poly_bound` if the atom is unimodular. An unbounded atom makes every floating-point residual read as nonzero, so `is_canonical` false-negatives. |
| numeric | `_fold_const` and `_numeric_average`. Both throw on an unknown head. |
| display | `SymbolicUtils.show_call` *and* `Symbolics._toexpr_op`. They are independent renderers; registering one leaves the other printing the raw head name. Both reach the coefficient through one lowering entry point in `printing/printing.jl`, so an atom that needs folding on the way to the page (a conjugate phase pair back to `cos`/`sin`) is handled there once rather than per renderer. |

Three rules make the atom canonical. The argument is `expand`ed before interning, so `(ω + 2J)·t` and `ω·t + 2J·t` are one atom. The sign is oriented, so `expim(-x)` interns at exponent `-1`, without which `ω·a†a - ω·b†b` carries unrelated atoms and never cancels; `_recognize` re-orients any `expim` it meets. And an integral scale must stay integral, because **atom identity is `===`**: `(-g)·(-1.0)` and `g` are `isequal` and hash alike without being the same object, `_fkey` is `objectid`, so a `Float64` scale of magnitude one mints a second `expim(g·t)` that no merge ever combines. Nothing is *wrong* afterwards, which is what makes it worth a rule: `to_num` prints the sum folded, the residual is still bounded, and the only symptom is that `Σ|scalar|` counts `0.4999999999999999` and `-0.5` as two terms of size `½` rather than one of size `10⁻¹⁶`, so a frame that is exactly canonical reads as `false`. Unimodularity also bounds the residual for free, so `_bounded_by` never lowers to SymbolicUtils, and display folds a conjugate pair back to `cos`/`sin` so a rotation still reads as one.

The line is unimodularity, not staticness: every factor of modulus one is a phase atom, `Rotation(a, θ)` and the `exp(imϕ)` of a squeeze included, which is why those certify exactly at `atol = 0` with no relation declared. What keeps trig are the factors that are not phases. A real orthogonal mixing of two Hermitian generators (`cos θ·x + sin θ·p`, a spin rotation about an axis) has genuinely real coefficients, and rewriting them as `(p + p⁻¹)/2` would trade a readable answer for exactness it does not need, since `cos(θ)` on a bare symbol is already a polynomial atom that folds. `cosh`/`sinh` are not phases at all and keep their `ParamRelation`.

**The gauge term.** For a time-dependent frame the transformed generator of motion is `U†HU + im·(∂ₜU†)U`; `transform` returns it and `conjugate` omits it, since the struct cannot know whether its argument is a Hamiltonian or an observable. For `U = exp(-im·θ·G)` the gauge is `-∂ₜθ·G`, but only when `[G, ∂ₜG] = 0`:

| Constructor | Time-dependent gauge |
|---|---|
| `Rotation` (all signatures except `(σ, W)`) | `-∂ₜθ·G` |
| `Squeeze(x, p, r, t)`, `Squeeze(a, b, r, t)`, `Squeeze(a, r, ϕ, t)` with `t`-free `ϕ` | `-∂ₜθ·G` |
| `Squeeze(a, r, ϕ, t)` with `t`-dependent `ϕ` | su(1,1) closed form below |
| `Displace(a, α, t)`, `Displace(x, p, dx, dp, t)` | closed form below |
| `Displace(a, A::QAdd, t)` | the same closed form, with `[A, Ȧ]/2` an operator rather than a phase |
| `Rotation(σ, W, t)` | `im·Ẇ†W`, entrywise |
| `UnitaryTransform(G, θ, t)` | `-∂ₜθ·G`, exact since `G` is time-independent |
| `RotatingFrame(H0, t)` | `-H0`, the `θ = t` case of the above |
| `Bogoliubov` | static only; there is no `(G, θ)` pair to differentiate |
| any of the above on an indexed family | the single-site gauge, summed over the wildcard |

The displacement generator moves, with `[A, Ȧ]` a nonzero c-number: from `e^{-A}∂ₜe^{A} = Ȧ - [A,Ȧ]/2` the Fock gauge is `-im(α̇a† - conj(α̇)a) - im(conj(α)α̇ - α·conj(α̇))/2`. That c-number contributes only a global phase, so no observable test sees it and both displacement gauges are pinned against a central-difference matrix oracle in a truncated Fock basis.

For a moving squeeze phase `[A, Ȧ]` is an operator, so no `-∂ₜθ·G` form exists and BCH does not terminate. Imposing `∂ₜ(U†gU) = -im[K, U†gU]` for `K = im(∂ₜU†)U` on the Bogoliubov pair `(cosh r, e^{iϕ}sinh r)` determines `K` up to a c-number, which su(1,1) fixes too since it closes on `a†a + 1/2` rather than `a†a`:

```
K = ∂ₜϕ·sinh²r·(a†a + 1/2) + β·a†²/2 + conj(β)·a²/2,   β = e^{iϕ}(∂ₜϕ·sinh r·cosh r - im·∂ₜr)
```

collapsing to `-∂ₜr·G` at constant `ϕ`. None of it is a global phase, so the oracle pins the whole thing, at a deeper truncation since a squeeze spreads far up the ladder. `Rotation(σ, W, t)` needs no scalar `θ` at all: the gauge is the operator of `im·Ẇ†W`, Hermitian for free since `∂ₜ(W†W) = 0`, reducing to `-Σ_k E_k σ^{kk}` for a level-diagonal `W`.

**`is_canonical` needs no representation of `U`.** A linear map on the generators of a site is induced by a unitary exactly when it preserves the defining algebra and Hermiticity, and that is sufficient rather than merely necessary: CCR plus Hermiticity makes a Fock map Bogoliubov, unitarily implementable by Stone-von Neumann, and preserving the matrix-unit algebra makes an N-level map inner by Skolem-Noether. Since the algebra constrains only the forward map, the residuals also include the round trip `conjugate(conjugate(g, U), inv(U)) - g`; without it a transform whose `inv_rules` is not the inverse passes every algebraic test while `inv` returns garbage. They run with ungated reduction and carry an `atol`, since floating-point coefficients cancel only to rounding.

Hermiticity is load-bearing on both sides: the phase-space CCR admits `x -> x + im`, which no unitary induces, so the pair carries its own Hermiticity residual beyond `[p̃, x̃] = -im`. N-level products are tested against the first row and column only, since `σ^{ij} = σ^{i1}σ^{1j}` makes the rest collapse by associativity (651 residuals down to 76 at `n = 5`).

Preserving each site's algebra separately is not the whole condition, because operators on *different* sites commute and a mode-mixing map can break that while passing every same-site test: flip the sign of one cross term of a two-mode squeeze and both `[ã, ã†]` and `[b̃, b̃†]` stay right while `[ã, b̃]` becomes `-2·cosh r·sinh r`. So a quadratic sweep over generator pairs on distinct sites runs too, gated on whether any rule actually sends a generator off its own site: a map that does not cannot break it, and that gate is what keeps the sweep off every single-site constructor.

**`RotatingFrame` derives its map from `ad_{H0}`, not from term shapes.** For `U = exp(-im·H0·t)` each generator obeys `dg̃/dt = im·[H0, g̃]`, so where the commutator closes on the generators of the sites `H0` touches it is a matrix `A` and the rule map is `exp(im·A·t)`. Classifying terms syntactically is what used to leave a beamsplitter, `(x² + p²)/2` and a parametric `a² + a'²` unreachable while `a'a` worked; they are all the same construction. The gauge stops being a table too, since `im·(∂ₜU')U = -H0` for any time-independent `H0`.

The exponential factors over the connected components of `A`. Size one is a phase. Size two is written `w·B` with `B` numeric, so the rate is read off rather than square-rooted out of the determinant, and `B² = I` versus `B² = -I` is the difference between a rotation and a squeeze. Adjoint-related components are derived with `adjoint` rather than exponentiated twice, since two independent runs can pick opposite rate signs and break Hermiticity.

**Any number of generators may couple**, by the same move that makes an N-level frame work: split the component as `A = d·I + w·C` with a real rate `w` and a *numeric* Hermitian `C`, diagonalize `C` through `eigen`, and give each eigenvalue one phase atom. `M·M† = I` then reduces to the numeric orthonormality of the eigenvectors, so nothing has to know that the eigenvalues are related, which is exactly what the trig spelling could not survive past two. The rate stays symbolic, so a hopping chain or ring of identical modes qualifies with symbolic `ω` and `g`, as does a spin rotation about an arbitrary axis. Three details earn their keep:

* **The shared diagonal comes out first.** `d` is what every generator of the component has in common, and on its own phase it cancels out of a number-conserving observable instead of riding inside one atom per eigenvalue. Pulling it out is also what admits `ω·Σa†a + g·(hopping)`, whose block is not traceless and so has no `2x2` closed form either.
* **The spectrum is snapped onto shared magnitudes, and onto integers.** `eigen` returns `±√2` as two floats that are not exact negatives, so one atom per raw eigenvalue would make `p` and `p⁻¹` unrelated atoms instead of one at opposite exponents, and a degeneracy would split in two. Integral magnitudes are then rounded exactly, so the rate can be scaled by an `Int` and satisfy the atom-identity rule above; a circulant ring returns its degenerate `-1` as `1.0000000000000002`, and its spectrum `{2, -1, -1}` collapses to the two atoms `expim(g·t)` and `expim(2g·t)`.
* **The exact `2x2` form is still tried first.** It reads the rate straight off the block instead of through numbers, so a beamsplitter keeps an exact `cos(g·t)` rather than a `0.4999999999999999` weight, and it is the only path that covers a hyperbolic block, whose `C` is not Hermitian and has no unitary exponential to diagonalize.

The eigenvectors are numbers, so these frames cancel to rounding rather than exactly (`atol`), and `transform(H0, U)` leaves rounding *proportional to the parameters of `H0`*, which no sup-norm bound can wave through. The atoms themselves stay exact.

`DressedFrame(H0)` is the static counterpart, diagonalizing numeric energies and couplings through `eigen` at any number of levels (hence `LinearAlgebra` as a hard dependency, free since Symbolics already loads it). Symbolically the level matrix splits into the blocks `H0` actually couples: one level is already diagonal, two is the rotation by its mixing angle `atan(2Ω, Δ)/2`, and three or more works when `_common_factor` splits the block as `d·I + w·C` with a numeric Hermitian `C`, since a dressed basis is *only* the eigenvector matrix and `C`'s eigenvectors are numbers while the energies `d + w·μ` stay symbolic. A uniformly coupled chain or ring qualifies with both its spacing and its coupling free; two independent symbolic couplings do not, and are refused naming the levels. The two-level path deliberately runs first: routing it through the split would trade its exact symbolic mixing angle for numeric eigenvectors whenever the coupling ratio happened to fold to a literal. What is left refused is a *general* symbolic eigendecomposition, and the usual citation for that is wrong below `n = 5`: cubics and quartics are solvable by radicals, and what actually blocks `n = 3` here is casus irreducibilis, since three real eigenvalues of a real symmetric matrix cannot be written in real radicals and the complex cube roots that result would poison every residual. Abel-Ruffini only applies from `n = 5`. Structured matrices are a different matter and `eigen` is a fallback rather than a floor here. Two structures are detected and answered exactly. A uniform nearest-neighbour chain is symmetric tridiagonal Toeplitz and diagonalizes by the discrete sine transform (`V[j,k] = sin(jkπ/(n+1))`, `E[k] = a + 2b·cos(kπ/(n+1))`, Chebyshev polynomials of the second kind); any circulant, not only the periodic band, diagonalizes by the DFT (`V[j,k] = exp(2πi·jk/n)`, `E[k] = Σⱼ cⱼ·exp(2πi·jk/n)` over the first row). In both the eigenvector matrix depends on the size alone, so it is exact whatever the entries are, and only the eigenvalues read them. Both closed forms are classical; SymbolicDiagonalization.jl catalogues these and other structure classes, and is GPL v3, so nothing may be copied from it into this MIT package and it cannot be a dependency either.

**The exact values are a rational multiple of one squarefree root, and that normal form is the coverage limit.** `Symbolics` combines neither `√(1//2)·√2` nor `√2·√2`, so leaving the arithmetic to it builds entries that are correct and unreadable; `_RootVal` keeps `q√r` in lowest terms instead and lowers to a `Num` only at the end. The table of `sinpi` values that fit that form covers denominators 1, 2, 3, 4 and 6, which is chains of 2, 3 and 5 levels and rings of 3, 4, 6, 8 and 12. `sinpi(1//12)` and `sinpi(1//10)` are sums of two roots and `sinpi(1//8)` is a nested one; a size needing them declines and falls back to `eigen` rather than answering with floats, which is the property that matters. Widening the form to a *sum* of `q√r` terms would add denominators 5, 10 and 12; nested radicands cannot be reached at all.

**Exactness and phase cancellation are mutually exclusive for a radical coefficient, so a frame keeps the numeric spectrum.** A phase cancels against its conjugate on the polynomial tier, where it is an atom carrying an exponent, while an exact rational survives only on the `Complex{Num}` tier, since `Monomial`'s scalar is a `ComplexF64`. A radical entry forces the coefficient onto the symbolic tier and the residual then stalls at `-1/4 + (1/4)·exp(ix)·exp(-ix)`; admitting the radical as a polynomial atom instead restores the cancellation but lowers `g√2` to a float and respells `expim` as `cos`/`sin`. Nothing bridges that without giving the monomial scalar an exact type. So the exact spectrum is taken only where no phase has to cancel: a static dressed basis keeps it, and a frame keeps it only when its weights are rational. The two-level frame is exact, but its off-diagonal is not reduced to zero: that needs `cos(θ)²Ω - sin(θ)²Ω - cos(θ)sin(θ)Δ = 0`, a relation between the mixing angle and its own arguments that a `ParamRelation` cannot state. Its dressed *energies* need a `sqrt` that no relation closes either, which is why a symbolic two-level block is a `DressedFrame` and not a `RotatingFrame`.

**A single N-level site is exactly solvable whatever the shape of `H0`**, since `ad_{H0}` always closes on the `n²` transitions. Off the level diagonal that is one component of size `n²`, which the general path above would take, but `RotatingFrame` factors it instead: `exp(-im·H0·t) = V·exp(-im·Λ·t)·V†` is `DressedFrame(H0)` wrapped around the diagonal frame of the dressed energies. The factored form is not a fallback but the better answer, since it phases the `n` levels where the general path would phase up to `n²` eigenvalues of `ad`.

**An N-level frame carries one phase per level, not one per transition**, which is the whole reason three or more levels work. `σ^{ij}` rotates at `E_i - E_j`, so phasing per generator mints an atom for every difference and `E₁-E₃ = (E₁-E₂) + (E₂-E₃)` becomes invisible, leaving residuals that need angle-addition. Carrying `P_i·conj(P_j)` on `n` level atoms makes every such relation exponent arithmetic, and a symbolic level-diagonal frame comes out exactly zero with no tolerance. Generally: when the generators are a tensor of two indices, phase the indices rather than the pairs.

**The scalar the exponentiation runs against is arbitrary.** `_exponentiate(G, θ)` is the whole driver, and `θ` is a pure multiplicand everywhere in it: nothing differentiates it, so a static `exp(-im·θ·G)` and a frame `exp(-im·H0·t)` differ only in the gauge their callers install. `RotatingFrame(H0, t)` is literally the `θ = t` case, since `_gauge(H0, t, t) = -∂ₜt·H0 = -H0`. What the frame must *not* be given is a moving `H0`: `exp(-im·H0·t)` is not its propagator and `-H0` is not its gauge, and before the guard `RotatingFrame(ω·cos(t)·a†a, t)` returned `exp(-im·cos(t)·t·ω)·a` where the phase is `exp(-im·ω·sin(t))`, built the wrong gauge, and certified. The commuting case `H0 = f(t)·G` is exactly `UnitaryTransform(G, θ, t)` with `θ = ∫f`; there is no symbolic integrator, so the user supplies it.

The hand-written constructors are deliberately *not* reimplemented on top of it. Their spellings are load-bearing (`_quad_squeeze`'s `1/exp(r)` is one atom at exponent `-1` and cancels exactly where the generic path would write `cosh ± sinh`), and the equivalence is a test rather than a shared implementation. That test has to compare numerically: the generic path spells a rotation on phase atoms and the hand-written one on `cos`/`sin`, and closing `(p + p⁻¹)/2 = cos θ` needs Euler's formula, which no `ParamRelation` states.

**A drive makes `ad` affine, and the shift absorbs it.** `[η(a + a†), a] = -η` is a c-number, so `[G, g_k] = Σ_j A[j,k] g_j + c_k`, and `g_k -> g_k + s_k` is homogeneous exactly when `Aᵀs = c`. The image is then the homogeneous one plus `Σ_j M[j,k] s_j - s_k`, with `M` read straight back out of the rules already built rather than exponentiated again, which is what keeps the forward and inverse shifts in step with their own directions. A component with `A = 0` has no shift and integrates to `im·c·θ` instead, the displacement a pure drive produces. A component where `c` is outside the range of `Aᵀ` is resonant, grows without bound, and throws: a dark mode of a hopping chain driven on resonance is the reachable instance. Any particular solution of `Aᵀs = c` will do, since two differ by a kernel element and `M` fixes those.

The affine part is the one piece of a frame **`is_canonical` cannot pin**. A constant commutes out of the CCR residual, and the round trip only catches a forward and inverse that disagree, so a shift wrong by the same factor in both directions passes every residual. Its magnitude is pinned by the `to_numeric` oracle instead.

**An involution exponentiates in closed form at any size.** `C² = κ·I` makes the series close on `I` and `C`, so `exp(im·w·θ·C) = cos(x)·I + im·(sin(x)/√κ)·C` with the rate still symbolic, where `_frame_multi_rule!` would need a Hermitian `C` and bake numeric eigenvectors into the rules. That is what admits a detuned parametric amplifier, whose block has an unshared diagonal (no `2x2` form) and is not Hermitian (no `eigen`). The exact `2x2` path still runs first: its sign orientation is what the printed two-mode forms are pinned to.

**`κ` is never rooted numerically, and that is what makes the whole block exact.** No rate is divided out of `C` first, so `κ` absorbs any scale and is whatever the coefficients make it, symbolic included; deciding `C² = κ·I` is then polynomial cancellation rather than a tolerance. The root is left as an unevaluated `ssqrt` node, which is the load-bearing choice. `sqrt` squares back to `κ` algebraically but raises a `DomainError` the moment a substitution makes `κ` negative; `ssqrt` substitutes on either branch but SymbolicUtils will not square it back, so the one identity `ssqrt(x)² = x` is supplied where residuals are tested. Neither spelling has both halves on its own. Building the node with `term` rather than by calling `ssqrt` avoids two traps: `ssqrt` returns `nothing` on a `Num`-wrapped constant, and on a bare positive one it folds to a `Float64` and puts the rounding straight back.

The payoff is that one `cos`/`sin` form covers every regime of a detuned parametric amplifier, including the ones no float could certify. Both `cos(√κ·θ)` and `sin(√κ·θ)/√κ` are *even* in `√κ`, hence single-valued entire functions of `κ`, so the branch cannot affect the value and `cos² + sin² = 1` closes on either side. Above threshold, at threshold and below it all certify with no tolerance at all, and so does a fully symbolic `Δ·a'a + g(a² + a'²)` whose sign is not even decidable. A real root is still preferred when one exists exactly, since `√(g²)` reduces to `g` by halving a `Rational` exponent on the polynomial tier and a real rate belongs on a single phase atom rather than split across `cos` and `sin`. An exact *rational* `κ` is deliberately left unrooted: a rational rate makes only part of the pair fold back, so the two images end up spelled differently and stop cancelling.

**`κ = 0` is reachable, and it is the threshold.** A parametric amplifier driven exactly at `Δ = 2g` has a nilpotent block, where the series stops after one term and the frame is a polynomial in `θ` carrying no atom whatsoever.

**Exactness is inherited from the input, so `3//10` and `0.3` are not the same request.** A coefficient stays exact when it does not round-trip through `ComplexF64`, so `3//10` and `1//3` are held as rationals while `0.3` is a native float that is not three tenths. Integers and dyadic rationals (`2`, `1//2`, `1//4`) are native *and* exact, which is why choosing `0.5` over `1//2` buys nothing. Scaling a Hamiltonian to integers is better still: it can make `κ` a perfect square and collapse the frame onto one phase atom. Float input is not rejected, but its residuals then cancel only to rounding and need the `atol` path, which is why that path applies the Pythagorean identities before bounding: two separately bounded terms that happen to cancel bound only by their sum, where the collapsed form bounds by the rounding itself.

**An operator-valued displacement amplitude gates its own site rather than covering it.** `U'gU = g + A` truncates at first order for any `A` commuting with the site it displaces, which is the polaron class, and `A` must also be normal or `[A, A']` survives in the commutator of the images. But `A`'s *own* site is then neither covered nor free: `U'XU = X` for every `X` commuting with the exponent, and for every `X` that does not the image is `X` times a displacement operator, which is no polynomial in the generators. Passing those through untouched is wrong with no diagnostic (a matrix oracle puts `U'σ¹²U` a distance `4.9` from `σ¹²`), and `is_canonical` cannot see it, since it only tests covered sites. So the amplitude is recorded and `conjugate` gates every term against it. Per *sub-product*, not per operator: `a` and `a†` each fail against an `a†a` amplitude while `a†a` passes, and `a†a` is the whole point. The check is memoized per call and skipped entirely when no amplitude is operator-valued.

**A collective site is the matrix-unit contraction without the matrix-unit algebra.** `U'S^{ij}U = Σ_kl conj(W[i,k])·W[j,l]·S^{kl}` holds verbatim for `S^{ij} = Σ_a |i⟩_a⟨j|` under `⊗_a U`, because conjugation is linear in the atom sum; it never uses the product law `S^{ij}S^{kl} = δ_{jk}S^{il}`, which collective operators do not obey. Only the residuals change. Hermiticity and the bracket law `[S^{ij}, S^{kl}] = δ_{jk}S^{il} - δ_{li}S^{kj}` carry over, the latter probed against a pivot row and column only (`E_ij = [E_ip, E_pj] + δ_ij E_pp` plus Jacobi extends it to every pair, `2|L|-1` probes instead of `|L|²`). What replaces the completeness identity is the invariance of the atom number `Σ_i S^{ii}`, which is *not* the identity here, and that residual is load-bearing rather than a consolation prize: the contragredient map `S^{ij} -> -S^{ji}` satisfies the bracket exactly, is Hermitian-consistent exactly, and is involutive so it round-trips exactly, yet no unitary induces it, and number invariance is the only residual that rejects it. The level count is the other difficulty: it lives on the space and a collective `Op` does not reference it, so it comes from `size(W)` for a rotation and from `1:maximum` of the levels `H0` names for a frame.

**An indexed family is `⊗ᵢUᵢ`, and site locality is what makes that sound.** An abstract index is the wildcard and a per-slot one (`i(3)`) stays an ordinary site, so both readings exist with no second spelling and no index minted on the user's behalf. Application is instantiate-then-substitute: the targets are collected during the coverage walk, the rules are `change_index`ed onto them, and the existing substitution engine runs unchanged. The relations have to move too, or a family squeeze declares `cosh(r(i))² - sinh(r(i))² = 1` and then false-negatives at every index but `i`. The load-bearing restriction is that a family rule image stay on the site it acts on. Allowing an off-family site is not merely unsupported but unsound: `a_i -> cosh(r)·a_i + sinh(r)·b†` passes every same-site residual while `[ã_i, ã_j†] = -sinh²r`, because the instances then share `b`. With the restriction the instances have disjoint support, so they commute for free and no cross-index residual is needed, which is fortunate: writing one would need a second `Index` the algebra may not mint and an `assume_distinct_index` wrapper to escape `Undetermined`. It also keeps substitution clear of diagonal splitting, since the multiset of `(space, index)` pairs per term is unchanged, so no same-space index pair appears that was not already there and substitution can end in `_canonicalize!` rather than `_accumulate_with_diag!`. Residuals at the wildcard certify every instance because `change_index` is a `*`-algebra isomorphism between the sites that the family map commutes with by construction. A time-dependent family gauge must be `Σ_w`, since `im(∂ₜU')U = Σᵢ im(∂ₜUᵢ')Uᵢ`, and that wrap lives in `_with_gauge`, which every timed constructor routes through and `RotatingFrame` bypasses with its already-summed `-H0`.

Two extensions are deliberately left open. A family whose rules span two spaces (`⊗ᵢ` of a beamsplitter, a two-mode squeeze, a `Bogoliubov`) needs site locality relaxed to index-uniformity plus an application-time guard against creating a same-space index pair inside a `Σ` with no diagonal split. An off-diagonal collective `H0` needs the dressed factorization parameterized on the operator builder, which touches `DressedFrame` too. Both throw today, naming the escape hatch.

**Mixed sites are outside the class** for a reason no implementation effort reaches: for `Δσ²² + g(a†σ¹² + aσ²¹)`, `[H0, σ¹²]` produces `aσ²²` and iterating never closes on a finite set of polynomial generators (the Jaynes-Cummings propagator involves `cos(gt√(a†a+1))`). `_ad_matrix` throws naming that product.

**`_residual_iszero` bounds rather than samples.** Floating-point frames leave rounding on `cos`/`sin` atoms, which are not atoms in the `Poly` sense and never collapse to a bare number, so the check computes a sup-norm bound over the expanded coefficient instead of evaluating at sample times: `cos` and `sin` contribute `1`, sums and products combine by the triangle and product inequalities, anything else makes the bound `nothing`. A residual bounded below `atol` is zero everywhere, so the verdict stays sound and unbounded factors keep exact transforms on the exact path.


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

**`Index` struct (`isbits`):**
```julia
struct Index
    name_id::Int32     # interned display name (intern.jl); 0 == anonymous site / NO_INDEX
    range_id::Int32    # interned range Num (intern.jl); 0 == no range
    space_index::Int32 # which space
    slot::Int32        # 0 == abstract; k>0 == concrete site k
end
```
`Index` is `isbits` so that the `Op` embedding it can be `isbits` (see *Single concrete `Op`*). The two `Num` fields it used to carry (`range`, `sym`) cannot be `isbits`, so they are interned/reconstructed instead:

- **`range`** is interned to a `range_id::Int32` (a module-global `Num` table in `intern.jl`) and recovered with `index_range(idx)::Num`, which returns the user's original `Num` (so a symbolic range `N` stays usable in coefficients). Range lives on the index, not lifted to the `Σ` scope, so a bare `a_i` knows its range without consulting the enclosing sum (term-locality), and the downstream `QuantumCumulants.jl` reads it off operator-attached indices.
- **`sym`** is dropped and rebuilt by `index_sym(idx)::Num` as `Sym{SymReal}(name; type=Int)` (plus the `IndexSlot` metadata when `slot != 0`). SymbolicUtils hashconsing makes the reconstruction the *same* object as the originally minted symbol, so substitution and `get_variables` are unaffected. The name-only `Sym` is the expensive part (a hashcons lookup, roughly 260 ns and 5 allocations), so it is **cached per name id** in `_SYM_BY_ID` (parallel to `_NAME_BY_ID`); `index_sym` on an abstract index is then a cached read. The cache is filled lazily on first use, *not* eagerly at intern time, so no `Sym` is baked into the precompile image. A baked `Sym` would be a stale duplicate of the new session's hashconsed canonical one and would break the `===` guarantee. An anonymous concrete site (`name_id == 0`, `slot == k`) reconstructs to the integer `Num(k)`. `to_numeric`'s indexed sites path uses this anonymous form *only* for the coefficient substitution (`ω(i) → ω(k)`); the resolved *operator* index keeps its real name (with `slot == k`), so resolved ops still print, satisfy `has_index`, and match a `d` override key.

`index_name(idx)::Symbol` resolves the name; `index_slot(idx)` reads the `slot` directly (the sym-metadata `index_slot(x)` method is kept for back-compat). Equality/hash compare interned ids and exclude `slot` (just as the old equality excluded `sym`); ordering uses the lexicographic name-rank table.

**Why interned ids?** Measured: interning the name to an `Int32` (vs storing bytes inline via `InlineStrings`) is what shrinks `Index`/`Op` enough to realize the performance win; bytes-inline reaches `isbits` but not the smaller struct. The tables are written only at construction (cold path) under a `ReentrantLock`, and read bounds-checked on the hot path, so a stale/out-of-range id throws `BoundsError` rather than corrupting memory. Construction is the only writer and is *not* thread-safe: populate the tables from one thread before any concurrent canonicalization (the lock guards writers against each other, not against the lock-free readers). Ids are insertion-order and are never serialized (not portable across sessions).

**`change_index(expr, from, to)`** performs symbolic substitution, replacing the index `from` with `to` throughout an expression tree (operator indices and the reconstructed `index_sym` of symbolic prefactors). Used for diagonal splitting and renaming sum indices.

**`NotIdentical` metadata.** `DoubleIndexedVariable(:g, i, j; identical=false)` creates `g(i,j)` with metadata `NotIdentical = true`. Two mechanisms enforce `g(i,i) = 0`: (1) at construction time, if `i == j` the function immediately returns `Num(0)`; (2) after a later substitution makes both arguments equal (e.g. `change_index` with `i → j`), `_check_not_identical` detects the equality and returns zero. Together these enforce `g(i,i) = 0` for off-diagonal coupling constants.


## Averaging

`average(expr)` converts operator expressions into symbolic scalars (SymbolicUtils `Term` nodes):

```julia
struct AvgFunc end                        # singleton callable, the "operation"
const sym_average = AvgFunc()

_average(op) = SymbolicUtils.Term{SymbolicUtils.SymReal}(sym_average, [op]; type = _avg_symtype(op))
```

**Why a custom `AvgFunc` instead of a `Sym`?** Using a custom struct lets us define `SymbolicUtils.show_call` for `⟨...⟩` display without type piracy on SymbolicUtils symbols.

**Hermitian averages are `Real`, everything else is `Number`.** A moment `⟨A⟩` is real iff `A` is Hermitian, so `_avg_symtype(op) = isequal(adjoint(op), op) ? Real : Number` is a structural test (both sides are canonical). The `Real` upgrade buys a faster `simplify` (its rewrite rules take the cheaper real path) and self-conjugation through `_conj_atom` (`conj(⟨A⟩) == ⟨A⟩`, versus a `Number` leaf that stays wrapped and relies on `inner_adjoint`). The `+`/`*` path is unaffected (it dispatches on the vartype `SymReal`, already used). Typing is per-leaf, so a distributed sum like `a + a'` → `⟨a⟩ + ⟨a'⟩` stays `Number`; `sym_sum` nodes inherit `Real` from a real body, and `make_time_dependent` lifts a Hermitian average to a `Real` unknown.

Because `substitute`/rewriting rebuild through `TermInterface.maketerm`, whose default `type` comes from `promote_symtype` (which cannot see Hermiticity behind the opaque operator symtype), the average and sum nodes define `maketerm` overrides that re-derive the symtype from the operator/body *value*; without them a `substitute` would silently revert a `Real` leaf to `Number`. `_avg_symtype_from_arg` falls back to `Number` when the argument is not a bare `QField` (e.g. a `conj`-wrapped operator from `qadjoint`).

### Why `average` returns `BasicSymbolic`, not `Num`

Symbolics.jl conventionally wraps public-API outputs in `Symbolics.Num`, a `<: Real` adapter over the raw `SymbolicUtils.BasicSymbolic` tree. `average` deliberately breaks that convention because `Num` wrapping is not type-stable across `average(::QAdd)`: pulling c-number prefactors out of the average routes the internal arithmetic through SymbolicUtils promotion rules, which produce `BasicSymbolic`, and a `QAdd` that contains a symbolic-prefactor term then disagrees in wrapper type with a `QAdd` that does not.

Returning `BasicSymbolic` uniformly across `average(::QSym)`, `average(::QAdd)`, and `average(::BasicSymbolic)` (scalars pass through unchanged) keeps the return type stable, lets every downstream function (`is_average`, `acts_on`, `get_indices`, `numeric_average`, `undo_average`) dispatch on one type without unwrapping, and avoids fragile re-wrapping. `Num`-accepting dispatches exist as defensive fallbacks for callers (or Symbolics arithmetic) that re-wrap.

**Indexed sums become a dedicated node.** When averaging an indexed `QAdd`, each index-dependent term is wrapped in a moment-layer sum node `sym_sum(body, scope)` (operation `SumFunc`, mirroring the `AvgFunc` average node) rather than stamped with metadata. The `body` is the term's averaged scalar (numeric prefactor and multi-factor c-number coefficient included), so an index-dependent coefficient stays *inside* the sum scope; the previous metadata-on-leaf form could not express ``\sum_i (\text{numeric} \cdot f_i \cdot g_i)`` because `Add`/`Mul` discard metadata on any composite or numerically-scaled node. Terms sharing the same `(indices, ne)` are grouped into one node with a multi-term body, mirroring display grouping (`_group_dep_terms`).

The scope rides as a `SumScope` *argument* of the `Term`, not as metadata, because SymbolicUtils `isequal`/`hash` ignore metadata: two sums over the same body but different scope would otherwise be `isequal` and wrongly cancel in a subtraction. `SumScope` defines value-based `==`/`isequal`/`hash`; SymbolicUtils stores it as an opaque scalar argument, the same way an `AvgFunc` node stores its `QField` (a bare `Vector` argument would instead be array-ified into a symbolic array, losing the `Index` objects). `is_indexed_sum` recognises the node, `get_sum_indices`/`get_sum_non_equal` (and the retained alias `has_sum_metadata`) read the scope, and `undo_average` rebuilds the summed `QAdd` by folding the node's `ne` into each term and attaching its `indices`. Every consumer that walks averaged trees (`get_indices`, `acts_on`, `make_time_dependent`, `inner_adjoint`, `numeric_average`) recurses into `body` and reconstructs the node around the rewritten body, so the scope is preserved without special-casing each operation.

**No `Complex{Num}` intermediates in the result.** `average(::QAdd)` walks `(term, c::CNum)` pairs and folds each one into a single result. It must never let `result` become a `Complex{Num}`, because `SymbolicUtils.unwrap(::Complex{<:Num})` materialises a literal `Term(complex, [re, img])` node whenever both `re` and `img` are symbolic, and that node is opaque to `simplify` / `expand` and also generates a runtime `complex(::Real, ::Complex)` call when consumers like ModelingToolkit codegen the equation (Base defines no such method). Each per-term contribution is unwrapped to `BasicSymbolic{SymReal}` (the single concrete Moshi type backing `BasicSymbolic`) at the boundary, so the accumulator and the grouped sum bodies are concretely typed rather than `Any`: every contribution reduces to `SymReal` because the averaged leaf, `Symbolics.IM`, and the `CNum` real/imaginary parts are all real-symtype. The imaginary unit enters the chain through `Symbolics.IM`, the `BasicSymbolic{SymReal}` symbol that prints as `im`, instead of through `1im::Complex{Bool}`. Concretely, the constant branch contributes `r + i * Symbolics.IM` and the operator branch contributes `r * avg + i * Symbolics.IM * avg`, skipping the `r` or `i` add when zero. This keeps every result a clean polynomial in `Symbolics.IM` that downstream `simplify(...; expand = true)` can reduce to zero on identically-equal differences.

**`undo_average` always returns `QAdd`.** It accepts both `BasicSymbolic` and `Num` inputs (unwrapping `Num` first). Scalars become single-term `QAdd`s with an empty operator sequence, and lone `QSym`s become single-term `QAdd`s with unit prefactor. This uniform return type makes `undo_average` type-stable despite walking SymbolicUtils expression trees (where `operation(x)` returns `Any`).

**Singleton single-op terms average to the bare op.** `average(::QAdd)` wraps a term in `_single_qadd(_CNUM_ONE, ops, ne)` only when cross-op canonicalization needs the constraint: `length(term.ops) > 1`. A single-op term averages to `only(term.ops)` regardless of `term.ne`, because the enclosing scope now rides on the `sym_sum` node rather than on the average leaf. This relies on the dead-NE invariant from *QAdd internals*: any `ne` reaching `average` references some op index or sum-scope index, but for a single-op term there is no second op to constrain, so a leaf-level `ne` would be observationally inert and would only obstruct downstream dict lookups (notably QC's MTK variable map).


## Numeric conversion

The backend-neutral core in `src/numeric/` never names a concrete numeric type.
QuantumOpticsBase and QuantumToolbox support lives in package extensions under `ext/`, so
the symbolic algebra can be loaded without either package.

### [Adding a numeric backend](@id numeric-backend-interface)

Define a concrete singleton subtype of `NumericBackend` and extend the exported hooks. A
static-only backend needs the following methods:

| Hook | Responsibility |
|:--|:--|
| `numeric_basis(be, h, dims)` | Build the full backend basis or dimension descriptor. |
| `numeric_num_subsystems(be, basis)` | Return the number of tensor factors. |
| `numeric_subbasis(be, basis, slot)` | Select one tensor factor. |
| `numeric_operator(be, op, subbasis)` | Convert a symbolic leaf. |
| `numeric_embed(be, basis, slot, leaf)` | Place a leaf in the full space. |
| `numeric_identity(be, basis)` | Build the full-space identity. |
| `numeric_assemble(be, basis, terms)` | Combine `(ComplexF64, factors)` terms. |
| `numeric_materialize(be, assembled, op_type)` | Select the public representation. |

`numeric_assemble` may return a lazy object. `numeric_materialize` must give `op_type` these
backend-neutral meanings:

- `nothing`: return the backend's ordinary eager operator;
- `identity`: return `assembled` unchanged;
- any other callable: a backend-defined explicit conversion.

The ordinary eager type need not be sparse; that is only the choice made by both bundled
backends. Materialization happens once, after the whole expression is assembled, so a leaf,
product, and sum have the same representation for a given `op_type`.

A minimal static backend looks like this:

```julia
import SecondQuantizedAlgebra as SQA

struct MyBackend <: SQA.NumericBackend end

SQA.numeric_basis(::MyBackend, h::SQA.HilbertSpace, dims) = error("build the basis")
SQA.numeric_num_subsystems(::MyBackend, basis) = error("count subsystems")
SQA.numeric_subbasis(::MyBackend, basis, slot::Int) = error("select a subsystem")
SQA.numeric_operator(::MyBackend, op::SQA.Op, subbasis) = error("construct a leaf")
SQA.numeric_embed(::MyBackend, basis, slot::Int, leaf) = error("embed a leaf")
SQA.numeric_identity(::MyBackend, basis) = error("construct the identity")
SQA.numeric_assemble(::MyBackend, basis, terms) = error("assemble terms")
SQA.numeric_materialize(::MyBackend, assembled, ::Nothing) =
    error("make an eager operator")
SQA.numeric_materialize(::MyBackend, assembled, ::typeof(identity)) = assembled
```

Optional capabilities add three small groups of methods:

| Capability | Additional hooks |
|:--|:--|
| time-dependent conversion | `numeric_assemble_td` |
| conversion from backend states | `numeric_backend(state)`, `numeric_basis(state)` |
| `numeric_average` / `expect` | the state hooks plus `numeric_expect` |

Third-party backends are passed explicitly as `backend=MyBackend()`; automatic discovery is
limited to the bundled extensions. `numeric_operator` covers the existing closed `OpKind`
roles. Packages can support a new backend/basis combination, but cannot currently define a
new symbolic operator role through this interface.

### Implementation notes

- `NumericContext` holds the concrete backend singleton, opaque basis, substitutions, and
  indexed-site map. The core emits `(ComplexF64, factors)` terms and performs no backend
  operator arithmetic.
- The internal assembled type must not depend on the number of terms. QuantumOptics uses
  the five-argument, vector-backed `LazySum`; QuantumToolbox uses the extension's
  vector-backed `VecSum` inside `QobjEvo`. Tests enforce inference stability.
- `numeric_average` calls `numeric_expect` on the lazy assembly directly. It therefore does
  not build an eager matrix merely to compute an expectation value.
- A non-empty `time_parameter` uses `numeric_assemble_td` and returns the native
  `TimeDependentSum` or `QobjEvo`. Only `op_type=nothing` and `identity` are accepted because
  a time-varying operator cannot be materialized once during conversion.
- **p-aware coefficients (differentiable control).** A `time_parameter` value may be a
  two-argument function `(p, t) -> value` of the solver parameter vector `p`. Detection is by
  arity: if any value function reads `p`, the conversion is p-aware and all value functions are
  lifted to `(p, t)` form. `p` stays opaque to the symbolic layer (a symbol maps to a value
  function; whether it consults `p` is that function's business), so the symbolic prefactor
  machinery is untouched. The QuantumToolbox backend threads `p` into the coefficient, giving a
  `QobjEvo` differentiable with respect to `p` (SciMLSensitivity with Enzyme/Mooncake, not
  Zygote). QuantumOptics has no solver `p` and rejects the form; the asymmetry is inherent to
  the backends.
- **`mesolve` superoperator path.** QuantumToolbox's superoperator builders fall back to a
  generic lazy tensor for operator types they do not recognize, which is not AD-traversable and
  emits a warning whose `try/catch` breaks Mooncake. The QuantumToolbox extension therefore
  teaches the `VecSum` to participate directly, so `liouvillian` builds a `VecSum` superoperator
  of concrete matrices that stays differentiation-friendly.
- Indexed conversion validates sites, unrolls sums, and emits the same term format as the
  non-indexed path. A site must be unique and within the backend basis.
- QuantumToolbox levels are zero-based internally, so symbolic `Transition(i, j)` maps to
  `projection(N, i-1, j-1)`.
- Lazy Hamiltonians are useful when materializing a large tensor product is prohibitive,
  but can slow superoperator construction. Prefer the eager default for small and medium
  solver problems; request `op_type=identity` deliberately.


## Hermitian conjugation (operators.jl)

Two exported helpers handle Hermitian conjugation on mixed operator/symbolic expression trees, including averaged operator nodes (`⟨...⟩`):

- **`qadjoint(x)`** (aliased as `qconj` and `dagger`): Hermitian conjugate that distributes through `SymbolicUtils.BasicSymbolic` trees and dispatches to `adjoint` on `QField` and `Number`. Distinct from `Base.conj`, which on a `BasicSymbolic` returns an opaque `conj(...)` wrapper instead of recursing into arguments; the distributed form is needed by downstream hashing and substitution machinery.
- **`inner_adjoint(x)`**: Pushes the adjoint *inside* `AvgFunc` nodes, rewriting `conj(⟨X⟩)` as `⟨X†⟩`. Used when building equations of motion where both sides must share the canonical "average-of-operator" form. Also collapses nested `conj(avg(...))` by recursing into the argument.

These cannot be wired directly to `Base.conj`/`Base.adjoint` on `SymbolicUtils.BasicSymbolic` because that would be type piracy — defining methods on `Base` functions for a type we don't own. Downstream packages (QuantumCumulants.jl, QuantumInputOutput.jl) call these explicitly.


## Ordering keys: three distinct orders

SQA carries three orderings that must not be conflated, because each answers a different question.

**`_site_compare` (partial, commutation).** The `SiteCmp` three-way comparator drives `_partial_sort!` during normal ordering. It is deliberately *not* identity-faithful: it returns `Equal` for operators that differ only in `axis` (Pauli/Spin) or levels `i,j` (Transition/CollectiveTransition), because it encodes which adjacent factors may be reordered, not whether two operators are the same. It is also partial (`Undetermined` for free indices with no resolving `ne`).

**`_full_op_key`/`_sort_key` (display sort).** Used by `sorted_arguments` for deterministic *display* order only. It carries the packed level/axis fields, without which two axes of one spin triple tie and display order falls back to dict iteration order. Still not a substitute for `order_key`, which is the identity-faithful one.

**`order_key`/`term_order_key`/`qadd_order_key` (total, identity-faithful).** The public ordering used by downstream packages (QuantumCumulants.jl) to pick canonical representatives and compare expressions reproducibly. The contract is `order_key(a) == order_key(b)` iff `isequal(a, b)`: the key carries every identity field, so it is a strict total order that never ties distinct operators. `order_key(o::Op)` is a single method returning a uniform tuple `(space_index, kind, name, index_key, l1, l2, g, nlev)`; the leading `space_index` groups operators by subspace, `kind` gives the cross-family order within a subspace, and the packed `l1..nlev` keep distinct levels/axes distinct (zero for roles that do not use them). It is preferred over hashing for ordering because Julia's `hash` is not stable across versions/platforms (which would make canonical-representative choice, and therefore generated equations, irreproducible) and collides; hashing remains correct for equality-only dedup via `Dict`/`Set`, which use `hash` plus `isequal`.

## Printing and LaTeX

**Terminal printing** uses Unicode: `†` for dagger, subscript digits (`₀`-`₉`) for Transition and CollectiveTransition levels, `σx`/`σy`/`σz` for Pauli axes. Summations render as `Σ(i=1:N)`.

**`sorted_arguments`** ensures deterministic output order. The sort key is `(length(ops), full_op_keys...)` where `_full_op_key(op) = (_sort_key(op)..., _type_order(op), _name_rank(op.name_id), op.l1, op.l2, op.g, op.nlev)`. This gives: shorter terms first, then by site, then by type (Destroy < Create < Transition < Pauli < Spin < Position < Momentum < CollectiveTransition), then by name.

**LaTeX** uses Latexify.jl's `@latexrecipe` macro. `transition_superscript(::Bool)` toggles the global `transition_idx_script` `Ref` between `:^` and `:_`, controlling whether Transition and CollectiveTransition level indices render as superscripts (`{name}^{{ij}}`) or subscripts (`{name}_{{ij}}`).
