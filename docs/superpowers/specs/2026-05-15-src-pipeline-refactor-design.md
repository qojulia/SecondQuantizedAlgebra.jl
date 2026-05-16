# src/ Rewrite: Ground-Up Algebra Redesign

**Date:** 2026-05-15
**Branch base:** `redesign-v2`
**Target:** new branch, big-bang rewrite

## Goal

Rewrite the operator algebra from first principles. Replace the current tangled worklist in `src/arithmetics/` with a streaming per-term pipeline of small, named, contract-bearing passes. Eliminate the `phys_ops` shadow field and the eager ground-state completeness expansion that together account for the algebra's two historical correctness bugs. End with a codebase where every file fits on one screen, every function has one job, and the whole `*` data flow is readable top-to-bottom in one short file.

## Architectural principle

Both historical bugs share one root cause: **applying a rewrite eagerly, before the surrounding context provides the information that determines whether the rewrite gives the right answer.**

- The original `phys_ops` substitution bug: site-sorting destroys physical adjacency that later index substitution needs. The fix was to memo the discarded order in a parallel `phys_ops` field. That patch outlived its diagnosis.
- The PR #99 dissipator bug: ground-state expansion `σᵍᵍ → 1 - Σσᵏᵏ` distributes a single atom into N+1 terms, entangling local algebra with global summation machinery. Cancellations that should happen at the `σᵍᵍ`-atom level fire too late, or not at all.

The unifying fix is one architectural rule:

> **Do not apply a rewrite when the information needed to apply it correctly is undetermined.**

Three categories of rewrite follow:

1. **Unconditionally safe** — fire eagerly. Local reductions on operators provably on the same site (same space, syntactically-equal index).
2. **Safe when site relationship is known** — fire eagerly if known, defer otherwise. Sorting two operators by site is safe iff their sites cannot merge under later substitution.
3. **Mathematically exact but structurally destructive** — fire only on user request. The completeness identity `σᵍᵍ = 1 - Σ σᵏᵏ` is true but turns one atom into N+1, defeating like-term collection.

The rewrite implements this rule by storing ops in a partial-canonical form that respects undetermined-site relationships, and by exposing completeness expansion as an opt-in public function.

## Constraints

- All existing `test/*.jl` pass after the rewrite. Tests that asserted eager GS expansion are updated to call the explicit `expand_completeness` (intended API change). All others pass unchanged.
- Public exports preserved: `normal_order`, `simplify`, `commutator`, `substitute`, operator constructors, Hilbert space constructors.
- One new public export: `expand_completeness(expr)`.
- One internal removal visible at the type level: `QTerm.phys_ops` field. No user code in the repo reaches into this field directly.
- `make benchlocal` medians and allocation counts within noise of the `redesign-v2` baseline.
- All quality gates green: JET zero issues, Aqua clean, ExplicitImports clean, CheckConcreteStructs clean.

## Out of scope

- Operator semantics. The algebra rules per operator type are ported verbatim.
- Hilbert space types. `FockSpace`, `NLevelSpace`, `PauliSpace`, `SpinSpace`, `PhaseSpace`, `ProductSpace` unchanged.
- Numeric layer. `to_numeric`, `numeric_average` unchanged.
- Average layer. `average`, `undo_average` unchanged.
- Printing. `printing.jl` and `latexify_recipes.jl` unchanged. The visible printing difference (symbolic-indexed ops in physical order) falls out of the storage change.
- Performance work beyond preserving the current bar.
- Public API additions other than `expand_completeness`.

---

## Storage: QTerm

```julia
struct QTerm
    ops::Vector{QSym}            # canonical partial-sort form
    ne::Vector{NonEqualPair}     # pairwise α ≠ β constraints
end
```

Gone: the `phys_ops::Vector{QSym}` field, the `_uses_phys_key` predicate, the conditional hash/equality logic.

Hash and equality are direct:

```julia
Base.hash(t::QTerm, h::UInt) = hash(t.ops, hash(t.ne, h))
Base.isequal(a::QTerm, b::QTerm) = isequal(a.ops, b.ops) && isequal(a.ne, b.ne)
```

Two `QTerm`s are equal as dict keys iff their `(ops, ne)` are structurally equal. The canonical-form rule below guarantees that two products which differ only by permutations of provably-commuting operators produce identical `ops` vectors.

## Canonical form: the partial-sort rule

`ops` is sorted by site, with one exception: **operators whose pairwise site relationship is undetermined keep their physical order with one another.**

Site relationship between two operators `a` and `b`:

- **Distinct** — different `space_index`; or different numeric indices on the same space; or `(a.index, b.index)` ∈ `ne`.
- **Equal** — same `space_index` AND syntactically-equal `index`.
- **Undetermined** — same `space_index`, indices not syntactically equal, no `ne` constraint resolves them.

The comparator returns a typed enum:

```julia
@enum SiteCmp::UInt8 Less Equal Undetermined Greater

_site_compare(a::QSym, b::QSym, ne::Vector{NonEqualPair})::SiteCmp
```

- `Less` / `Greater` — distinct sites; the sort places them in canonical order.
- `Equal` — same site; sort never reorders these.
- `Undetermined` — sort leaves them in physical order with one another.

`Equal` and `Undetermined` are both "don't reorder," but distinguishing them in the type eliminates the "I returned 0 but meant Equal" class of bugs.

`_partial_sort!(ops, ne)` is the stable sort using this comparator. It is the only function that ever reorders `ops`.

Example: `[A_fock, σ²¹ᵢ, B_fock, σ¹²ⱼ]` with `i, j` symbolic on the same atom space and `(i,j) ∉ ne` sorts to `[A_fock, B_fock, σ²¹ᵢ, σ¹²ⱼ]`. The two fock operators move to the front (distinct sites from the σs); the two σs stay in their physical order (undetermined relationship with each other).

---

## The pipeline: streaming passes

Default canonicalization runs two passes. Each pass has signature `(ops::Vector{QSym}, c::CNum, sink::F) where {F}`. Passes don't see `ne` — they only act on provably-same-site adjacencies, which is a property of two operators, not of the surrounding term.

The `where {F}` on every sink parameter forces closure specialization. The nested do-block chain inlines into one fused function with zero closure allocation.

```julia
@inline function _stream!(out, ops::Vector{QSym}, c::CNum, ne::Vector{NonEqualPair})
    _commute_ops(ops, c) do ops1, c1
        _reduce_ops(ops1, c1) do ops2, c2
            _canonicalize_to_dict!(out, ops2, c2, ne)
        end
    end
end
```

Each pass starts with an O(length(ops)) "any work to do?" pre-scan and calls `sink(ops, c)` once if not.

### The two default passes

**`_commute_ops(ops, c, sink)`** — apply commutation swaps for adjacent operators on provably the same site that need them (e.g. `aa† → a†a + 1`). Operators with undetermined site relationship are not commuted. Forks emit multiple terms through `sink`.

**`_reduce_ops(ops, c, sink)`** — fold adjacent provably-same-site pairs via `_reduce_pair`. Single output per input (possibly with zero coefficient, in which case nothing is emitted).

**`_canonicalize_to_dict!(out, ops, c, ne)`** — terminal sink. Constructs `QTerm(ops, ne)`, looks up in `out`, sums `c`, drops zero-coefficient entries. The only place a `QTerm` value is constructed during a pipeline run.

Neither default pass applies the completeness identity. `σᵍᵍ` is a valid canonical-form atom.

### Opt-in passes

**`_expand_gs_ops(ops, c, sink)`** — apply `σᵍᵍ → 1 - Σ_{k≠g} σᵏᵏ` to every `σᵍᵍ` in `ops`. May fork by `n_levels` per `σᵍᵍ`.

**`_substitute_ops(ops, c, d, sink)`** — substitute operators per dict `d` (operator→scalar, operator→single-term QAdd, or symbolic-substitution dict for prefactors). Forks when an operator splices in a multi-op replacement.

### Pass contracts

- **`_commute_ops`** — Pre: `ops` is in canonical partial-sort order (post-`_partial_sort!`). Post: every emitted `ops` has no adjacent provably-same-site pair `(a, b)` for which `_can_commute(a, b)` is false. Undetermined-site and distinct-site adjacents are untouched.
- **`_reduce_ops`** — Pre: input satisfies `_commute_ops`'s post-condition. Post: every emitted `ops` has no adjacent provably-same-site pair `(a, b)` for which `_reduce_pair(a, b)` returns non-`nothing`.
- **`_canonicalize_to_dict!`** — Pre: `ops` is in canonical form (post-passes). Post: `out` contains no zero-coefficient entries; like terms combined; every dict key is a freshly constructed `QTerm(ops, ne)`.

The output of `_stream!` satisfies the **canonical-form invariant**: `ops` in canonical partial-sort order, no adjacent provably-same-site pair with a remaining commutation residual or reduction, like terms combined. The invariant says nothing about `σᵍᵍ` atoms — they are allowed.

### Aliasing rules

Forks emit more than one output term. Single-output passes don't.

- **Single-output passes** mutate `ops` in place; the sink receives the same Vector that was passed in.
- **Fork passes** call `copy(ops)` before mutating each branch except possibly one. **No two branches ever share a mutable `ops` Vector.**
- **The terminal sink** takes ownership of the `ops` it receives. Callers must not mutate after passing.

Each pass file documents which rule it follows.

---

## The canonicalize primitive

Every operation that perturbs a term's `ops` (multiplication concatenation, substitution, adjoint reversal, index relabeling) ends with the same primitive:

```julia
@inline function _canonicalize!(out, ops::Vector{QSym}, c::CNum, ne::Vector{NonEqualPair})
    _partial_sort!(ops, ne)
    _stream!(out, ops, c, ne)
    return nothing
end
```

One call: re-sort, then run the default passes, then insert into `out`. Used by `*`, `Base.adjoint`, `change_index`, `normal_order`, and the off-diagonal branch of diagonal splitting.

---

## Public entry points

`Base.:*` is implemented as four dedicated fast paths (see "Multiplication" below). The other entry points:

```julia
function normal_order(q::QAdd)
    out = _new_buffer(length(q))
    for (t, c) in q
        _stream!(out, copy(t.ops), c, t.ne)
    end
    QAdd(out, copy(q.indices))
end

simplify(q::QAdd) = normal_order(q)

commutator(a, b) = a * b - b * a

function expand_completeness(q::QAdd)
    out = _new_buffer(length(q))
    for (t, c) in q
        _expand_gs_ops(copy(t.ops), c) do ops1, c1
            _stream!(out, ops1, c1, t.ne)
        end
    end
    QAdd(out, copy(q.indices))
end

function substitute(q::QAdd, d)
    out = _new_buffer(length(q))
    for (t, c) in q
        _substitute_ops(copy(t.ops), c, d) do ops1, c1
            _stream!(out, ops1, c1, t.ne)
        end
    end
    QAdd(out, copy(q.indices))
end
```

`substitute` re-canonicalizes through `_stream!` because substituting an index may turn previously-undetermined pairs into determined ones, allowing further commutation and reduction to fire.

---

## Multiplication

Four operator-operator `Base.:*` overloads, all sharing one inner helper. Dedicated fast paths from day one — promotion-to-QAdd would add 2 Dict allocations per call on the `QSym` hot paths.

```julia
@inline function _emit_product!(
    out::QTermDict,
    ta_ops::Vector{QSym}, ca::CNum, ta_ne::Vector{NonEqualPair},
    tb_ops::Vector{QSym}, cb::CNum, tb_ne::Vector{NonEqualPair},
    sum_indices::Vector{Index}, needs_diag_split::Bool,
)
    n = length(ta_ops) + length(tb_ops)
    ops = Vector{QSym}(undef, n)
    copyto!(ops, 1, ta_ops, 1, length(ta_ops))
    copyto!(ops, length(ta_ops) + 1, tb_ops, 1, length(tb_ops))
    ne = _merge_ne(ta_ne, tb_ne)
    c  = _mul_cnum(ca, cb)
    if needs_diag_split
        _accumulate_with_diag!(out, ops, c, sum_indices, ne)
    else
        _canonicalize!(out, ops, c, ne)
    end
    return nothing
end

function Base.:*(a::QSym, b::QSym)
    out = _new_buffer(1)
    _emit_product!(out,
        QSym[a], _CNUM_ONE, _EMPTY_NE,
        QSym[b], _CNUM_ONE, _EMPTY_NE,
        _EMPTY_INDICES, false)
    QAdd(out, Index[])
end

function Base.:*(a::QAdd, b::QSym)
    out = _new_buffer(length(a))
    needs = !isempty(a.indices)
    for (ta, ca) in a
        _emit_product!(out, ta.ops, ca, ta.ne,
                            QSym[b], _CNUM_ONE, _EMPTY_NE,
                            a.indices, needs)
    end
    QAdd(out, copy(a.indices))
end

function Base.:*(a::QSym, b::QAdd)
    out = _new_buffer(length(b))
    needs = !isempty(b.indices)
    for (tb, cb) in b
        _emit_product!(out, QSym[a], _CNUM_ONE, _EMPTY_NE,
                            tb.ops, cb, tb.ne,
                            b.indices, needs)
    end
    QAdd(out, copy(b.indices))
end

function Base.:*(a::QAdd, b::QAdd)
    out = _new_buffer(length(a) * length(b))
    sum_indices = _merge_unique(a.indices, b.indices)
    needs = !isempty(sum_indices)
    for (ta, ca) in a, (tb, cb) in b
        _emit_product!(out, ta.ops, ca, ta.ne,
                            tb.ops, cb, tb.ne,
                            sum_indices, needs)
    end
    QAdd(out, sum_indices)
end
```

The scalar overloads (`QSym*Number`, `QAdd*Number`, etc.) and `Base.:^` chain through these unchanged from the current implementation.

---

## Diagonal splitting

Diagonal splitting handles indexed-summation collisions. The redesign operates on `ops` directly: substitution preserves physical adjacency because `ops` already preserves physical order for symbolic-indexed pairs (that is exactly what the partial-sort guarantees).

```julia
function _accumulate_with_diag!(
    out::QTermDict, ops::Vector{QSym}, c::CNum,
    sum_indices::Vector{Index}, ne::Vector{NonEqualPair},
)
    _canonicalize!(out, copy(ops), c, ne)         # off-diagonal contribution

    distinct = _distinct_op_indices(ops)
    length(distinct) < 2 && return nothing

    for sum_idx in sum_indices
        _depends_on_index_ops(c, ops, sum_idx) || continue
        for ext_idx in distinct
            ext_idx == sum_idx && continue
            ext_idx.space_index == sum_idx.space_index || continue
            _ne_contains(ne, sum_idx, ext_idx) && continue
            sub_ops = QSym[change_index(o, sum_idx, ext_idx) for o in ops]
            sub_c   = change_index(c, sum_idx, ext_idx)
            sub_ne  = _drop_ne_with(ne, sum_idx)
            _emit_diagonal!(out, sub_ops, sub_c, sub_ne, sum_idx, ext_idx)
        end
    end
    return nothing
end
```

### Helper signatures

Five helpers are ported from `src/arithmetics/qadd_arithmetic.jl` with `phys_ops` arguments surgically removed. Bodies otherwise unchanged.

```julia
_distinct_op_indices(ops::Vector{QSym})::Vector{Index}
_depends_on_index_ops(c::CNum, ops::Vector{QSym}, idx::Index)::Bool
_emit_diagonal!(
    out::QTermDict, ops::Vector{QSym}, c::CNum, ne::Vector{NonEqualPair},
    α::Index, β::Index,
)::Nothing
_drop_ne_with(ne::Vector{NonEqualPair}, idx::Index)::Vector{NonEqualPair}     # unchanged
_ne_contains(ne::Vector{NonEqualPair}, α::Index, β::Index)::Bool              # unchanged
```

---

## Operations on existing QAdds

### `Base.adjoint`

```julia
function Base.adjoint(q::QAdd)
    out = _new_buffer(length(q))
    for (t, c) in q
        rev = QSym[adjoint(o) for o in Iterators.reverse(t.ops)]
        _canonicalize!(out, rev, conj(c), t.ne)
    end
    QAdd(out, copy(q.indices))
end
```

One canonicalize call handles both the post-reverse re-sort and the pipeline.

### `change_index(q, from, to)`

Transform indices in `ops` and `c`, transform `ne` via `_substitute_ne`, then re-canonicalize. The diagonal-split path inside `change_index` (when relabeling onto an existing index) calls `_accumulate_with_diag!`.

### `substitute(q, dict)`

Detailed in "Public entry points." `_substitute_ops` handles operator substitution (with scalar or single-term QAdd values) and symbolic-coefficient substitution (`SymbolicUtils.substitute` on `c`) in one pass.

---

## Operator hooks

Each concrete `QSym` subtype (`Destroy`, `Create`, `Transition`, `Pauli`, `Spin`, `Position`, `Momentum`) implements five methods:

- `_can_commute(a, b)::Bool` — true iff swapping `a` and `b` requires no residual commutator. Called only on provably-same-site pairs.
- `_commute_pair(a, b) -> (QSym, QSym, CNum)` — for same-site non-commuting pairs, returns `(b, a, residual_coeff)`.
- `_reduce_pair(a, b) -> Union{Nothing, Tuple{QSym, CNum}, CNum}` — local algebraic identity for same-site pairs. `nothing` if none.
- `_ground_state_expand(op) -> Union{Nothing, Tuple{Int, Int, Int}}` — for a `Transition` `σᵍᵍ`, returns `(g, n_levels, site_index)`; else `nothing`. Only `Transition` overrides non-trivially.
- `_site_compare(a, b, ne)::SiteCmp` — three-way comparison driving the partial-sort. See "Canonical form" above.

### Cross-type fallbacks

Operators of different concrete types are always on distinct sites (every `QSym` subtype lives in its own Hilbert space variant). Generic fallbacks cover cross-type pairs:

```julia
_can_commute(::QSym, ::QSym)         = true
_commute_pair(::QSym, ::QSym)        = error("unreachable: cross-type _commute_pair")
_reduce_pair(::QSym, ::QSym)         = nothing
_ground_state_expand(::QSym)         = nothing
_site_compare(a::QSym, b::QSym, _)   = _type_order(typeof(a)) < _type_order(typeof(b)) ? Less : Greater
```

These five hooks are the entire interface between the algebra and operator types.

---

## Performance posture

- One Dict allocation per `*` (output buffer).
- One Vector per `QTerm` (down from two).
- Pre-sized `ops` allocation per pair (`Vector{QSym}(undef, n+m)` + `copyto!`).
- Zero intermediate term-list materialization (streaming pipeline).
- Zero closure allocation (`where {F}` specialization).
- Early-out pre-scans on each pass.
- `CNum` fast paths untouched (`expressions/cnum.jl` unchanged).
- Canonical-form sort cost paid once at term construction.

### Known regression risk: chained operations

`substitute → change_index → simplify` re-canonicalizes every term three times. Current code skips passes when known unnecessary; the rewrite always runs the full pipeline.

**Mitigation reserved, not shipped:** a `_known_canonical::Bool` flag on QTerm (or on the QAdd container) set to `true` after pipeline output. Re-canonicalization short-circuits on known-canonical input. Cost: 1 byte per term, one branch per call. Added only if benchmarks demand it.

---

## File layout

```
src/SecondQuantizedAlgebra.jl       # module + exports (adds expand_completeness)
src/types.jl                        # QField, QSym abstract types
src/precompile.jl

src/expressions/
  cnum.jl                           # unchanged
  index_types.jl                    # unchanged
  index.jl                          # unchanged (change_index dispatch lives here)
  qterm.jl                          # SIMPLIFIED — one Vector, direct hash/eq
  qadd.jl                           # mostly unchanged

src/operators/                      # each gains five hook methods
  hilbertspace.jl
  fock.jl
  nlevel.jl
  pauli.jl
  spin.jl
  phase_space.jl
  operators.jl

src/algebra/                        # NEW — replaces arithmetics/
  canonical_sort.jl                 # SiteCmp enum, _site_compare driver, _partial_sort!
  commute.jl                        # _commute_ops pass
  reduce.jl                         # _reduce_ops pass
  expand_completeness.jl            # _expand_gs_ops pass + public expand_completeness
  substitute.jl                     # _substitute_ops pass + public substitute
  canonicalize.jl                   # _canonicalize_to_dict! terminal sink
  weyl.jl                           # normal ↔ symmetric ordering (ported as-is)
  pipelines.jl                      # _stream!, _canonicalize!, _emit_product!,
                                    # _accumulate_with_diag!, Base.:*, Base.adjoint,
                                    # normal_order, simplify, commutator
  README.md                         # pass-debugging idiom (see "Debugging a pass")

src/average.jl                      # unchanged; split into src/average/ if >200 LOC
src/numeric.jl                      # unchanged

src/printing/
  printing.jl                       # unchanged
  latexify_recipes.jl               # unchanged
```

Deleted: the entire `src/arithmetics/` directory.

Target: every source file under 200 LOC.

---

## Debugging a pass

Pass-as-callback nesting with `where {F}` specialization produces dense stack traces. To make passes debuggable, each pass file documents its contracts at the top, and the development idiom for inspection is:

```julia
using SecondQuantizedAlgebra: _commute_ops, CNum, _CNUM_ONE

# Vector-pushing sink for inspection.
emitted = Tuple{Vector{QSym}, CNum}[]
_commute_ops(ops_input, _CNUM_ONE) do ops, c
    push!(emitted, (copy(ops), c))   # copy because of aliasing rules
end
display(emitted)
```

`src/algebra/README.md` documents this idiom plus the contracts table for every pass.

---

## Migration steps

Single new branch off `redesign-v2`, committed as a sequence of reviewable steps. Each step ends green on all verification gates before the next begins.

1. **QTerm storage** — drop `phys_ops` field, `_uses_phys_key`, conditional hash/eq, `_term_key` phys_ops arg, `_addto!` phys_ops arg, `_copy_key` phys_ops field. Test suite breaks at this point; subsequent steps fix it.

2. **Operator hooks** (additive) — implement the five hooks per operator type by porting logic from `_try_reduction!` / `_try_swap!` in `src/arithmetics/ordering.jl`. Define cross-type fallbacks in `operators/operators.jl`. Old worklist still runs at this point; hooks dormant.

3. **`canonical_sort.jl`** — define `@enum SiteCmp`, `_type_order`, `_site_compare` dispatching to per-type hooks, `_partial_sort!`. Unit tests against hand-built operator vectors with mixed site relationships.

4. **Scaffold `src/algebra/`** — create empty stub files for each pass; wire includes into `src/SecondQuantizedAlgebra.jl`.

5. **Implement passes** — `_commute_ops`, `_reduce_ops`, `_canonicalize_to_dict!`, `_expand_gs_ops`, `_substitute_ops`. Each gets unit tests using the Vector-pushing sink idiom against hand-constructed `(ops, c, sink)` cases. Aliasing rules documented and tested.

6. **Implement `pipelines.jl`** — `_stream!`, `_canonicalize!`, `_emit_product!`, `_accumulate_with_diag!`, the four `Base.:*` overloads, `Base.adjoint`, `normal_order`, `simplify`, `commutator`, `expand_completeness`. Switch the module's algebra entry points to these. Existing tests now run on the new pipeline.

7. **Update tests for behavior changes:**
   - `test/nlevel_test.jl:50` — `simplify(σ * σ')` no longer expands; either compare against atomic form or wrap with `expand_completeness`.
   - `test/nlevel_test.jl:52` — `simplify(σgg, h)` no longer takes `h`; rewrite as `expand_completeness(σgg)`.
   - `test/nlevel_test.jl:63-64` — adjust as above.
   - `test/normal_order_test.jl:85-89` — `normal_order(σ11, h)` becomes `expand_completeness(σ11)`.
   - `test/normal_order_test.jl:97-99` — same.
   - `test/integration_test.jl:545+` — verify; update any direct `QTerm` field accesses to drop `phys_ops`.

8. **Add `test/canonical_invariant_test.jl`** — assert the canonical-form invariant on outputs of `*`, `normal_order`, `simplify`, `commutator`, `substitute`, `expand_completeness` against hand-built and random inputs.

9. **Add PR #99 regression test** — Lindblad dissipator on `op = σ(2,2,j)` with jump `J = σ(2,1,i)` produces `R - (R+γ)·σ(2,2,j)`, etc. The full test matrix from ChristophHotter's bug report.

10. **Move Weyl conversion** to `algebra/weyl.jl`. Implementation unchanged.

11. **Delete `src/arithmetics/`** — `ordering.jl`, `normal_order.jl`, `simplify.jl`, `qadd_arithmetic.jl`, `commutator.jl`, `substitute.jl`. Update includes.

12. **Final cleanup** — update precompile workload signatures, refresh exports list, write `src/algebra/README.md`, `make format`. Split `average.jl` into `src/average/{types.jl,transform.jl}.jl` only if it remains >200 LOC.

---

## Verification gates

Every migration step must satisfy all of these before being claimed complete:

- `make test` — all `test/*.jl` pass.
- `make format` — clean.
- JET (`jet_test.jl`) — zero issues.
- Aqua (`aqua_test.jl`) — clean.
- ExplicitImports (`explicit_imports_test.jl`) — clean.
- CheckConcreteStructs (`concrete_test.jl`) — passes.
- `make benchlocal` — medians and allocation counts within noise of the `redesign-v2` baseline.
- `test/canonical_invariant_test.jl` (added step 8) — passes.
- `test/dissipator_regression_test.jl` (added step 9) — passes.

The final merge to `redesign-v2` is gated on all gates green plus a holistic benchmark comparison.

---

## Behavior changes (API-visible surface)

Exhaustive list of every visible change:

1. **`σᵍᵍ` is no longer expanded by `*`, `normal_order`, or `simplify`.** Use `expand_completeness(expr)` to materialize the completeness identity. Fixes PR #99.
2. **`expand_completeness` added to exports.**
3. **`normal_order(expr, h)` and `simplify(expr, h)` no longer exist.** The `h::HilbertSpace` argument was the opt-in for GS expansion; that operation is now `expand_completeness(expr)`.
4. **`simplify` and `normal_order` are aliased.** They were semantically identical before but separately defined; now they're the same function.
5. **Operators with undetermined-site relationship stay in physical order in stored `ops`.** Visible in printing of indexed expressions. The dict-key identity is unchanged (canonical form already accounted for this via `phys_ops`).
6. **`QTerm.phys_ops` field removed.** No user code in the repo reads it.

All other public behavior is identical.
