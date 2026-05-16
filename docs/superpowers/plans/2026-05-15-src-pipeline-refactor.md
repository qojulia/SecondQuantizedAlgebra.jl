# src/ Pipeline Refactor Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use `superpowers:subagent-driven-development` (recommended) or `superpowers:executing-plans` to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.
>
> **Commits are user actions.** Per project CLAUDE.md, Claude must not run `git commit`, `git push`, or any history-modifying git command. Each task ends with a "Pause" step that hands control back to the user to review and commit.

**Goal:** Rewrite the operator algebra in `src/` as a streaming per-term pipeline of small named passes, dropping `phys_ops` and deferring ground-state expansion. End with every file <200 LOC and a clear top-to-bottom data flow for `*`.

**Architecture:** `QTerm` stores one `Vector{QSym}` in partial-canonical form (distinct sites sorted, undetermined sites stay in physical order). Multiplication concatenates physical ops, runs `_partial_sort!`, then streams through `_commute_ops → _reduce_ops → _canonicalize_to_dict!`. Ground-state expansion becomes opt-in via `expand_completeness(expr)`.

**Tech Stack:** Julia ≥ 1.10, `Pkg.test`, JET.jl, Aqua.jl, ExplicitImports.jl, CheckConcreteStructs.jl, BenchmarkTools.jl, Runic.jl (formatting).

**Spec:** [docs/superpowers/specs/2026-05-15-src-pipeline-refactor-design.md](../specs/2026-05-15-src-pipeline-refactor-design.md)

---

## Pre-flight

### Task 0: Branch and benchmark baseline

**Files:**
- Read: `.git/HEAD`, `redesign-v2` is current branch
- Run: `make benchlocal`, save output to `benchmark/baseline-redesign-v2.txt`

- [ ] **Step 1: Verify current branch and clean tree**

Run:
```bash
git status
git branch --show-current
```
Expected: clean tree, on `redesign-v2`.

- [ ] **Step 2: Pause — user creates the rewrite branch**

User runs:
```bash
git checkout -b src-pipeline-rewrite
```

- [ ] **Step 3: Capture benchmark baseline**

Run:
```bash
make benchlocal 2>&1 | tee benchmark/baseline-redesign-v2.txt
```
Expected: all benchmarks complete; file written. This is the reference for "within noise" comparisons at every verification gate.

- [ ] **Step 4: Pause — user inspects baseline file, commits it**

User reviews `benchmark/baseline-redesign-v2.txt` and commits if they want it tracked.

---

## Migration Step 1: QTerm Storage Simplification

Goal: drop `phys_ops` from `QTerm` and the conditional hash/eq logic. Tests will fail at the end of this step — that is expected.

### Task 1.1: Simplify QTerm struct and hash/eq

**Files:**
- Modify: [src/expressions/qterm.jl:20-54](src/expressions/qterm.jl#L20-L54)

- [ ] **Step 1: Replace the `QTerm` struct definition**

Current `qterm.jl` lines 20-26:
```julia
struct QTerm
    ops::Vector{QSym}
    ne::Vector{NonEqualPair}
    phys_ops::Vector{QSym}
end

QTerm(ops::Vector{QSym}, ne::Vector{NonEqualPair}) = QTerm(ops, ne, copy(ops))
```

Replace with:
```julia
struct QTerm
    ops::Vector{QSym}
    ne::Vector{NonEqualPair}
end
```

(Delete the second constructor — there is now only one.)

- [ ] **Step 2: Delete `_uses_phys_key` (lines 28-42)**

The entire function:
```julia
function _uses_phys_key(term::QTerm)
    ...
end
```
gets removed.

- [ ] **Step 3: Replace `Base.isequal` and `Base.hash` (lines 44-54)**

Replace the conditional logic with the direct form:
```julia
Base.isequal(a::QTerm, b::QTerm) = isequal(a.ops, b.ops) && isequal(a.ne, b.ne)
Base.:(==)(a::QTerm, b::QTerm) = isequal(a, b)
Base.hash(term::QTerm, h::UInt) = hash(:QTerm, hash(term.ops, hash(term.ne, h)))
```

- [ ] **Step 4: Run module load to catch syntax errors**

Run:
```bash
julia --project=. -e 'using SecondQuantizedAlgebra'
```
Expected: errors complaining about `_uses_phys_key`, `phys_ops` references in other files. Note them — they're fixed in subsequent tasks.

- [ ] **Step 5: Pause**

User reviews the QTerm struct change before continuing.

### Task 1.2: Update `_term_key`, `_addto!`, `_copy_key`

**Files:**
- Modify: [src/expressions/qterm.jl:134-172](src/expressions/qterm.jl#L134-L172)

- [ ] **Step 1: Replace `_copy_key` (lines 134-136)**

Old:
```julia
function _copy_key(term::QTerm)
    return QTerm(copy(term.ops), _copy_ne(term.ne), copy(term.phys_ops))
end
```

New:
```julia
function _copy_key(term::QTerm)
    return QTerm(copy(term.ops), _copy_ne(term.ne))
end
```

- [ ] **Step 2: Replace `_term_key` (lines 145-152)**

Old:
```julia
function _term_key(
        ops::Vector{QSym},
        ne::Vector{NonEqualPair} = _EMPTY_NE,
        phys_ops::Vector{QSym} = ops
    )
    actual_phys = _needs_phys_tracking(ops) ? phys_ops : ops
    return QTerm(copy(ops), _canonical_ne(ne), copy(actual_phys))
end
```

New:
```julia
function _term_key(
        ops::Vector{QSym},
        ne::Vector{NonEqualPair} = _EMPTY_NE,
    )
    return QTerm(copy(ops), _canonical_ne(ne))
end
```

(Delete `_needs_phys_tracking` if it exists in this file — verify with grep.)

- [ ] **Step 3: Replace `_addto!` (lines 166-172)**

Old:
```julia
function _addto!(
        d::QTermDict, ops::Vector{QSym}, c::CNum,
        ne::Vector{NonEqualPair} = _EMPTY_NE,
        phys_ops::Vector{QSym} = ops
    )
    return _addto_key!(d, _term_key(ops, ne, phys_ops), c)
end
```

New:
```julia
function _addto!(
        d::QTermDict, ops::Vector{QSym}, c::CNum,
        ne::Vector{NonEqualPair} = _EMPTY_NE,
    )
    return _addto_key!(d, _term_key(ops, ne), c)
end
```

- [ ] **Step 4: Grep for any remaining `phys_ops` in qterm.jl**

Run:
```bash
grep -n phys_ops src/expressions/qterm.jl
```
Expected: no matches. If any remain, address them.

- [ ] **Step 5: Pause**

### Task 1.3: Strip `phys_ops` from `qadd.jl` adjoint and helpers

**Files:**
- Modify: [src/expressions/qadd.jl:45-56,71-81](src/expressions/qadd.jl)

- [ ] **Step 1: Replace `_ordered_qadd` (lines 45-56)**

Old:
```julia
function _ordered_qadd(
        c::CNum, ops::Vector{QSym},
        ne::Vector{NonEqualPair} = _EMPTY_NE,
        phys_ops::Vector{QSym} = ops
    )
    _iszero_cnum(c) && return QAdd(QTermDict(), Index[])
    d = QTermDict()
    for t in _apply_ordering(c, ops)
        _addto!(d, t.ops, t.prefactor, ne, phys_ops)
    end
    return QAdd(d, Index[])
end
```

New (interim — still uses old `_apply_ordering`; will be replaced in Step 6):
```julia
function _ordered_qadd(
        c::CNum, ops::Vector{QSym},
        ne::Vector{NonEqualPair} = _EMPTY_NE,
    )
    _iszero_cnum(c) && return QAdd(QTermDict(), Index[])
    d = QTermDict()
    for t in _apply_ordering(c, ops)
        _addto!(d, t.ops, t.prefactor, ne)
    end
    return QAdd(d, Index[])
end
```

- [ ] **Step 2: Simplify `Base.adjoint` (lines 71-81)**

Old:
```julia
function Base.adjoint(q::QAdd)
    d = QTermDict()
    for (term, c) in q.arguments
        adj_ops = QSym[adjoint(op) for op in reverse(term.ops)]
        _site_sort!(adj_ops)
        adj_phys_ops = QSym[adjoint(op) for op in reverse(term.phys_ops)]
        _phys_sort!(adj_phys_ops)
        _addto!(d, adj_ops, conj(c), term.ne, adj_phys_ops)
    end
    return QAdd(d, copy(q.indices))
end
```

New (interim — uses old `_site_sort!`; will be replaced in Step 6 to route through `_canonicalize!`):
```julia
function Base.adjoint(q::QAdd)
    d = QTermDict()
    for (term, c) in q.arguments
        adj_ops = QSym[adjoint(op) for op in reverse(term.ops)]
        _site_sort!(adj_ops)
        _addto!(d, adj_ops, conj(c), term.ne)
    end
    return QAdd(d, copy(q.indices))
end
```

- [ ] **Step 3: Grep `qadd.jl` for any remaining `phys_ops`**

Run:
```bash
grep -n phys_ops src/expressions/qadd.jl
```
Expected: no matches.

- [ ] **Step 4: Pause**

### Task 1.4: Strip `phys_ops` from `qadd_arithmetic.jl`

**Files:**
- Modify: [src/arithmetics/qadd_arithmetic.jl](src/arithmetics/qadd_arithmetic.jl)

- [ ] **Step 1: Remove `phys_ops` arguments from `_addto!` call sites and helper signatures**

In `_apply_diag_terms!` (lines 28-36), drop the `phys_ops` parameter:
```julia
function _apply_diag_terms!(
        d::QTermDict, terms::Vector{OrderedTerm}, ne::Vector{NonEqualPair},
    )
    for t in terms
        _addto!(d, t.ops, t.prefactor, ne)
    end
    return d
end
```

In `_emit_diagonal!` (lines 77-97), drop `diag_phys_ops::Vector{QSym}` parameter and any uses. Update the body to remove `phys_ops` plumbing.

In `_accumulate_with_diag!` (lines 106-148), drop `phys_ops` setup (lines 111-112) and any reference. Drop the `phys_ops` argument from `_substitute_unsorted` (which is being deleted in step 11 anyway).

In `_substitute_unsorted` (lines 48-55), drop the trailing arguments related to phys_ops.

In all `Base.:*` overloads (`QSym*QSym` line 5-11, `QAdd*QSym` 150-160, `QSym*QAdd` 162-172, `QAdd*QAdd` 174-199), remove `unsorted`/`phys_ops` plumbing.

- [ ] **Step 2: Run tests to catch what's broken**

Run:
```bash
make test 2>&1 | head -80
```
Expected: tests fail. Note the failure patterns — most should be about missing methods or shape mismatches. The actual fix happens in Step 6.

- [ ] **Step 3: Pause**

User confirms the QTerm storage is simplified and the codebase loads (even if tests fail). This is the bottom of the migration — every subsequent step adds or replaces.

---

## Migration Step 2: Operator Hooks (Additive)

Goal: implement the five operator hooks per operator type. The old worklist still runs the algebra at this stage — the hooks are dormant until step 6. This step is purely additive; no behavior changes.

### Task 2.1: Define `SiteCmp` enum and cross-type fallbacks

**Files:**
- Create: [src/algebra/canonical_sort.jl](src/algebra/canonical_sort.jl) (new file, included later in Step 4; for now define the enum somewhere it can be referenced — put it in [src/types.jl](src/types.jl))
- Modify: [src/types.jl](src/types.jl)
- Modify: [src/operators/operators.jl](src/operators/operators.jl)

- [ ] **Step 1: Add `SiteCmp` enum to `src/types.jl`**

Append to `types.jl`:
```julia
"""
    SiteCmp

Three-way site comparison for the canonical partial-sort. `Equal` and
`Undetermined` both mean "do not reorder," but distinguishing them in the type
prevents the "I returned 0 but meant Equal" bug class.
"""
@enum SiteCmp::UInt8 Less Equal Undetermined Greater
```

- [ ] **Step 2: Add `_type_order` (per-concrete-type total ordering) in `src/operators/operators.jl`**

Append:
```julia
# Total ordering across QSym concrete types — used by _site_compare cross-type fallback.
_type_order(::Type{Destroy})    = 0
_type_order(::Type{Create})     = 1
_type_order(::Type{Transition}) = 2
_type_order(::Type{Pauli})      = 3
_type_order(::Type{Spin})       = 4
_type_order(::Type{Position})   = 5
_type_order(::Type{Momentum})   = 6
```

- [ ] **Step 3: Add cross-type hook fallbacks in `src/operators/operators.jl`**

Append:
```julia
# Generic fallbacks for cross-type operator pairs (always distinct sites).
_can_commute(::QSym, ::QSym)   = true
_commute_pair(::QSym, ::QSym)  = error("unreachable: cross-type _commute_pair")
_reduce_pair(::QSym, ::QSym)   = nothing
_ground_state_expand(::QSym)   = nothing

function _site_compare(a::QSym, b::QSym, ne::Vector{NonEqualPair})::SiteCmp
    ta = typeof(a); tb = typeof(b)
    ta === tb && error("unreachable: same-type _site_compare must be overridden by $ta")
    return _type_order(ta) < _type_order(tb) ? Less : Greater
end
```

- [ ] **Step 4: Load the module to verify no syntax errors**

Run:
```bash
julia --project=. -e 'using SecondQuantizedAlgebra'
```
Expected: loads cleanly. The new functions are unused at this point.

- [ ] **Step 5: Pause**

### Task 2.2: Implement Fock hooks (Destroy, Create)

**Files:**
- Modify: [src/operators/fock.jl](src/operators/fock.jl)

- [ ] **Step 1: Add Fock hooks**

Append to `fock.jl`:
```julia
# Same-type same-site site relationship is always Equal (same name/space/index).
# Different name or different space → distinct (sorted by name lex order, then space).
function _site_compare(a::Destroy, b::Destroy, ne::Vector{NonEqualPair})::SiteCmp
    a.space_index == b.space_index || return a.space_index < b.space_index ? Less : Greater
    a.name == b.name || return a.name < b.name ? Less : Greater
    a.index == b.index && return Equal
    # Different indices on same name/space — Fock ops with explicit indices are
    # not currently supported in mainline; conservative Undetermined return.
    return Undetermined
end
_site_compare(a::Create, b::Create, ne) = _site_compare(Destroy(a.name, a.space_index, a.index),
                                                          Destroy(b.name, b.space_index, b.index), ne)

# Destroy < Create within the same site (canonical normal order).
function _site_compare(a::Create, b::Destroy, ne::Vector{NonEqualPair})::SiteCmp
    cmp = _site_compare(Destroy(a.name, a.space_index, a.index),
                        Destroy(b.name, b.space_index, b.index), ne)
    cmp === Equal && return Greater   # Create > Destroy in canonical order
    return cmp
end
function _site_compare(a::Destroy, b::Create, ne::Vector{NonEqualPair})::SiteCmp
    cmp = _site_compare(Destroy(a.name, a.space_index, a.index),
                        Destroy(b.name, b.space_index, b.index), ne)
    cmp === Equal && return Less      # Destroy < Create in canonical order
    return cmp
end

# Same-site commutation: aa† does not commute (1 unit residual).
_can_commute(a::Destroy, b::Create) = false
_can_commute(a::Create,  b::Destroy) = true
_can_commute(a::Destroy, b::Destroy) = true
_can_commute(a::Create,  b::Create)  = true

# _commute_pair(a, b) returns (b, a, residual_coeff) — the swapped product
# plus a CNum factor times the contracted form (which is the "+1" branch).
# Convention: residual_coeff is the prefactor on the contracted-form output.
_commute_pair(a::Destroy, b::Create) = (b, a, _CNUM_ONE)   # aa† = a†a + 1·(identity branch)

# Reductions: ladder operators on the same site don't reduce locally.
# (a·a is not a closed-form simplification; only a·a† triggers the commutation residual.)
_reduce_pair(a::Destroy, b::Create) = nothing
_reduce_pair(a::Create,  b::Destroy) = nothing
_reduce_pair(a::Destroy, b::Destroy) = nothing
_reduce_pair(a::Create,  b::Create)  = nothing
```

- [ ] **Step 2: Pause**

### Task 2.3: Implement Transition hooks

**Files:**
- Modify: [src/operators/nlevel.jl](src/operators/nlevel.jl)

Look at [src/arithmetics/ordering.jl:43-55](src/arithmetics/ordering.jl#L43-L55) for the current Transition composition logic.

- [ ] **Step 1: Add Transition hooks**

Append to `nlevel.jl`:
```julia
function _site_compare(a::Transition, b::Transition, ne::Vector{NonEqualPair})::SiteCmp
    a.space_index == b.space_index || return a.space_index < b.space_index ? Less : Greater
    a.name == b.name || return a.name < b.name ? Less : Greater
    a.index == b.index && return Equal
    # Same name/space but different indices → check ne constraint.
    _ne_contains(ne, a.index, b.index) && return a.index < b.index ? Less : Greater
    return Undetermined
end

# Transitions on the same site never commute freely (composition fires).
_can_commute(a::Transition, b::Transition) = false

# Composition: σⁱʲ · σᵏˡ = δⱼₖ σⁱˡ.
function _reduce_pair(a::Transition, b::Transition)
    a.name == b.name || return nothing
    a.space_index == b.space_index || return nothing
    a.index == b.index || return nothing
    if a.j == b.i
        new = Transition(a.name, a.i, b.j, a.space_index, a.index, a.ground_state, a.n_levels)
        return (new, _CNUM_ONE)
    else
        return _CNUM_ZERO    # δ_{j,k} = 0
    end
end

# _commute_pair is never called for Transition (because _can_commute is always false
# and only out-of-order pairs would call it; here we use _reduce_pair instead).
_commute_pair(a::Transition, b::Transition) = error("unreachable: Transition uses _reduce_pair, not _commute_pair")

# Ground-state expansion: σᵍᵍ → 1 - Σ_{k≠g} σᵏᵏ.
function _ground_state_expand(op::Transition)
    op.i == op.ground_state && op.j == op.ground_state || return nothing
    return (op.ground_state, op.n_levels, op.space_index)
end
```

- [ ] **Step 2: Pause**

### Task 2.4: Implement Pauli hooks

**Files:**
- Modify: [src/operators/pauli.jl](src/operators/pauli.jl)

Look at [src/arithmetics/ordering.jl:58-67](src/arithmetics/ordering.jl#L58-L67) for the current Pauli composition logic and the Levi-Civita lookup at lines 11-12.

- [ ] **Step 1: Add Pauli hooks**

Append to `pauli.jl`:
```julia
# Per-site Levi-Civita lookup (3×3 antisymmetric). _levi_civita[j][k] = ε_{jk(6-j-k)}.
const _levi_civita = ((0, 1, -1), (-1, 0, 1), (1, -1, 0))

function _site_compare(a::Pauli, b::Pauli, ne::Vector{NonEqualPair})::SiteCmp
    a.space_index == b.space_index || return a.space_index < b.space_index ? Less : Greater
    a.name == b.name || return a.name < b.name ? Less : Greater
    a.index == b.index && return Equal
    _ne_contains(ne, a.index, b.index) && return a.index < b.index ? Less : Greater
    return Undetermined
end

_can_commute(a::Pauli, b::Pauli) = false   # always compose on same site

# σⱼ·σₖ = δⱼₖI + iϵⱼₖₗσₗ
function _reduce_pair(a::Pauli, b::Pauli)
    a.name == b.name || return nothing
    a.space_index == b.space_index || return nothing
    a.index == b.index || return nothing
    if a.axis == b.axis
        return _CNUM_ONE     # σⱼ² = 1
    else
        eps = _levi_civita[a.axis][b.axis]
        new = Pauli(a.name, 6 - a.axis - b.axis, a.space_index, a.index)
        return (new, _mul_cnum(_to_cnum(im * eps), _CNUM_ONE))
    end
end

_commute_pair(a::Pauli, b::Pauli) = error("unreachable: Pauli uses _reduce_pair")
```

- [ ] **Step 2: Pause**

### Task 2.5: Implement Spin hooks

**Files:**
- Modify: [src/operators/spin.jl](src/operators/spin.jl)

Look at [src/arithmetics/ordering.jl:143-149](src/arithmetics/ordering.jl#L143-L149) for the current Spin commutation logic.

- [ ] **Step 1: Add Spin hooks**

Append to `spin.jl`:
```julia
function _site_compare(a::Spin, b::Spin, ne::Vector{NonEqualPair})::SiteCmp
    a.space_index == b.space_index || return a.space_index < b.space_index ? Less : Greater
    a.name == b.name || return a.name < b.name ? Less : Greater
    a.index == b.index || (
        _ne_contains(ne, a.index, b.index) && return a.index < b.index ? Less : Greater
    )
    a.index == b.index || return Undetermined
    # Same site: canonical axis order is ascending (Sx=1 < Sy=2 < Sz=3).
    a.axis == b.axis && return Equal
    return a.axis < b.axis ? Less : Greater
end

# Same-site spins commute only when in canonical axis order (ascending).
_can_commute(a::Spin, b::Spin) = a.axis <= b.axis

# [Sj, Sk] = iϵⱼₖₗSl. The pair (Sj·Sk) with j>k commutes to (Sk·Sj + iε·Sl).
function _commute_pair(a::Spin, b::Spin)
    a.name == b.name || error("unreachable")
    a.space_index == b.space_index || error("unreachable")
    a.index == b.index || error("unreachable")
    a.axis > b.axis || error("unreachable: _commute_pair called on in-order pair")
    eps = _levi_civita[a.axis][b.axis]
    contracted = Spin(a.name, 6 - a.axis - b.axis, a.space_index, a.index)
    return (b, a, _mul_cnum(_to_cnum(im * eps), _CNUM_ONE))
end

_reduce_pair(::Spin, ::Spin) = nothing
```

(Note: the Spin design represents the commutator residual differently from Fock — the residual contracts to a single op, not the identity. The `_commute_pair` returns `(b, a, residual_coeff_for_contracted_op)`. The pipeline needs to handle this — clarified in `_commute_ops` design.)

- [ ] **Step 2: Pause**

### Task 2.6: Implement PhaseSpace hooks (Position, Momentum)

**Files:**
- Modify: [src/operators/phase_space.jl](src/operators/phase_space.jl)

Look at [src/arithmetics/ordering.jl:152-156](src/arithmetics/ordering.jl#L152-L156) for current PhaseSpace logic.

- [ ] **Step 1: Add PhaseSpace hooks**

Append to `phase_space.jl`:
```julia
function _site_compare(a::Position, b::Position, ne::Vector{NonEqualPair})::SiteCmp
    a.space_index == b.space_index || return a.space_index < b.space_index ? Less : Greater
    a.name == b.name || return a.name < b.name ? Less : Greater
    a.index == b.index ? Equal : Undetermined
end
_site_compare(a::Momentum, b::Momentum, ne) =
    _site_compare(Position(a.name, a.space_index, a.index),
                  Position(b.name, b.space_index, b.index), ne)

# Position < Momentum in canonical order.
function _site_compare(a::Position, b::Momentum, ne::Vector{NonEqualPair})::SiteCmp
    cmp = _site_compare(Position(a.name, a.space_index, a.index),
                        Position(b.name, b.space_index, b.index), ne)
    cmp === Equal && return Less
    return cmp
end
function _site_compare(a::Momentum, b::Position, ne::Vector{NonEqualPair})::SiteCmp
    cmp = _site_compare(Position(a.name, a.space_index, a.index),
                        Position(b.name, b.space_index, b.index), ne)
    cmp === Equal && return Greater
    return cmp
end

# Same-site P·X = X·P − i (commute residual is −i times identity).
_can_commute(a::Position, b::Momentum) = true
_can_commute(a::Momentum, b::Position) = false
_can_commute(a::Position, b::Position) = true
_can_commute(a::Momentum, b::Momentum) = true

_commute_pair(a::Momentum, b::Position) = (b, a, _to_cnum(-im))   # P·X = X·P − i·I

_reduce_pair(::Position, ::Momentum) = nothing
_reduce_pair(::Momentum, ::Position) = nothing
_reduce_pair(::Position, ::Position) = nothing
_reduce_pair(::Momentum, ::Momentum) = nothing
```

- [ ] **Step 2: Pause**

### Task 2.7: Smoke-test the hooks

**Files:**
- Create: `test/operator_hooks_test.jl` (temporary scratch test, deleted in step 12)

- [ ] **Step 1: Write a smoke test that exercises every hook**

Create `test/operator_hooks_test.jl`:
```julia
using Test
using SecondQuantizedAlgebra
using SecondQuantizedAlgebra: _site_compare, _can_commute, _commute_pair,
    _reduce_pair, _ground_state_expand, SiteCmp, Less, Greater, Equal, Undetermined,
    _CNUM_ONE, _CNUM_ZERO, _EMPTY_NE

@testset "Operator hooks smoke" begin
    h = FockSpace(:f)
    a = Destroy(h, :a)
    ad = adjoint(a)
    @test _site_compare(a, ad, _EMPTY_NE) === Less
    @test _site_compare(ad, a, _EMPTY_NE) === Greater
    @test _can_commute(a, ad) === false
    @test _can_commute(ad, a) === true
    sw = _commute_pair(a, ad)
    @test sw[1] === ad && sw[2] === a

    ha = NLevelSpace(:atom, 2)
    σ12 = Transition(ha, :σ, 1, 2)
    σ21 = Transition(ha, :σ, 2, 1)
    @test _site_compare(σ12, σ21, _EMPTY_NE) === Equal
    @test _can_commute(σ12, σ21) === false
    red = _reduce_pair(σ21, σ12)   # σ²¹·σ¹² = σ²²
    @test red isa Tuple
    @test red[1].i == 2 && red[1].j == 2
    @test _ground_state_expand(Transition(ha, :σ, 1, 1)) !== nothing
end
```

- [ ] **Step 2: Run the smoke test**

Run:
```bash
julia --project=. test/operator_hooks_test.jl
```
Expected: all tests pass.

- [ ] **Step 3: Run the main test suite**

Run:
```bash
make test 2>&1 | tail -30
```
Expected: many failures from Step 1 (QTerm changes). The hooks themselves are dormant.

- [ ] **Step 4: Pause**

User reviews the operator hooks; they're additive (old worklist still runs).

---

## Migration Step 3: Canonical Sort

### Task 3.1: Create `src/algebra/` and implement `_partial_sort!`

**Files:**
- Create: `src/algebra/canonical_sort.jl`

- [ ] **Step 1: Create the algebra/ directory**

Run:
```bash
mkdir -p src/algebra
```

- [ ] **Step 2: Write `src/algebra/canonical_sort.jl`**

```julia
"""
    _partial_sort!(ops::Vector{QSym}, ne::Vector{NonEqualPair})

In-place stable partial-sort using `_site_compare`. Distinct-site adjacents
are placed in canonical order. `Equal` and `Undetermined` pairs are not
reordered (stable sort preserves their physical order).
"""
function _partial_sort!(ops::Vector{QSym}, ne::Vector{NonEqualPair})
    n = length(ops)
    n < 2 && return ops
    # Insertion sort: stable, O(n²) worst case, but n is small in practice.
    for i in 2:n
        j = i
        while j > 1
            cmp = _site_compare(ops[j - 1], ops[j], ne)
            cmp === Greater || break
            ops[j - 1], ops[j] = ops[j], ops[j - 1]
            j -= 1
        end
    end
    return ops
end
```

- [ ] **Step 3: Add include in `src/SecondQuantizedAlgebra.jl`**

After the line `include("operators/operators.jl")` (line 25), add:
```julia
include("algebra/canonical_sort.jl")
```
The include order matters: `canonical_sort.jl` needs operator hooks defined first.

- [ ] **Step 4: Verify module loads**

Run:
```bash
julia --project=. -e 'using SecondQuantizedAlgebra'
```
Expected: loads cleanly.

- [ ] **Step 5: Pause**

### Task 3.2: Unit-test the partial sort

**Files:**
- Create: `test/canonical_sort_test.jl`

- [ ] **Step 1: Write the test**

```julia
using Test
using SecondQuantizedAlgebra
using SecondQuantizedAlgebra: _partial_sort!, _EMPTY_NE, Index, NO_INDEX

@testset "Partial sort: distinct sites" begin
    h = FockSpace(:c) ⊗ NLevelSpace(:a, 2)
    a = Destroy(h, :a, 1)
    σ = Transition(h, :σ, 1, 2)
    ops = QSym[σ, a]
    _partial_sort!(ops, _EMPTY_NE)
    # Destroy < Transition by type order
    @test ops[1] isa Destroy
    @test ops[2] isa Transition
end

@testset "Partial sort: same-site preserved" begin
    h = NLevelSpace(:a, 3)
    σ12 = Transition(h, :σ, 1, 2)
    σ21 = Transition(h, :σ, 2, 1)
    ops = QSym[σ12, σ21]
    pre = copy(ops)
    _partial_sort!(ops, _EMPTY_NE)
    @test ops == pre   # same site never reordered
end

@testset "Partial sort: undetermined preserved" begin
    ha = NLevelSpace(:a, 2)
    i = Index(ha, :i, 5, ha)
    j = Index(ha, :j, 5, ha)
    σ_i = Transition(ha, :σ, 1, 2, i)
    σ_j = Transition(ha, :σ, 1, 2, j)
    ops = QSym[σ_j, σ_i]   # physical order: j first
    _partial_sort!(ops, _EMPTY_NE)
    @test ops[1] === σ_j   # preserved because (i, j) relationship is Undetermined
    @test ops[2] === σ_i
end

@testset "Partial sort: ne resolves undetermined to distinct" begin
    ha = NLevelSpace(:a, 2)
    i = Index(ha, :i, 5, ha)
    j = Index(ha, :j, 5, ha)
    σ_i = Transition(ha, :σ, 1, 2, i)
    σ_j = Transition(ha, :σ, 1, 2, j)
    ops = QSym[σ_j, σ_i]
    _partial_sort!(ops, [(i, j)])
    # Now distinct: sort by index, i < j alphabetically
    @test ops[1] === σ_i
    @test ops[2] === σ_j
end
```

- [ ] **Step 2: Run the test**

Run:
```bash
julia --project=. test/canonical_sort_test.jl
```
Expected: all four testsets pass.

- [ ] **Step 3: Pause**

---

## Migration Step 4: Scaffold `src/algebra/`

### Task 4.1: Create empty stub files

**Files:**
- Create: `src/algebra/commute.jl`, `reduce.jl`, `expand_completeness.jl`, `substitute.jl`, `canonicalize.jl`, `pipelines.jl`, `weyl.jl`

- [ ] **Step 1: Create empty stubs**

Run:
```bash
for f in commute reduce expand_completeness substitute canonicalize pipelines weyl; do
    echo "# $f.jl — see docs/superpowers/specs/2026-05-15-src-pipeline-refactor-design.md" \
        > src/algebra/$f.jl
done
```

- [ ] **Step 2: Wire stubs into the module**

In `src/SecondQuantizedAlgebra.jl`, after the `include("algebra/canonical_sort.jl")` line, add:
```julia
include("algebra/canonicalize.jl")
include("algebra/commute.jl")
include("algebra/reduce.jl")
include("algebra/expand_completeness.jl")
include("algebra/substitute.jl")
include("algebra/pipelines.jl")
include("algebra/weyl.jl")
```
The order matters: passes before pipelines.

- [ ] **Step 3: Verify module loads with empty stubs**

Run:
```bash
julia --project=. -e 'using SecondQuantizedAlgebra'
```
Expected: loads cleanly (stubs contain only comments).

- [ ] **Step 4: Pause**

---

## Migration Step 5: Implement Passes

### Task 5.1: `_canonicalize_to_dict!` (terminal sink)

**Files:**
- Modify: `src/algebra/canonicalize.jl`

- [ ] **Step 1: Write the test**

Append to `test/canonical_sort_test.jl` (or create `test/passes_test.jl`):
```julia
using Test
using SecondQuantizedAlgebra
using SecondQuantizedAlgebra: _canonicalize_to_dict!, QTermDict, _CNUM_ONE, _CNUM_ZERO, _EMPTY_NE

@testset "canonicalize_to_dict! basic insert" begin
    h = FockSpace(:f)
    a = Destroy(h, :a)
    out = QTermDict()
    _canonicalize_to_dict!(out, QSym[a], _CNUM_ONE, _EMPTY_NE)
    @test length(out) == 1
end

@testset "canonicalize_to_dict! like-term collection" begin
    h = FockSpace(:f)
    a = Destroy(h, :a)
    out = QTermDict()
    _canonicalize_to_dict!(out, QSym[a], _CNUM_ONE, _EMPTY_NE)
    _canonicalize_to_dict!(out, QSym[a], _CNUM_ONE, _EMPTY_NE)
    @test length(out) == 1
    @test first(values(out)) == 2 + 0im
end

@testset "canonicalize_to_dict! zero coeff dropped" begin
    h = FockSpace(:f)
    a = Destroy(h, :a)
    out = QTermDict()
    _canonicalize_to_dict!(out, QSym[a], _CNUM_ZERO, _EMPTY_NE)
    @test isempty(out)
end
```

- [ ] **Step 2: Run test to verify it fails (function not defined)**

Run:
```bash
julia --project=. test/passes_test.jl 2>&1 | head -10
```
Expected: error about `_canonicalize_to_dict!` not defined.

- [ ] **Step 3: Implement `_canonicalize_to_dict!`**

In `src/algebra/canonicalize.jl`:
```julia
"""
    _canonicalize_to_dict!(out::QTermDict, ops::Vector{QSym}, c::CNum, ne::Vector{NonEqualPair})

Terminal sink for the pipeline. Constructs `QTerm(ops, ne)`, looks up in `out`,
sums `c`, drops zero-coefficient entries. **Takes ownership of `ops`** — the
caller must not mutate after this call.

Aliasing: the input `ops` is moved into the new `QTerm`. Pre-condition: `ops`
is in canonical form (post-passes).
"""
@inline function _canonicalize_to_dict!(
        out::QTermDict, ops::Vector{QSym}, c::CNum, ne::Vector{NonEqualPair},
    )
    _iszero_cnum(c) && return out
    return _addto_key!(out, QTerm(ops, _canonical_ne(ne)), c)
end
```

- [ ] **Step 4: Run the test**

Run:
```bash
julia --project=. test/passes_test.jl
```
Expected: pass.

- [ ] **Step 5: Pause**

### Task 5.2: `_reduce_ops` pass

**Files:**
- Modify: `src/algebra/reduce.jl`
- Modify: `test/passes_test.jl`

- [ ] **Step 1: Write tests**

Append to `test/passes_test.jl`:
```julia
using SecondQuantizedAlgebra: _reduce_ops

@testset "_reduce_ops: Transition composition" begin
    h = NLevelSpace(:a, 3)
    σ12 = Transition(h, :σ, 1, 2)
    σ23 = Transition(h, :σ, 2, 3)
    emitted = Tuple{Vector{QSym}, CNum}[]
    _reduce_ops(QSym[σ12, σ23], _CNUM_ONE) do o, c
        push!(emitted, (copy(o), c))
    end
    @test length(emitted) == 1
    @test length(emitted[1][1]) == 1
    @test emitted[1][1][1].i == 1 && emitted[1][1][1].j == 3
end

@testset "_reduce_ops: zero from incompatible composition" begin
    h = NLevelSpace(:a, 3)
    σ12 = Transition(h, :σ, 1, 2)
    σ31 = Transition(h, :σ, 3, 1)
    emitted = Tuple{Vector{QSym}, CNum}[]
    _reduce_ops(QSym[σ12, σ31], _CNUM_ONE) do o, c
        push!(emitted, (copy(o), c))
    end
    @test isempty(emitted)   # zero coefficient → nothing emitted
end

@testset "_reduce_ops: no-op input passes through" begin
    h = FockSpace(:f)
    a = Destroy(h, :a)
    ad = adjoint(a)
    emitted = Tuple{Vector{QSym}, CNum}[]
    _reduce_ops(QSym[ad, a], _CNUM_ONE) do o, c
        push!(emitted, (copy(o), c))
    end
    @test length(emitted) == 1
    @test emitted[1][1] == QSym[ad, a]
end
```

- [ ] **Step 2: Run tests, verify they fail**

Run:
```bash
julia --project=. test/passes_test.jl
```
Expected: errors for `_reduce_ops` not defined.

- [ ] **Step 3: Implement `_reduce_ops`**

In `src/algebra/reduce.jl`:
```julia
"""
    _reduce_ops(ops::Vector{QSym}, c::CNum, sink::F) where {F}

Fold adjacent provably-same-site pairs via `_reduce_pair`. Single output per
input. Aliasing: mutates `ops` in place; sink receives the same Vector.

Pre: `ops` is in canonical partial-sort order. Post: no adjacent provably-
same-site pair `(a, b)` for which `_reduce_pair(a, b)` returns non-`nothing`.
"""
@inline function _reduce_ops(ops::Vector{QSym}, c::CNum, sink::F) where {F}
    n = length(ops)
    n < 2 && (sink(ops, c); return)
    i = 1
    while i < length(ops)
        a, b = ops[i], ops[i + 1]
        res = _reduce_pair(a, b)
        if res === nothing
            i += 1
        elseif res isa CNum
            # Pair contracts to a scalar (e.g. σⁱʲσʲᵏ when j≠k yielding 0, or Pauli δ).
            _iszero_cnum(res) && return   # zero — emit nothing
            c = _mul_cnum(c, res)
            deleteat!(ops, i:(i + 1))
            i = max(i - 1, 1)
        else
            new_op, factor = res
            c = _mul_cnum(c, factor)
            ops[i] = new_op
            deleteat!(ops, i + 1)
            i = max(i - 1, 1)
        end
    end
    sink(ops, c)
    return
end
```

- [ ] **Step 4: Run tests**

Run:
```bash
julia --project=. test/passes_test.jl
```
Expected: all `_reduce_ops` tests pass.

- [ ] **Step 5: Pause**

### Task 5.3: `_commute_ops` pass

**Files:**
- Modify: `src/algebra/commute.jl`
- Modify: `test/passes_test.jl`

- [ ] **Step 1: Write tests**

Append to `test/passes_test.jl`:
```julia
using SecondQuantizedAlgebra: _commute_ops

@testset "_commute_ops: Fock aa† → a†a + 1" begin
    h = FockSpace(:f)
    a = Destroy(h, :a)
    ad = adjoint(a)
    emitted = Tuple{Vector{QSym}, CNum}[]
    _commute_ops(QSym[a, ad], _CNUM_ONE) do o, c
        push!(emitted, (copy(o), c))
    end
    # Two terms: a†a (the swap) and 1 (the residual, ops empty)
    @test length(emitted) == 2
    sort!(emitted, by = e -> length(e[1]))
    @test isempty(emitted[1][1]) && emitted[1][2] == _CNUM_ONE   # the "1" branch
    @test emitted[2][1] == QSym[ad, a]                            # the swapped branch
end

@testset "_commute_ops: no-op on already-ordered pair" begin
    h = FockSpace(:f)
    a = Destroy(h, :a)
    ad = adjoint(a)
    emitted = Tuple{Vector{QSym}, CNum}[]
    _commute_ops(QSym[ad, a], _CNUM_ONE) do o, c
        push!(emitted, (copy(o), c))
    end
    @test length(emitted) == 1
    @test emitted[1][1] == QSym[ad, a]
end
```

- [ ] **Step 2: Verify they fail**

Run:
```bash
julia --project=. test/passes_test.jl
```
Expected: error.

- [ ] **Step 3: Implement `_commute_ops`**

In `src/algebra/commute.jl`:
```julia
"""
    _commute_ops(ops::Vector{QSym}, c::CNum, sink::F) where {F}

Apply commutation swaps for adjacent provably-same-site pairs that need them.
Undetermined-site and distinct-site adjacents are untouched. May fork: emits
the swapped branch through the sink, plus a residual branch with the
contracted form.

Aliasing: the contracted branch mutates `ops` in place and continues. The
swapped branch receives a `copy(ops)`. The terminal sink takes ownership.
"""
@inline function _commute_ops(ops::Vector{QSym}, c::CNum, sink::F) where {F}
    _commute_ops_at(ops, c, 1, sink)
    return
end

function _commute_ops_at(ops::Vector{QSym}, c::CNum, start::Int, sink::F) where {F}
    n = length(ops)
    i = start
    while i < length(ops)
        a, b = ops[i], ops[i + 1]
        # Only act on pairs known to be same-site (Equal) where _can_commute is false.
        cmp = _site_compare(a, b, _EMPTY_NE)   # ne not threaded through passes
        if cmp === Equal && !_can_commute(a, b)
            sw, contracted_form_a, residual = _commute_pair(a, b)
            # Forked branch: ops with (a, b) → (sw_b, sw_a) — i.e. the swap.
            swapped = copy(ops)
            swapped[i]     = sw           # = b (per _commute_pair convention)
            swapped[i + 1] = contracted_form_a   # = a
            # Recursively process the swapped branch from position i (now in canonical order).
            _commute_ops_at(swapped, c, max(i, 1), sink)
            # Contracted branch: pair removed, multiply coefficient by residual.
            new_c = _mul_cnum(c, residual)
            deleteat!(ops, i:(i + 1))
            _iszero_cnum(new_c) && return
            c = new_c
            i = max(i - 1, 1)
        else
            i += 1
        end
    end
    sink(ops, c)
    return
end
```

- [ ] **Step 4: Run tests**

Run:
```bash
julia --project=. test/passes_test.jl
```
Expected: all `_commute_ops` tests pass.

- [ ] **Step 5: Pause**

### Task 5.4: `_expand_gs_ops` pass

**Files:**
- Modify: `src/algebra/expand_completeness.jl`
- Modify: `test/passes_test.jl`

- [ ] **Step 1: Write tests**

Append to `test/passes_test.jl`:
```julia
using SecondQuantizedAlgebra: _expand_gs_ops

@testset "_expand_gs_ops: σ¹¹ → 1 - σ²²" begin
    h = NLevelSpace(:a, 2)
    σ11 = Transition(h, :σ, 1, 1)
    emitted = Tuple{Vector{QSym}, CNum}[]
    _expand_gs_ops(QSym[σ11], _CNUM_ONE) do o, c
        push!(emitted, (copy(o), c))
    end
    # Expect two emissions: (1, +c) and (σ²², -c)
    @test length(emitted) == 2
    sort!(emitted, by = e -> length(e[1]))
    @test isempty(emitted[1][1]) && emitted[1][2] == _CNUM_ONE
    @test length(emitted[2][1]) == 1 && emitted[2][2] == -_CNUM_ONE
    @test emitted[2][1][1].i == 2 && emitted[2][1][1].j == 2
end

@testset "_expand_gs_ops: passthrough when no σᵍᵍ" begin
    h = FockSpace(:f)
    a = Destroy(h, :a)
    emitted = Tuple{Vector{QSym}, CNum}[]
    _expand_gs_ops(QSym[a], _CNUM_ONE) do o, c
        push!(emitted, (copy(o), c))
    end
    @test length(emitted) == 1
    @test emitted[1][1] == QSym[a]
end
```

- [ ] **Step 2: Verify failure**

Run:
```bash
julia --project=. test/passes_test.jl
```
Expected: error for `_expand_gs_ops`.

- [ ] **Step 3: Implement `_expand_gs_ops` and `expand_completeness`**

In `src/algebra/expand_completeness.jl`:
```julia
"""
    _expand_gs_ops(ops::Vector{QSym}, c::CNum, sink::F) where {F}

Apply `σᵍᵍ → 1 - Σ_{k≠g} σᵏᵏ` to every `σᵍᵍ` in `ops`. May fork by `n_levels`
per `σᵍᵍ`. Recurses to handle multiple `σᵍᵍ` atoms in one term.

Aliasing: each branch receives its own `copy(ops)`; original is not mutated.
"""
function _expand_gs_ops(ops::Vector{QSym}, c::CNum, sink::F) where {F}
    # Find first σᵍᵍ.
    idx = 0
    for k in eachindex(ops)
        _ground_state_expand(ops[k]) === nothing || (idx = k; break)
    end
    idx == 0 && (sink(ops, c); return)

    op = ops[idx]
    (g, n_levels, site) = _ground_state_expand(op)
    # Emit identity branch: ops with op removed, coefficient unchanged.
    id_ops = QSym[ops[k] for k in eachindex(ops) if k != idx]
    _expand_gs_ops(id_ops, c, sink)
    # Emit -σᵏᵏ branch for each k ≠ g.
    neg_c = _neg_cnum(c)
    for k in 1:n_levels
        k == g && continue
        new_ops = copy(ops)
        new_ops[idx] = Transition(op.name, k, k, op.space_index, op.index, g, n_levels)
        _expand_gs_ops(new_ops, neg_c, sink)
    end
    return
end

"""
    expand_completeness(q::QAdd) -> QAdd

Apply the ground-state completeness identity `σᵍᵍ = 1 - Σ_{k≠g} σᵏᵏ` to every
`σᵍᵍ` in `q`, then re-canonicalize.
"""
function expand_completeness(q::QAdd)
    out = QTermDict()
    for (t, c) in q
        _expand_gs_ops(copy(t.ops), c) do ops1, c1
            _stream!(out, ops1, c1, t.ne)
        end
    end
    return QAdd(out, copy(q.indices))
end
```

(Note: `_stream!` is defined in pipelines.jl in Task 6.1.)

- [ ] **Step 4: Pause** — full test of expand_completeness deferred to Task 6.1 when `_stream!` exists.

### Task 5.5: `_substitute_ops` pass

**Files:**
- Modify: `src/algebra/substitute.jl`
- Modify: `test/passes_test.jl`

- [ ] **Step 1: Implement `_substitute_ops`**

In `src/algebra/substitute.jl`:
```julia
"""
    _substitute_ops(ops::Vector{QSym}, c::CNum, d, sink::F) where {F}

Walk `ops` applying substitutions from `d`. `d` may contain:
- QSym → CNum or Number (operator → scalar; folds into coefficient)
- QSym → QAdd (operator → expression; splices the QAdd's single term in)
- Symbolic → Symbolic (applied to coefficient via SymbolicUtils.substitute)

For multi-term QAdd substitutions, forks per term.

Aliasing: each fork receives its own copy.
"""
function _substitute_ops(ops::Vector{QSym}, c::CNum, d, sink::F) where {F}
    # Walk ops, find first op present in d, substitute.
    for i in eachindex(ops)
        op = ops[i]
        haskey(d, op) || continue
        val = d[op]
        if val isa Number
            new_c = _mul_cnum(c, _to_cnum(val))
            _iszero_cnum(new_c) && return
            new_ops = QSym[ops[k] for k in eachindex(ops) if k != i]
            return _substitute_ops(new_ops, new_c, d, sink)
        elseif val isa QAdd
            # Splice each term of val into position i.
            for (vt, vc) in val
                spliced_ops = copy(ops)
                deleteat!(spliced_ops, i)
                for (k, vop) in enumerate(vt.ops)
                    insert!(spliced_ops, i + k - 1, vop)
                end
                new_c = _mul_cnum(c, vc)
                _substitute_ops(spliced_ops, new_c, d, sink)
            end
            return
        end
    end
    # No operator substitutions matched. Apply symbolic substitution to c.
    new_c = _substitute_cnum(c, d)
    _iszero_cnum(new_c) && return
    sink(ops, new_c)
    return
end

# Helper: apply symbolic substitution to a CNum coefficient.
# (Stub — port logic from current src/arithmetics/substitute.jl line 37.)
function _substitute_cnum(c::CNum, d)
    rep, imp = real(c), imag(c)
    new_rep = SymbolicUtils.substitute(Symbolics.value(rep), d)
    new_imp = SymbolicUtils.substitute(Symbolics.value(imp), d)
    return Num(new_rep) + im * Num(new_imp)
end
```

- [ ] **Step 2: Write a basic test**

Append to `test/passes_test.jl`:
```julia
using SecondQuantizedAlgebra: _substitute_ops

@testset "_substitute_ops: operator → scalar" begin
    h = FockSpace(:f)
    a = Destroy(h, :a)
    emitted = Tuple{Vector{QSym}, CNum}[]
    _substitute_ops(QSym[a, adjoint(a)], _CNUM_ONE, Dict(a => 2)) do o, c
        push!(emitted, (copy(o), c))
    end
    @test length(emitted) == 1
    @test emitted[1][2] == 2 + 0im
    @test emitted[1][1] == QSym[adjoint(a)]
end
```

- [ ] **Step 3: Run test**

Run:
```bash
julia --project=. test/passes_test.jl
```
Expected: pass.

- [ ] **Step 4: Pause**

---

## Migration Step 6: Wire Entry Points

### Task 6.1: `_stream!` and `_canonicalize!` primitive

**Files:**
- Modify: `src/algebra/pipelines.jl`

- [ ] **Step 1: Implement `_stream!` and `_canonicalize!`**

In `src/algebra/pipelines.jl`:
```julia
"""
    _stream!(out::QTermDict, ops::Vector{QSym}, c::CNum, ne::Vector{NonEqualPair})

The default canonicalization pipeline: commute → reduce → insert.
"""
@inline function _stream!(
        out::QTermDict, ops::Vector{QSym}, c::CNum, ne::Vector{NonEqualPair},
    )
    _commute_ops(ops, c) do ops1, c1
        _reduce_ops(ops1, c1) do ops2, c2
            _canonicalize_to_dict!(out, ops2, c2, ne)
        end
    end
    return nothing
end

"""
    _canonicalize!(out::QTermDict, ops::Vector{QSym}, c::CNum, ne::Vector{NonEqualPair})

Re-establish canonical form: sort, then run the pipeline.
"""
@inline function _canonicalize!(
        out::QTermDict, ops::Vector{QSym}, c::CNum, ne::Vector{NonEqualPair},
    )
    _partial_sort!(ops, ne)
    _stream!(out, ops, c, ne)
    return nothing
end
```

- [ ] **Step 2: Add a smoke test**

Append to `test/passes_test.jl`:
```julia
using SecondQuantizedAlgebra: _stream!, _canonicalize!

@testset "_stream!: idempotent on canonical input" begin
    h = FockSpace(:f)
    a = Destroy(h, :a)
    out = QTermDict()
    _stream!(out, QSym[adjoint(a), a], _CNUM_ONE, _EMPTY_NE)
    @test length(out) == 1   # a†a stays as a†a (canonical form)
end

@testset "_canonicalize!: aa† → a†a + 1" begin
    h = FockSpace(:f)
    a = Destroy(h, :a)
    out = QTermDict()
    _canonicalize!(out, QSym[a, adjoint(a)], _CNUM_ONE, _EMPTY_NE)
    @test length(out) == 2   # a†a term and identity term
end
```

- [ ] **Step 3: Run**

Run:
```bash
julia --project=. test/passes_test.jl
```
Expected: pass.

- [ ] **Step 4: Pause**

### Task 6.2: `_emit_product!` and the four `Base.:*` overloads

**Files:**
- Modify: `src/algebra/pipelines.jl`

- [ ] **Step 1: Implement `_emit_product!` and `_accumulate_with_diag!` (stub)**

Append to `src/algebra/pipelines.jl`:
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

# Stub — full implementation in Task 6.3.
function _accumulate_with_diag!(out, ops, c, sum_indices, ne)
    _canonicalize!(out, copy(ops), c, ne)
    return nothing
end
```

- [ ] **Step 2: Add `Base.:*` overloads in `algebra/pipelines.jl`**

The new methods live in their permanent home — `src/algebra/pipelines.jl`. Append:

```julia
function Base.:*(a::QSym, b::QSym)
    out = QTermDict()
    _emit_product!(out,
        QSym[a], _CNUM_ONE, _EMPTY_NE,
        QSym[b], _CNUM_ONE, _EMPTY_NE,
        Index[], false)
    return QAdd(out, Index[])
end

function Base.:*(a::QAdd, b::QSym)
    out = QTermDict()
    needs = !isempty(a.indices)
    for (ta, ca) in a
        _emit_product!(out, ta.ops, ca, ta.ne,
                            QSym[b], _CNUM_ONE, _EMPTY_NE,
                            a.indices, needs)
    end
    return QAdd(out, copy(a.indices))
end

function Base.:*(a::QSym, b::QAdd)
    out = QTermDict()
    needs = !isempty(b.indices)
    for (tb, cb) in b
        _emit_product!(out, QSym[a], _CNUM_ONE, _EMPTY_NE,
                            tb.ops, cb, tb.ne,
                            b.indices, needs)
    end
    return QAdd(out, copy(b.indices))
end

function Base.:*(a::QAdd, b::QAdd)
    if !isempty(a.indices) && !isempty(b.indices)
        for idx in a.indices
            idx in b.indices && throw(
                ArgumentError(
                    "Summation index $(idx.name) appears in both factors."
                )
            )
        end
    end
    out = QTermDict()
    sum_indices = _merge_unique(a.indices, b.indices)
    needs = !isempty(sum_indices)
    for (ta, ca) in a, (tb, cb) in b
        _emit_product!(out, ta.ops, ca, ta.ne,
                            tb.ops, cb, tb.ne,
                            sum_indices, needs)
    end
    return QAdd(out, sum_indices)
end
```

- [ ] **Step 3: Delete the old `Base.:*` definitions from `src/arithmetics/qadd_arithmetic.jl`**

Julia would warn about method redefinition because `algebra/pipelines.jl` is included after `arithmetics/qadd_arithmetic.jl`. To avoid that, delete the old definitions now:

Open [src/arithmetics/qadd_arithmetic.jl](src/arithmetics/qadd_arithmetic.jl) and remove lines that define:
- `Base.:*(a::QSym, b::QSym)` (around line 5)
- `Base.:*(a::QAdd, b::QSym)` (around line 150)
- `Base.:*(a::QSym, b::QAdd)` (around line 162)
- `Base.:*(a::QAdd, b::QAdd)` (around line 174)

Keep the scalar overloads (`Base.:*(a::QSym, b::Number)` etc.) and `Base.:^` for now — they survive into step 11 untouched, OR move them to `pipelines.jl` too if convenient. The scalar overloads do not touch the algebra pipeline.

- [ ] **Step 4: Run the test suite to see how many pass**

Run:
```bash
make test 2>&1 | tail -40
```
Expected: most non-indexed tests pass; indexed tests may fail due to the stub `_accumulate_with_diag!`.

- [ ] **Step 5: Pause**

### Task 6.3: Full `_accumulate_with_diag!` and helpers

**Files:**
- Modify: `src/algebra/pipelines.jl`

- [ ] **Step 1: Port helpers from `qadd_arithmetic.jl`**

Look at lines 57-66 (`_distinct_op_indices`), 77-97 (`_emit_diagonal!`), and 192-200 of `src/expressions/index.jl` (`_depends_on_index_term`).

Append to `src/algebra/pipelines.jl`:
```julia
function _distinct_op_indices(ops::Vector{QSym})
    out = Index[]
    for op in ops
        idx = op.index
        has_index(idx) || continue
        idx in out && continue
        push!(out, idx)
    end
    return out
end

function _depends_on_index_ops(c::CNum, ops::Vector{QSym}, idx::Index)
    # Operator dependence: any op carries the index.
    for op in ops
        op.index == idx && return true
    end
    # Coefficient dependence: delegate to current index.jl helper.
    return _depends_on_index_term(c, ops, idx)
end

function _emit_diagonal!(
        out::QTermDict, ops::Vector{QSym}, c::CNum,
        ne::Vector{NonEqualPair}, α::Index, β::Index,
    )
    new_ne = _merge_ne_pair(ne, α, β)
    _canonicalize!(out, ops, c, new_ne)
    return nothing
end
```

- [ ] **Step 2: Replace the `_accumulate_with_diag!` stub with the full implementation**

Replace the stub (added in Task 6.2 Step 1) with:
```julia
function _accumulate_with_diag!(
        out::QTermDict, ops::Vector{QSym}, c::CNum,
        sum_indices::Vector{Index}, ne::Vector{NonEqualPair},
    )
    _canonicalize!(out, copy(ops), c, ne)   # off-diagonal

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

- [ ] **Step 3: Run indexed tests**

Run:
```bash
julia --project=. -e 'using Pkg; Pkg.test(SecondQuantizedAlgebra; test_args=["indexing_test"])'
```

If `test_args` filtering isn't supported, run:
```bash
julia --project=. test/indexing_test.jl
```
Expected: most tests pass; some may still fail due to test assertions touching `phys_ops` or eager-GS.

- [ ] **Step 4: Pause**

### Task 6.4: Rewrite `Base.adjoint` to use `_canonicalize!`

**Files:**
- Modify: [src/expressions/qadd.jl:71-81](src/expressions/qadd.jl#L71-L81)

- [ ] **Step 1: Replace `Base.adjoint`**

Replace the current (already cleaned in Task 1.3) body with:
```julia
function Base.adjoint(q::QAdd)
    out = QTermDict()
    for (t, c) in q
        rev = QSym[adjoint(o) for o in Iterators.reverse(t.ops)]
        _canonicalize!(out, rev, conj(c), t.ne)
    end
    return QAdd(out, copy(q.indices))
end
```

- [ ] **Step 2: Run adjoint-related tests**

Run:
```bash
julia --project=. test/qadd_test.jl
```
Expected: pass.

- [ ] **Step 3: Pause**

### Task 6.5: Wire `change_index` to use `_canonicalize!`

**Files:**
- Modify: [src/expressions/index.jl](src/expressions/index.jl) lines around 224-260 (`_diagonal_split!`) and 128, 241-244, 257 (per audit Item A.5).

- [ ] **Step 1: Locate `change_index(q::QAdd, ...)` and `_diagonal_split!`**

Run:
```bash
grep -n "function change_index\|function _diagonal_split" src/expressions/index.jl
```

- [ ] **Step 2: Rewrite `change_index(q::QAdd, from::Index, to::Index)`**

Replace the body so it transforms `t.ops` and `t.ne`, then routes through `_canonicalize!`:
```julia
function change_index(q::QAdd, from::Index, to::Index)
    out = QTermDict()
    needs = !isempty(q.indices)
    sum_indices = q.indices
    for (t, c) in q
        new_ops = QSym[change_index(o, from, to) for o in t.ops]
        new_c   = change_index(c, from, to)
        new_ne  = _substitute_ne(t.ne, from, to)
        if needs
            _accumulate_with_diag!(out, new_ops, new_c, sum_indices, new_ne)
        else
            _canonicalize!(out, new_ops, new_c, new_ne)
        end
    end
    return QAdd(out, copy(q.indices))
end
```

(Note: `change_index(o, from, to)` for QSym leaves operators alone if their index doesn't match; this is the existing per-operator `change_index` method.)

- [ ] **Step 3: Delete `_diagonal_split!` if it became unused**

Check usage:
```bash
grep -n "_diagonal_split" src/
```
If only `change_index` referenced it, remove the function. Otherwise keep and update callers similarly.

- [ ] **Step 4: Run indexed tests**

Run:
```bash
julia --project=. test/indexing_test.jl
```
Expected: pass.

- [ ] **Step 5: Pause**

### Task 6.6: Replace `simplify`, `normal_order`, `substitute`, `commutator`

**Files:**
- Modify: [src/arithmetics/simplify.jl](src/arithmetics/simplify.jl)
- Modify: [src/arithmetics/normal_order.jl](src/arithmetics/normal_order.jl)
- Modify: [src/arithmetics/substitute.jl](src/arithmetics/substitute.jl)
- Modify: [src/arithmetics/commutator.jl](src/arithmetics/commutator.jl)

These files are being deleted in step 11. We replace their contents with the new entry points first; deletion is just removing dead files.

- [ ] **Step 1: Replace `simplify.jl`**

Replace the entire contents of `src/arithmetics/simplify.jl` with:
```julia
# This file is deleted in migration step 11; logic moved to algebra/pipelines.jl.
# Until then, simplify is defined here as an alias to normal_order.
simplify(q::QAdd) = normal_order(q)
```

- [ ] **Step 2: Replace `normal_order.jl`**

Read [src/arithmetics/normal_order.jl](src/arithmetics/normal_order.jl) lines 1-50 to find the existing `normal_order(q::QAdd)` and Weyl-conversion functions.

Replace `normal_order(q::QAdd)`:
```julia
function normal_order(q::QAdd)
    out = QTermDict()
    for (t, c) in q
        _stream!(out, copy(t.ops), c, t.ne)
    end
    return QAdd(out, copy(q.indices))
end
```

Remove `normal_order(expr, h::HilbertSpace)` overloads and `_apply_ground_state` (lines 241-249) — they are replaced by `expand_completeness`.

Keep `normal_to_symmetric` and `symmetric_to_normal` (Weyl conversion); these move to `algebra/weyl.jl` in step 10.

- [ ] **Step 3: Replace `substitute.jl`**

Replace the entire contents with:
```julia
function substitute(q::QAdd, d)
    out = QTermDict()
    for (t, c) in q
        _substitute_ops(copy(t.ops), c, d) do ops1, c1
            _stream!(out, ops1, c1, t.ne)
        end
    end
    return QAdd(out, copy(q.indices))
end

# Also accept QSym inputs: promote to QAdd.
substitute(s::QSym, d) = substitute(_single_qadd(_CNUM_ONE, QSym[s]), d)
```

- [ ] **Step 4: Replace `commutator.jl`**

Replace with:
```julia
commutator(a, b) = a * b - b * a
anticommutator(a, b) = a * b + b * a
```

- [ ] **Step 5: Run the full test suite**

Run:
```bash
make test 2>&1 | tail -50
```
Expected: most pass. Failures should now only be tests that assert eager GS expansion or access `phys_ops` directly — those get fixed in Step 7.

- [ ] **Step 6: Pause**

### Task 6.7: Export `expand_completeness`

**Files:**
- Modify: [src/SecondQuantizedAlgebra.jl:89-106](src/SecondQuantizedAlgebra.jl#L89-L106)

- [ ] **Step 1: Add `expand_completeness` to the export list**

Find the `export` block (lines 89-106). Add `expand_completeness` after `simplify`:
```julia
    simplify, expand, expand_completeness, commutator, anticommutator,
```

- [ ] **Step 2: Verify the module exports**

Run:
```bash
julia --project=. -e 'using SecondQuantizedAlgebra; @assert isdefined(SecondQuantizedAlgebra, :expand_completeness)'
```
Expected: no error.

- [ ] **Step 3: Pause**

---

## Migration Step 7: Update Existing Tests

### Task 7.1: Update `test/nlevel_test.jl`

**Files:**
- Modify: [test/nlevel_test.jl:50,52,63-64](test/nlevel_test.jl)

- [ ] **Step 1: Update line 50 (eager GS expansion)**

Find:
```julia
@test isequal(simplify(normal_order(σ * σ')), simplify(1 - σee))
```

Replace with:
```julia
@test isequal(simplify(expand_completeness(σ * σ')), simplify(1 - σee))
```

(Or assert against atomic form: `@test isequal(σ * σ', σgg)` — pick the one that matches user intent for this test.)

- [ ] **Step 2: Update line 52 (`simplify(σgg, h)`)**

Find:
```julia
@test isequal(simplify(σgg, h), simplify(1 - σee, h))
```

Replace with:
```julia
@test isequal(expand_completeness(σgg), simplify(1 - σee))
```

- [ ] **Step 3: Update lines 63-64**

Read the lines and update along the same pattern.

- [ ] **Step 4: Run nlevel tests**

Run:
```bash
julia --project=. test/nlevel_test.jl
```
Expected: pass.

- [ ] **Step 5: Pause**

### Task 7.2: Update `test/normal_order_test.jl`

**Files:**
- Modify: [test/normal_order_test.jl:85-99](test/normal_order_test.jl)

- [ ] **Step 1: Update lines 85-89 (normal_order with h argument)**

Find:
```julia
@test isequal(normal_order(σ11, h), ...)
```

Replace with:
```julia
@test isequal(expand_completeness(σ11), ...)
```

- [ ] **Step 2: Update lines 97-99 (same pattern for product space)**

Same transformation.

- [ ] **Step 3: Run**

Run:
```bash
julia --project=. test/normal_order_test.jl
```
Expected: pass.

- [ ] **Step 4: Pause**

### Task 7.3: Update `test/integration_test.jl` and audit for `phys_ops`

**Files:**
- Modify: [test/integration_test.jl](test/integration_test.jl)

- [ ] **Step 1: Grep for `phys_ops` references**

Run:
```bash
grep -rn "phys_ops" test/
```

- [ ] **Step 2: For each match, decide on the fix**

Patterns to expect:
- Direct `term.phys_ops` access → delete the line if it was asserting a structural detail, or replace with `term.ops` if it was asserting the operator string.
- `QTerm(ops, ne, phys_ops)` constructor → drop the third argument.

Edit each match.

- [ ] **Step 3: Run the full test suite**

Run:
```bash
make test 2>&1 | tail -30
```
Expected: all tests green.

- [ ] **Step 4: Pause**

---

## Migration Step 8: Canonical-Form Invariant Test

### Task 8.1: Create `test/canonical_invariant_test.jl`

**Files:**
- Create: `test/canonical_invariant_test.jl`
- Modify: [test/runtests.jl](test/runtests.jl)

- [ ] **Step 1: Write the invariant test file**

```julia
using Test
using SecondQuantizedAlgebra
using SecondQuantizedAlgebra: _site_compare, _can_commute, _reduce_pair,
    SiteCmp, Less, Equal, Undetermined, Greater, _EMPTY_NE

# Returns true if `ops` has no adjacent provably-same-site pair with an
# unresolved commute or reduce.
function _is_canonical(t)
    ops = t.ops
    ne = t.ne
    for i in 1:(length(ops) - 1)
        cmp = _site_compare(ops[i], ops[i + 1], ne)
        if cmp === Greater
            return false   # distinct sites in wrong order
        end
        if cmp === Equal
            if !_can_commute(ops[i], ops[i + 1])
                return false   # unresolved commutation
            end
            if _reduce_pair(ops[i], ops[i + 1]) !== nothing
                return false   # unresolved reduction
            end
        end
    end
    return true
end

@testset "Canonical-form invariant: `*`" begin
    h = FockSpace(:f) ⊗ NLevelSpace(:a, 3)
    a = Destroy(h, :a, 1)
    σ12 = Transition(h, :σ, 1, 2, 2)
    σ21 = Transition(h, :σ, 2, 1, 2)
    q = a * adjoint(a) * σ12 * σ21
    for (t, _) in q
        @test _is_canonical(t)
    end
end

@testset "Canonical-form invariant: normal_order" begin
    h = FockSpace(:f)
    a = Destroy(h, :a)
    q = a * adjoint(a) * a
    for (t, _) in normal_order(q)
        @test _is_canonical(t)
    end
end

@testset "Canonical-form invariant: commutator" begin
    h = FockSpace(:f)
    a = Destroy(h, :a)
    q = commutator(a, adjoint(a))
    for (t, _) in q
        @test _is_canonical(t)
    end
end

@testset "Canonical-form invariant: substitute" begin
    h = FockSpace(:f)
    a = Destroy(h, :a)
    q = adjoint(a) * a
    sub = substitute(q, Dict(a => 0.5))
    for (t, _) in sub
        @test _is_canonical(t)
    end
end

@testset "Canonical-form invariant: expand_completeness" begin
    h = NLevelSpace(:a, 3)
    σ11 = Transition(h, :σ, 1, 1)
    exp = expand_completeness(σ11)
    for (t, _) in exp
        @test _is_canonical(t)
    end
end
```

- [ ] **Step 2: Add to test runner**

Open [test/runtests.jl](test/runtests.jl). Add to the list of test files:
```julia
"canonical_invariant_test"
```

- [ ] **Step 3: Run**

Run:
```bash
julia --project=. test/canonical_invariant_test.jl
```
Expected: all pass.

- [ ] **Step 4: Pause**

---

## Migration Step 9: Dissipator Regression Test (PR #99)

### Task 9.1: Create `test/dissipator_regression_test.jl`

**Files:**
- Create: `test/dissipator_regression_test.jl`
- Modify: [test/runtests.jl](test/runtests.jl)

- [ ] **Step 1: Write the regression test**

Port the MWE from PR #99's comment by ChristophHotter:

```julia
using Test
using SecondQuantizedAlgebra
using Symbolics: @variables

@testset "PR #99 dissipator regression" begin
    hc = FockSpace(:cavity)
    ha = NLevelSpace(:atom, 2)
    h = hc ⊗ ha

    a = Destroy(h, :a, 1)
    σ(α, β, i) = IndexedOperator(Transition(h, :σ, α, β, 2), i)

    @variables N Δ κ γ R ν

    i = Index(h, :i, N, ha)
    j = Index(h, :j, N, ha)
    k = Index(h, :k, N, ha)

    J = [a, σ(1, 2, i), σ(2, 1, i), σ(2, 2, i)]
    Jd = adjoint.(J)
    rates = [κ, γ, R, ν]

    function decay(op)
        return simplify(
            0.5 * rates[1] * (2 * Jd[1] * op * J[1] - Jd[1] * J[1] * op - op * Jd[1] * J[1]) +
            sum(Σ(0.5 * rates[x] * (2 * Jd[x] * op * J[x] - Jd[x] * J[x] * op - op * Jd[x] * J[x]), i)
                for x in 2:length(J))
        )
    end

    # Expected results from PR #99:
    @test isequal(decay(a' * a), simplify(-κ * (a' * a)))
    @test isequal(decay(σ(2, 2, j)), simplify(R - (R + γ) * σ(2, 2, j)))
    @test isequal(decay(a' * σ(1, 2, j)), simplify(-(γ + κ + R + ν) / 2 * (a' * σ(1, 2, j))))
    @test isequal(decay(σ(1, 2, j) * σ(2, 1, k)), simplify(-(γ + R + ν) * (σ(1, 2, j) * σ(2, 1, k))))
end
```

- [ ] **Step 2: Add to test runner**

Add `"dissipator_regression_test"` to the test file list in [test/runtests.jl](test/runtests.jl).

- [ ] **Step 3: Run**

Run:
```bash
julia --project=. test/dissipator_regression_test.jl
```
Expected: all pass. If a test fails, the rewrite did not fix the regression — investigate before continuing.

- [ ] **Step 4: Pause**

---

## Migration Step 10: Move Weyl Conversion

### Task 10.1: Relocate `normal_to_symmetric` and `symmetric_to_normal`

**Files:**
- Move: from [src/arithmetics/normal_order.jl](src/arithmetics/normal_order.jl) (lines 71-193 per audit Item F)
- To: `src/algebra/weyl.jl`

- [ ] **Step 1: Copy the Weyl-conversion functions to the new file**

Open `src/algebra/weyl.jl` (stubbed in Task 4.1). Replace its single comment line with the copied contents of `normal_to_symmetric`, `symmetric_to_normal`, and `_convert_term!`, `_convert_ordering` from the original `normal_order.jl`.

- [ ] **Step 2: Delete the moved functions from `normal_order.jl`**

Open `src/arithmetics/normal_order.jl`. Delete the same lines that were copied (the file will be deleted entirely in step 11; we just keep it clean for the interim).

- [ ] **Step 3: Run Weyl-related tests**

Run:
```bash
grep -l "normal_to_symmetric\|symmetric_to_normal" test/
# For each match:
julia --project=. test/<match>.jl
```
Expected: pass.

- [ ] **Step 4: Pause**

---

## Migration Step 11: Delete `src/arithmetics/`

### Task 11.1: Remove the directory and update includes

**Files:**
- Delete: `src/arithmetics/ordering.jl`, `qadd_arithmetic.jl`, `simplify.jl`, `substitute.jl`, `commutator.jl`, `normal_order.jl`
- Modify: [src/SecondQuantizedAlgebra.jl:28-38](src/SecondQuantizedAlgebra.jl#L28-L38)

- [ ] **Step 1: Verify all functions defined in `src/arithmetics/` are now defined elsewhere**

Run:
```bash
grep -rn "^function\|^@inline function" src/arithmetics/ | awk '{print $0}'
```
Each function listed must have a replacement in `src/algebra/` or be intentionally deleted. Cross-check.

- [ ] **Step 2: Delete the `src/arithmetics/` directory**

Run:
```bash
rm -r src/arithmetics
```

- [ ] **Step 3: Update `src/SecondQuantizedAlgebra.jl` include list**

Remove these lines from `src/SecondQuantizedAlgebra.jl`:
```julia
include("arithmetics/ordering.jl")
include("arithmetics/qadd_arithmetic.jl")
include("arithmetics/simplify.jl")
include("arithmetics/substitute.jl")
include("arithmetics/commutator.jl")
include("arithmetics/normal_order.jl")
```

The `include("algebra/...")` lines (added in Step 4) replace them in correct order.

- [ ] **Step 4: Run the full test suite**

Run:
```bash
make test 2>&1 | tail -20
```
Expected: all tests pass.

- [ ] **Step 5: Pause**

---

## Migration Step 12: Final Cleanup

### Task 12.1: Update precompile workload

**Files:**
- Modify: [src/precompile.jl](src/precompile.jl)

- [ ] **Step 1: Update signatures of any `precompile(...)` calls that reference removed functions**

Open `src/precompile.jl`. Find any `precompile(_site_sort!, ...)`, `precompile(_addto!, (..., Vector{QSym}))` calls — update signatures to match the new shapes (drop `phys_ops` argument tuple).

Replace with concrete-method precompiles for the new hot path:
```julia
@setup_workload begin
    h = FockSpace(:f)
    @compile_workload begin
        a = Destroy(h, :a)
        ad = adjoint(a)
        # Exercise the hot paths
        commutator(a, ad)
        normal_order(a * ad * a)
    end
end
```

- [ ] **Step 2: Run with precompilation enabled**

Run:
```bash
julia --project=. -e 'using SecondQuantizedAlgebra; @time using SecondQuantizedAlgebra' 2>&1
```
Expected: first call slow, second call fast. (Or test TTFX with `julia-ttfx` skill if needed.)

- [ ] **Step 3: Pause**

### Task 12.2: Add `src/algebra/README.md` (debugging idiom)

**Files:**
- Create: `src/algebra/README.md`

- [ ] **Step 1: Write the README**

Create `src/algebra/README.md`:
```markdown
# `src/algebra/` — the operator algebra pipeline

This directory implements `*`, `simplify`, `normal_order`, `commutator`,
`substitute`, `expand_completeness`, and `Base.adjoint` as a streaming per-term
pipeline of small named passes.

## Pipeline

```
Base.:* → _emit_product! → _canonicalize! → _partial_sort! + _stream!
                                                    │
                                                    ↓
                                          _commute_ops → _reduce_ops → _canonicalize_to_dict!
```

## Pass files

| File                       | Defines                                          |
|----------------------------|--------------------------------------------------|
| `canonical_sort.jl`        | `_partial_sort!`, `_site_compare` driver         |
| `commute.jl`               | `_commute_ops` pass                              |
| `reduce.jl`                | `_reduce_ops` pass                               |
| `expand_completeness.jl`   | `_expand_gs_ops` pass + public `expand_completeness` |
| `substitute.jl`            | `_substitute_ops` pass + public `substitute`     |
| `canonicalize.jl`          | `_canonicalize_to_dict!` terminal sink           |
| `weyl.jl`                  | `normal_to_symmetric`, `symmetric_to_normal`     |
| `pipelines.jl`             | `_stream!`, `_canonicalize!`, `_emit_product!`, `_accumulate_with_diag!`, `Base.:*`, `Base.adjoint`, `normal_order`, `simplify`, `commutator` |

## Debugging a pass

Pass functions take `(ops, c, sink::F) where {F}`. To inspect what a pass emits
without running the full pipeline, pass a Vector-pushing sink:

```julia
using SecondQuantizedAlgebra: _commute_ops, CNum, _CNUM_ONE

emitted = Tuple{Vector{QSym}, CNum}[]
_commute_ops(ops_input, _CNUM_ONE) do ops, c
    push!(emitted, (copy(ops), c))    # copy! aliasing rules below
end
display(emitted)
```

## Aliasing rules

- **Single-output passes** (`_reduce_ops`, `_substitute_ops` non-splice path):
  mutate `ops` in place; sink receives the same Vector.
- **Fork passes** (`_commute_ops`, `_expand_gs_ops`, `_substitute_ops` splice
  path): every branch except possibly one calls `copy(ops)` before mutating.
  No two branches share a mutable `ops` Vector.
- **The terminal sink** (`_canonicalize_to_dict!`) takes ownership of the `ops`
  it receives. Callers must not mutate after passing.

When in doubt, `copy(ops)` in your sink. Passes specialize on the sink type via
`where {F}` so the cost is just the copy itself, not a closure allocation.

## Operator hooks

Each concrete `QSym` subtype implements five methods (defined next to the
operator in `src/operators/`):

| Hook                    | Purpose                                      |
|-------------------------|----------------------------------------------|
| `_can_commute(a, b)`    | true iff swap requires no residual           |
| `_commute_pair(a, b)`   | `(b, a, residual_coeff)` for non-commuting   |
| `_reduce_pair(a, b)`    | local identity for same-site pair            |
| `_ground_state_expand`  | for σᵍᵍ: returns `(g, n_levels, site_index)` |
| `_site_compare(a, b, ne)` | three-way `SiteCmp`                         |

Cross-type pairs use generic fallbacks (always distinct sites, always commute).
See `src/operators/operators.jl`.
```

- [ ] **Step 2: Pause**

### Task 12.3: Run formatter and full verification

**Files:**
- All `src/**.jl`

- [ ] **Step 1: Format**

Run:
```bash
make format
```
Expected: clean output, no diffs (or only formatting-driven diffs).

- [ ] **Step 2: Run full test suite**

Run:
```bash
make test 2>&1 | tail -30
```
Expected: all green.

- [ ] **Step 3: Run benchmarks vs baseline**

Run:
```bash
make benchlocal 2>&1 | tee benchmark/post-rewrite.txt
diff benchmark/baseline-redesign-v2.txt benchmark/post-rewrite.txt
```
Expected: medians and allocation counts within noise. Investigate any >10% regression.

- [ ] **Step 4: Run quality gates explicitly**

Run:
```bash
julia --project=. test/jet_test.jl
julia --project=. test/aqua_test.jl
julia --project=. test/explicit_imports_test.jl
julia --project=. test/concrete_test.jl
```
Expected: all pass.

- [ ] **Step 5: Remove the temporary smoke test from Task 2.7**

Run:
```bash
rm -f test/operator_hooks_test.jl
```

(The functionality is covered by `test/canonical_invariant_test.jl` and the per-pass tests in `test/passes_test.jl`.)

- [ ] **Step 6: Pause — final user review**

User reviews the entire rewrite branch:
```bash
git log redesign-v2..src-pipeline-rewrite --stat
git diff redesign-v2..src-pipeline-rewrite --stat
```

User decides whether to merge into `redesign-v2` or open a PR for separate review.

---

## Verification gates (run after every task that touches `src/`)

- `make test` — all `test/*.jl` pass.
- `make format` — clean.
- `julia --project=. test/jet_test.jl` — zero JET issues.
- `julia --project=. test/aqua_test.jl` — Aqua clean.
- `julia --project=. test/explicit_imports_test.jl` — clean.
- `julia --project=. test/concrete_test.jl` — passes.
- `make benchlocal` — medians and allocation counts within noise of `benchmark/baseline-redesign-v2.txt` (run after Task 12.3 Step 3 for the holistic check; per-task benchmarks not required).
- `test/canonical_invariant_test.jl` (after Task 8.1) — passes.
- `test/dissipator_regression_test.jl` (after Task 9.1) — passes.

Any gate failure blocks the task. If the bench gate fails by >10% on a single operation, decide whether to ship the `_known_canonical` flag mitigation described in the spec or investigate further.

---

## Out of plan: deferred follow-ups

These are documented in the spec but **not** executed by this plan. Each gets its own future spec if needed:

- `_known_canonical::Bool` flag on `QTerm` for chained-operation perf.
- Hash-consing of `QSym` leaves (per `TODOs.md`).
- UUID-based `i ≠ j` tracking (per `TODOs.md`, conflicts with hash-consing).
- `SmallDict` for short QAdd sums.
- Benchmark CI gates.
