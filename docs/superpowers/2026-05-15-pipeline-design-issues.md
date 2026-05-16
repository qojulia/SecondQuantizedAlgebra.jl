# Design memo: passes pipeline issues found during execution of `2026-05-15-src-pipeline-refactor.md`

**Status of execution:** Tasks 1.1 through 4.1 complete and verified. Tasks 5.1 through 6.1 implemented but block on the four issues below. `test/passes_test.jl` shows 10/13 testsets passing; the 3 failures all trace back to the same root causes.

The plan provides verbatim code for every pass and every hook. The issues below are not implementation mistakes; they are internal contradictions in the plan that surface only when the pipeline runs end to end. Each fix is small. The recommended fixes are listed under "Proposed resolution" for each issue.

---

## Issue 1: Canonical-direction direction is contradictory for Fock

### Evidence

Plan Task 2.2 comment:

> `# Destroy < Create within the same site (canonical normal order).`

Plan Task 5.3 test:

```julia
@testset "_commute_ops: Fock aa† → a†a + 1" begin
    h = FockSpace(:f)
    a = Destroy(h, :a)
    ad = adjoint(a)
    emitted = Tuple{Vector{QSym}, CNum}[]
    _commute_ops(QSym[a, ad], _CNUM_ONE) do o, c
        push!(emitted, (copy(o), c))
    end
    @test length(emitted) == 2
    sort!(emitted, by = e -> length(e[1]))
    @test isempty(emitted[1][1]) && emitted[1][2] == _CNUM_ONE
    @test emitted[2][1] == QSym[ad, a]                            # the swapped branch
end
```

The test asserts `emitted[2][1] == QSym[ad, a]` (Create first). That is **standard physics normal order**, i.e. Create < Destroy. The comment in 2.2 says the opposite.

Project `CLAUDE.md` confirms standard normal order is the intent:

> **Eager ordering produces canonical form**: under `NormalOrder` (default), ... applies all three transformations (commutation swaps, ...).

`NormalOrder` is named precisely after standard normal order. Result: Create on the left, Destroy on the right.

### Proposed resolution

Treat the test as the source of truth and flip the comments + the `Less`/`Greater` returns in Tasks 2.2 hooks:

```julia
# Create < Destroy within the same site (canonical normal order).
function _site_compare(a::Create, b::Destroy, ne::Vector{NonEqualPair})::SiteCmp
    cmp = _site_compare(Destroy(a.name, a.space_index, a.index),
                        Destroy(b.name, b.space_index, b.index), ne)
    cmp === Equal && return Less        # <-- was Greater
    return cmp
end
function _site_compare(a::Destroy, b::Create, ne::Vector{NonEqualPair})::SiteCmp
    cmp = _site_compare(Destroy(a.name, a.space_index, a.index),
                        Destroy(b.name, b.space_index, b.index), ne)
    cmp === Equal && return Greater     # <-- was Less
    return cmp
end
```

This makes `_site_compare(Destroy, Create) === Greater` for same-site, which lets `_partial_sort!` move `[a, ad]` to `[ad, a]`. Same fix applies to `PhaseSpace` cross-type comparisons (`Position` < `Momentum` is the convention there, so leave those as the plan wrote them and verify the test passes; if not, flip).

---

## Issue 2: `_commute_ops` trigger predicate is unreachable for cross-type Fock

### Evidence

Plan Task 5.3 implementation:

```julia
cmp = _site_compare(a, b, _EMPTY_NE)
if cmp === Equal && !_can_commute(a, b)
    ...
```

`_site_compare(::Destroy, ::Create)` and `_site_compare(::Create, ::Destroy)` are explicitly defined to return `Less` or `Greater`, never `Equal`. So the predicate `cmp === Equal && ...` never fires for the Fock ladder pair. Same for `PhaseSpace` (`Position` vs `Momentum`).

### Proposed resolution

Two acceptable shapes; pick one.

**(A)** Run `_partial_sort!` first so commute always sees in-canonical-order pairs, then fire commute on `cmp === Less && !_can_commute(a, b)`:

```julia
# Pre-condition: ops is already partial-sorted. A non-commuting same-site pair
# must therefore be in canonical order (Less); firing produces the swap branch.
if cmp === Less && !_can_commute(a, b)
    ...
end
```

This works because: after partial sort, the only way two same-site operators end up adjacent and `_can_commute` is false is the canonical-order case. The "in order, with residual" case.

**(B)** Separate "is same site?" from "what is the canonical order?" by giving each operator type a `_is_same_site(a, b)::Bool` hook and use that as the gate:

```julia
if _is_same_site(a, b) && !_can_commute(a, b)
    ...
```

`_site_compare` then only encodes ordering. This is cleaner but changes the hook surface area in Tasks 2.1 through 2.6.

**Recommended:** (A) is the smaller diff. Just one line in `_commute_ops_at`. (B) is the better long-term design and worth doing if we are happy to revise Step 2.

---

## Issue 3: Pipeline order `commute -> reduce` is wrong for Transition and Pauli

### Evidence

Plan Task 6.1 `_stream!`:

```julia
@inline function _stream!(out, ops, c, ne)
    _commute_ops(ops, c) do ops1, c1
        _reduce_ops(ops1, c1) do ops2, c2
            _canonicalize_to_dict!(out, ops2, c2, ne)
        end
    end
end
```

Plan Tasks 2.3, 2.4 hooks:

```julia
_can_commute(::Transition, ::Transition) = false
_commute_pair(::Transition, ::Transition) = error("unreachable: Transition uses _reduce_pair, not _commute_pair")

_can_commute(::Pauli, ::Pauli) = false
_commute_pair(::Pauli, ::Pauli) = error("unreachable: Pauli uses _reduce_pair")
```

If `_commute_ops` fires on any same-site Transition pair (which `!_can_commute === true` says it should), the pipeline hits `error(...)` and crashes. The plan's stated semantics are "Transition composition happens in reduce, not commute," but the pipeline runs commute first.

### Proposed resolution

Flip the pass order in `_stream!`:

```julia
@inline function _stream!(out, ops, c, ne)
    _reduce_ops(ops, c) do ops1, c1
        _commute_ops(ops1, c1) do ops2, c2
            _canonicalize_to_dict!(out, ops2, c2, ne)
        end
    end
end
```

Then `_reduce_ops` handles Transition/Pauli first (collapsing the pair), and `_commute_ops` only ever sees Fock/Spin/PhaseSpace ladder pairs.

**Caveat:** `_commute_ops` can produce a residual branch with empty ops or with a contracted op. The residual itself may then be reducible (or even contain a `σᵍᵍ` for the contracted-op case). The cleaner shape is:

```
reduce -> commute -> reduce -> sink
```

The terminal `reduce` catches any pair that the commute residual produced. For the test cases in the plan, the simpler `reduce -> commute` works; for fully general use the double-reduce is safer. Decide based on what regression coverage we want to claim. Recommend the double-reduce.

---

## Issue 4: `_commute_pair(Spin, Spin)` residual shape is undefined

### Evidence

Plan Task 2.5:

```julia
function _commute_pair(a::Spin, b::Spin)
    ...
    eps = _levi_civita[a.axis][b.axis]
    contracted = Spin(a.name, 6 - a.axis - b.axis, a.space_index, a.index)
    return (b, a, _mul_cnum(_to_cnum(im * eps), _CNUM_ONE))
end
```

The local `contracted` is computed and discarded. The returned tuple is `(b, a, residual_coeff)`. Plan note:

> The Spin design represents the commutator residual differently from Fock: the residual contracts to a single op, not the identity. The `_commute_pair` returns `(b, a, residual_coeff_for_contracted_op)`. The pipeline needs to handle this; clarified in `_commute_ops` design.

But `_commute_ops_at` in Task 5.3 treats the residual as a scalar applied to the *empty-ops* branch (where the pair is just deleted). That matches Fock and PhaseSpace where the residual is `coeff * 1`. For Spin the residual is `coeff * contracted_op`, which the pipeline does not currently emit.

Running `Sy * Sx` through the current pipeline would produce wrong physics (it would emit `Sx Sy + i*epsilon` instead of `Sx Sy + i*epsilon * Sz`).

### Proposed resolution

Change `_commute_pair` to return one of two shapes, and have `_commute_ops_at` dispatch:

```julia
# Shape A: (swap_b, swap_a, residual_coeff)            -- residual is coeff * identity
# Shape B: (swap_b, swap_a, residual_coeff, residual_op) -- residual is coeff * residual_op
```

Then in `_commute_ops_at`:

```julia
res = _commute_pair(a, b)
sw, contracted_form_a, residual = res[1], res[2], res[3]
# ... swap branch ...
new_c = _mul_cnum(c, residual)
if length(res) == 4
    # Spin-shaped residual: replace pair with contracted op
    ops[i] = res[4]
    deleteat!(ops, i + 1)
else
    # Identity-shaped residual: delete pair
    deleteat!(ops, i:(i + 1))
end
```

Cleaner: have `_commute_pair` always return a 4-tuple `(swap_b, swap_a, residual_coeff, residual_ops::Vector{QSym})` where Fock/PhaseSpace use the empty vector and Spin uses a one-element vector. The pipeline always splices.

**Recommended:** uniform 4-tuple with a `Vector{QSym}` residual. Makes the pipeline branch-free and lets future operator types add multi-op residuals without further pipeline changes.

---

## Summary of recommended plan edits

| Issue | Plan location | Edit |
|---|---|---|
| 1 | Task 2.2 | Flip `Less`/`Greater` in cross-type Fock `_site_compare` (and verify PhaseSpace). |
| 2 | Task 5.3 | Change trigger to `cmp === Less && !_can_commute(a, b)`. |
| 3 | Task 6.1 | Swap pass order to `reduce -> commute -> reduce -> sink`. |
| 4 | Tasks 2.5 + 5.3 | Uniform 4-tuple `(swap_b, swap_a, residual_coeff, residual_ops::Vector{QSym})` for `_commute_pair`; pipeline always splices the residual. |

All four edits are local (no other tasks ripple). After applying them, `passes_test.jl` should pass cleanly and Step 6 can land on a sound pass layer.

## Optional: stronger canonical-form invariant test up front

The bugs above would have been caught earlier by an explicit "after `*`, every term is canonical" test. Plan Task 8.1 adds one but only at the very end. Moving (a stripped-down version of) `canonical_invariant_test.jl` to land right after Task 6.1 would surface this class of issue much sooner during execution. Worth considering during the plan revision.

## What I will do after the plan is revised

1. Re-read the revised plan sections (1.x, 2.x, 5.x, 6.1).
2. Apply the corrected hooks and pipeline.
3. Re-run `test/passes_test.jl` (should pass cleanly).
4. Proceed through Steps 6 through 12 as written.

No code from the existing work needs to be rolled back. Tasks 1.1 through 4.1 are clean; the algebra/ files (commute.jl, reduce.jl, expand_completeness.jl, substitute.jl, pipelines.jl) need targeted edits only.
