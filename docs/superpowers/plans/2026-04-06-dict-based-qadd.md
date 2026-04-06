# Dict-Based QAdd Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace `QAdd`'s internal `Vector{QMul}` storage with `Dict{Vector{QSym}, CNum}` so that like-term collection is automatic, equality is order-independent, and `_collect_like_terms` becomes a no-op.

**Architecture:** `QAdd` stores a `Dict{Vector{QSym}, CNum}` mapping operator sequences to prefactors. The `Vector{QMul}` constructor auto-collects like terms and drops zeros. `terms(::QAdd)` yields `QMul` views for computation. `sorted_terms(::QAdd)` yields deterministically-ordered `Vector{QMul}` for printing and TermInterface. The `indices` and `non_equal` fields for symbolic sums are unchanged.

**Tech Stack:** Julia, SymbolicUtils, Symbolics, TermInterface, QuantumOpticsBase

**Known broken tests this resolves:**
- 3× QAdd `isequal` order-independence (`a+σ_i ≠ σ_i+a`)
- 1× `average(2*σ)` structural equality

---

## File Structure

| File | Responsibility | Change |
|------|---------------|--------|
| `src/qadd.jl` | QAdd struct, arithmetic, Σ, expand_sums | **Major rewrite** — new struct, all methods |
| `src/simplify.jl` | Worklist simplifier | Moderate — delete `_collect_like_terms`, adapt iteration |
| `src/commutator.jl` | Commutator computation | Moderate — replace `_ZERO_QADD`, adapt iteration |
| `src/normal_order.jl` | Normal ordering | Small — adapt iteration |
| `src/average.jl` | Average system | Moderate — adapt `average(QAdd)`, `acts_on(QAdd)` |
| `src/index.jl` | Index operations | Small — adapt `change_index`, `get_indices` |
| `src/numeric.jl` | Numeric conversion | Small — adapt iteration |
| `src/interface.jl` | TermInterface | Small — `arguments()` → `sorted_terms()` |
| `src/printing.jl` | Unicode display | Small — use `sorted_terms()` |
| `src/latexify_recipes.jl` | LaTeX rendering | Small — use `sorted_terms()` |
| `test/*.jl` | All test files | Moderate — 130 `.arguments` accesses |

---

### Task 1: Core QAdd struct, constructors, and helpers

**Files:**
- Modify: `src/qadd.jl:1-47`

- [ ] **Step 1: Replace QAdd struct and constructors**

Replace lines 1-22 of `src/qadd.jl` with:

```julia
"""
    QAdd <: QTerm

Sum of [`QMul`](@ref) terms, optionally with summation indices for symbolic sums.

Internally stores a `Dict{Vector{QSym}, CNum}` mapping operator sequences to prefactors.
Like terms are auto-collected on construction. Zero-prefactor terms are dropped.

Fields:
- `dict::Dict{Vector{QSym}, CNum}` — operator sequence → prefactor
- `indices::Vector{Index}` — summation indices (empty = regular sum)
- `non_equal::Vector{Tuple{Index,Index}}` — pairwise inequality constraints
"""
const QTermDict = Dict{Vector{QSym}, CNum}

struct QAdd <: QTerm
    dict::QTermDict
    indices::Vector{Index}
    non_equal::Vector{Tuple{Index, Index}}
end

# Primary constructor from QMul vector — auto-collects like terms, drops zeros
function QAdd(
        terms::Vector{QMul}, indices::Vector{Index},
        non_equal::Vector{Tuple{Index, Index}}
    )
    d = QTermDict()
    for t in terms
        key = t.args_nc
        d[key] = get(d, key, _CNUM_ZERO) + t.arg_c
    end
    filter!(p -> !_iszero_cnum(p.second), d)
    return QAdd(d, indices, non_equal)
end
QAdd(terms::Vector{QMul}) = QAdd(terms, Index[], Tuple{Index, Index}[])
```

- [ ] **Step 2: Add `_iszero_cnum` helper**

Add right after the `_CNUM_NEG_IM` constant in `src/qmul.jl` (around line 20):

```julia
_iszero_cnum(c::CNum) = isequal(real(c), _NUM_ZERO) && isequal(imag(c), _NUM_ZERO)
```

- [ ] **Step 3: Replace length, equality, hash, adjoint**

Replace lines 24-42 of `src/qadd.jl` with:

```julia
Base.length(a::QAdd) = length(a.dict)
Base.iszero(a::QAdd) = isempty(a.dict)

function Base.isequal(a::QAdd, b::QAdd)
    isequal(a.dict, b.dict) || return false
    a.indices == b.indices || return false
    a.non_equal == b.non_equal || return false
    return true
end
Base.:(==)(a::QAdd, b::QAdd) = isequal(a, b)
Base.hash(q::QAdd, h::UInt) = hash(:QAdd, hash(q.dict, hash(q.indices, hash(q.non_equal, h))))

function Base.adjoint(q::QAdd)
    d = QTermDict()
    for (ops, c) in q.dict
        adj_ops = QSym[adjoint(op) for op in reverse(ops)]
        _site_sort!(adj_ops)
        d[adj_ops] = get(d, adj_ops, _CNUM_ZERO) + conj(c)
    end
    filter!(p -> !_iszero_cnum(p.second), d)
    return QAdd(d, q.indices, q.non_equal)
end
```

- [ ] **Step 4: Add iteration helpers**

Add after the adjoint method:

```julia
"""
    terms(q::QAdd)

Iterate over `QMul` views of each term (unordered). For computation only.
"""
terms(q::QAdd) = (QMul(c, copy(ops)) for (ops, c) in q.dict)

"""
    sorted_terms(q::QAdd) -> Vector{QMul}

Return terms in deterministic order. For printing and TermInterface.
"""
function sorted_terms(q::QAdd)
    isempty(q.dict) && return QMul[]
    pairs = sort!(collect(q.dict); by = p -> map(_sort_key, p.first))
    return QMul[QMul(c, ops) for (ops, c) in pairs]
end

"""
    Base.getindex(q::QAdd, key::Vector{QSym}) -> CNum

Look up the prefactor for a given operator sequence. Returns zero if absent.
"""
Base.getindex(q::QAdd, key::Vector{QSym}) = get(q.dict, key, _CNUM_ZERO)
```

---

### Task 2: Addition operators

**Files:**
- Modify: `src/qadd.jl:49-115`

- [ ] **Step 1: Replace all addition methods**

Replace lines 49-115 (everything from `## Addition` through `## QAdd * ...`) with:

```julia
## Addition — always returns QAdd

# QSym + QSym
Base.:+(a::QSym, b::QSym) = QAdd(QMul[_to_qmul(a), _to_qmul(b)])

# QMul + QMul
Base.:+(a::QMul, b::QMul) = QAdd(QMul[a, b])

# QMul + QSym
Base.:+(a::QMul, b::QSym) = QAdd(QMul[a, _to_qmul(b)])
Base.:+(a::QSym, b::QMul) = b + a

# QAdd + QMul — merge into existing dict
function Base.:+(a::QAdd, b::QMul)
    d = copy(a.dict)
    key = b.args_nc
    d[key] = get(d, key, _CNUM_ZERO) + b.arg_c
    _iszero_cnum(d[key]) && delete!(d, key)
    return QAdd(d, a.indices, a.non_equal)
end
Base.:+(a::QMul, b::QAdd) = b + a

# QAdd + QSym
Base.:+(a::QAdd, b::QSym) = a + _to_qmul(b)
Base.:+(a::QSym, b::QAdd) = b + a

# QAdd + QAdd — merge dicts
function Base.:+(a::QAdd, b::QAdd)
    d = copy(a.dict)
    for (ops, c) in b.dict
        d[ops] = get(d, ops, _CNUM_ZERO) + c
    end
    filter!(p -> !_iszero_cnum(p.second), d)
    indices = vcat(a.indices, b.indices) |> unique
    non_equal = vcat(a.non_equal, b.non_equal) |> unique
    return QAdd(d, indices, non_equal)
end

# QField + Number
Base.:+(a::QSym, b::Number) = QAdd(QMul[_to_qmul(a), _scalar_qmul(b)])
Base.:+(a::Number, b::QSym) = b + a
Base.:+(a::QMul, b::Number) = QAdd(QMul[a, _scalar_qmul(b)])
Base.:+(a::Number, b::QMul) = b + a
Base.:+(a::QAdd, b::Number) = a + _scalar_qmul(b)
Base.:+(a::Number, b::QAdd) = b + a

# Subtraction
Base.:-(a::QAdd) = QAdd(QTermDict(ops => -c for (ops, c) in a.dict), a.indices, a.non_equal)
Base.:-(a::QField, b::QField) = a + (-b)
Base.:-(a::QField, b::Number) = a + (-b)
Base.:-(a::Number, b::QField) = a + (-b)
```

---

### Task 3: Multiplication operators and Σ

**Files:**
- Modify: `src/qadd.jl:117-253`

- [ ] **Step 1: Replace all multiplication methods**

Replace lines 117-153 (the `## QAdd * ...` section through `QAdd / Number`) with:

```julia
## QAdd * ... (distributive)

# QAdd * Number
function Base.:*(a::QAdd, b::Number)
    cb = _to_cnum(b)
    d = QTermDict(ops => c * cb for (ops, c) in a.dict)
    filter!(p -> !_iszero_cnum(p.second), d)
    return QAdd(d, a.indices, a.non_equal)
end
Base.:*(a::Number, b::QAdd) = b * a

# QAdd * QSym
function Base.:*(a::QAdd, b::QSym)
    result = QMul[QMul(c, ops) * b for (ops, c) in a.dict]
    return QAdd(result, a.indices, a.non_equal)
end
function Base.:*(a::QSym, b::QAdd)
    result = QMul[a * QMul(c, ops) for (ops, c) in b.dict]
    return QAdd(result, b.indices, b.non_equal)
end

# QAdd * QMul
function Base.:*(a::QAdd, b::QMul)
    result = QMul[QMul(c, ops) * b for (ops, c) in a.dict]
    return QAdd(result, a.indices, a.non_equal)
end
function Base.:*(a::QMul, b::QAdd)
    result = QMul[a * QMul(c, ops) for (ops, c) in b.dict]
    return QAdd(result, b.indices, b.non_equal)
end

# QAdd * QAdd
function Base.:*(a::QAdd, b::QAdd)
    result = QMul[]
    for (ops_a, c_a) in a.dict, (ops_b, c_b) in b.dict
        push!(result, QMul(c_a, ops_a) * QMul(c_b, ops_b))
    end
    indices = vcat(a.indices, b.indices) |> unique
    non_equal = vcat(a.non_equal, b.non_equal) |> unique
    return QAdd(result, indices, non_equal)
end

# QAdd / Number
Base.:/(a::QAdd, b::Number) = a * inv(b)
```

- [ ] **Step 2: Update Σ(::QAdd, ...) to use dict**

Replace the `Σ(expr::QAdd, ...)` method (around line 193) with:

```julia
function Σ(expr::QAdd, i::Index, non_equal::Vector{Index} = Index[])
    ne_pairs = Tuple{Index, Index}[(i, j) for j in non_equal]
    all_indices = vcat(expr.indices, [i])
    all_ne = vcat(expr.non_equal, ne_pairs)
    return QAdd(expr.dict, all_indices, all_ne)
end
```

- [ ] **Step 3: Update expand_sums to use dict iteration**

Replace the `expand_sums` function (around line 222) — change `for term in s.arguments` to iterate dict:

```julia
function expand_sums(s::QAdd)
    isempty(s.indices) && return s

    result_terms = QMul[]
    result_ne = copy(s.non_equal)

    for (ops, c) in s.dict
        term = QMul(c, ops)
        term_indices = get_indices(term)
        needs_split = false

        for idx in term_indices
            for sum_idx in s.indices
                if idx != sum_idx &&
                        idx.space_index == sum_idx.space_index &&
                        !((sum_idx, idx) in s.non_equal) &&
                        !((idx, sum_idx) in s.non_equal)
                    diag_term = change_index(term, sum_idx, idx)
                    push!(result_terms, diag_term)
                    push!(result_ne, (sum_idx, idx))
                    needs_split = true
                end
            end
        end

        push!(result_terms, term)
    end

    return QAdd(result_terms, copy(s.indices), result_ne)
end
```

---

### Task 4: Simplify integration

**Files:**
- Modify: `src/simplify.jl`

- [ ] **Step 1: Adapt `_qsimplify(::QAdd)` methods**

Replace the iteration patterns at lines 36-40 and 52-56:

```julia
# Line 36-40: change s.arguments → s.dict
function _qsimplify(s::QAdd, ord::OrderingConvention)
    all_terms = QMul[]
    for (ops, c) in s.dict
        append!(all_terms, _simplify_qmul(c, ops, ord))
    end
    return _collect_like_terms(QAdd(all_terms, s.indices, s.non_equal))
end

# Line 52-56: same change
function _qsimplify(s::QAdd, ord::OrderingConvention, h::HilbertSpace)
    all_terms = QMul[]
    for (ops, c) in s.dict
        append!(all_terms, _simplify_qmul(c, ops, ord))
    end
    return _collect_like_terms(_apply_ground_state(QAdd(all_terms, s.indices, s.non_equal), h))
end
```

Note: `_collect_like_terms` is now a no-op (the QAdd constructor auto-collects), but keep it temporarily for safety. Delete it in Step 2.

- [ ] **Step 2: Delete `_collect_like_terms`**

Delete the entire `_collect_like_terms` function (around lines 197-212). Then remove the `_collect_like_terms(...)` wrapper calls in `_qsimplify` from Step 1:

```julia
# Remove the wrapper — QAdd constructor auto-collects
function _qsimplify(s::QAdd, ord::OrderingConvention)
    all_terms = QMul[]
    for (ops, c) in s.dict
        append!(all_terms, _simplify_qmul(c, ops, ord))
    end
    return QAdd(all_terms, s.indices, s.non_equal)
end
```

- [ ] **Step 3: Adapt `Symbolics.expand(::QAdd)`**

Replace the iteration at line 220:

```julia
function Symbolics.expand(s::QAdd; kwargs...)
    result = QMul[
        QMul(_expand_prefactor(c; kwargs...), ops)
            for (ops, c) in s.dict
    ]
    isempty(result) && return QAdd(QMul[QMul(_to_cnum(0), QSym[])], s.indices, s.non_equal)
    return QAdd(result, s.indices, s.non_equal)
end
```

---

### Task 5: Commutator adaptation

**Files:**
- Modify: `src/commutator.jl`

- [ ] **Step 1: Replace `_ZERO_QADD` with function**

Replace line 10:
```julia
# Old: const _ZERO_QADD = QAdd(QMul[QMul(_to_cnum(0), QSym[])])
# New:
_zero_qadd() = QAdd(QTermDict(), Index[], Tuple{Index, Index}[])
```

Then find-and-replace all `_ZERO_QADD` → `_zero_qadd()` in the file (and in any test file that imports it).

- [ ] **Step 2: Replace all `.arguments` iteration**

Every `for a_ in a.arguments` becomes `for a_ in terms(a)`.
Every `QAdd(collapsed.arguments, new_indices, new_ne)` becomes `QAdd(collapsed.dict, new_indices, new_ne)`.
Every `QAdd(all_terms, ...)` stays (the `Vector{QMul}` constructor still works).

Specifically replace these patterns:
- Line 65: `QAdd(collapsed.arguments, ...)` → `QAdd(collapsed.dict, ...)`
- Line 72: `for a_ in a.arguments` → `for a_ in terms(a)`
- Line 85: `QAdd(collapsed.arguments, ...)` → `QAdd(collapsed.dict, ...)`
- Line 91: `for b_ in b.arguments` → `for b_ in terms(b)`
- Line 102: `for b_ in b.arguments` → `for b_ in terms(b)`
- Line 109: `for a_ in a.arguments, b_ in b.arguments` → `for a_ in terms(a), b_ in terms(b)`
- Line 127: `QAdd(collapsed.arguments, ...)` → `QAdd(collapsed.dict, ...)`
- Line 135: `for a_ in a.arguments` → `for a_ in terms(a)`
- Line 150: `QAdd(collapsed.arguments, ...)` → `QAdd(collapsed.dict, ...)`
- Line 158: `for b_ in b.arguments` → `for b_ in terms(b)`
- Line 165: `for t in result.arguments` → `for t in terms(result)` and fix the zero check: `all(iszero, ...)` → `iszero(result)`

---

### Task 6: All remaining src/ files

**Files:**
- Modify: `src/normal_order.jl`, `src/index.jl`, `src/average.jl`, `src/numeric.jl`

- [ ] **Step 1: Adapt `src/normal_order.jl`**

Line 17: `for term in expr.arguments` → `for term in terms(expr)`
Line 19: `append!(all_terms, expanded.arguments)` → `append!(all_terms, collect(terms(expanded)))`
Line 61: `append!(result, expanded.arguments)` → `append!(result, collect(terms(expanded)))`

- [ ] **Step 2: Adapt `src/index.jl`**

Line 93: `for t in s.arguments` → `for (ops, c) in s.dict` and `change_index(t, ...)` → `change_index(QMul(c, ops), ...)`
Line 121: `for t in s.arguments` → `for (ops, c) in s.dict` and `get_indices(t)` → `get_indices(QMul(c, ops))`

Full replacement for `change_index(::QAdd, ...)`:
```julia
function change_index(s::QAdd, from::Index, to::Index)
    new_terms = QMul[change_index(QMul(c, ops), from, to) for (ops, c) in s.dict]
    # ... rest of index/non_equal logic unchanged ...
    return QAdd(new_terms, new_indices, new_ne)
end
```

Full replacement for `get_indices(::QAdd)`:
```julia
function get_indices(s::QAdd)
    inds = copy(s.indices)
    for (ops, c) in s.dict
        for idx in get_indices(QMul(c, ops))
            idx ∉ inds && push!(inds, idx)
        end
    end
    return inds
end
```

- [ ] **Step 3: Adapt `src/average.jl`**

Line 96: `mapreduce(average, +, op.arguments)` → `mapreduce(t -> average(t), +, terms(op))`
Line 112: `QAdd(result.arguments, indices, non_equal)` → `QAdd(result.dict, indices, non_equal)`
Line 258: `for t in op.arguments` → `for t in terms(op)`

- [ ] **Step 4: Adapt `src/numeric.jl`**

Line 69: `for t in s.arguments` → `for (ops, c) in s.dict` and `to_numeric(t, ...)` → `to_numeric(QMul(c, ops), ...)`
Line 169: same pattern for the Dict-parameter version

---

### Task 7: Printing, LaTeX, TermInterface

**Files:**
- Modify: `src/printing.jl`, `src/latexify_recipes.jl`, `src/interface.jl`

- [ ] **Step 1: Adapt `src/printing.jl`**

Replace lines 194-207 (the `show(io, ::QAdd)` body):

```julia
    st = sorted_terms(x)
    isempty(st) && return write(io, "0")
    show(io, st[1])
    for i in 2:length(st)
        term = st[i]
        if _is_real_negative(term.arg_c)
            write(io, " - ")
            show(io, QMul(-term.arg_c, term.args_nc))
        else
            write(io, " + ")
            show(io, term)
        end
    end
```

- [ ] **Step 2: Adapt `src/latexify_recipes.jl`**

Line 105: `x.arguments...` → `sorted_terms(x)...`
Line 107: `x.arguments...` → `sorted_terms(x)...`

- [ ] **Step 3: Adapt `src/interface.jl`**

Line 39: `SymbolicUtils.arguments(a::QAdd) = a.arguments` → `SymbolicUtils.arguments(a::QAdd) = sorted_terms(a)`

The `maketerm` at line 42 stays unchanged — it takes `Vector{QMul}` args and passes to the auto-collecting constructor.

---

### Task 8: Update ALL test files

**Files:**
- Modify: all 11 test files with `.arguments` access (130 occurrences)

- [ ] **Step 1: Global replacement patterns**

Apply these mechanical replacements across all test files:

| Old pattern | New pattern |
|-------------|-------------|
| `length(result.arguments)` | `length(result)` |
| `all(iszero, result.arguments)` | `iszero(result)` |
| `all(iszero, simplify(...).arguments)` | `iszero(simplify(...))` |
| `result.arguments[1].arg_c` | `only(collect(terms(result))).arg_c` |
| `result.arguments[1].args_nc` | `only(collect(terms(result))).args_nc` |
| `result.arguments[1].args_nc[1]` | `only(collect(terms(result))).args_nc[1]` |
| `for term in sj.arguments` | `for term in terms(sj)` |
| `terms.arguments` (in average_test) | use `terms(result)` or iterate dict |

For multi-term results where test checks `result.arguments[1]`, use:
```julia
st = sorted_terms(result)
@test st[1].arg_c == expected
```

- [ ] **Step 2: Fix `test/qadd_test.jl`** (15 occurrences)

All `length(s.arguments)` → `length(s)`.
All `all(x -> ..., s.arguments)` → `all(x -> ..., terms(s))`.

- [ ] **Step 3: Fix `test/simplify_test.jl`** (19 occurrences)

Single-term results: `result.arguments[1]` → `only(collect(terms(result)))`.
Multi-term lengths: `length(result.arguments)` → `length(result)`.
Zero checks: `all(iszero, result.arguments)` → `iszero(result)`.

- [ ] **Step 4: Fix `test/normal_order_test.jl`** (32 occurrences)

Same patterns. Many tests check `result.arguments[1].arg_c` for Pauli products — use `only(collect(terms(result)))`.

- [ ] **Step 5: Fix `test/commutator_test.jl`** (9 occurrences)

Replace `_ZERO_QADD` imports with `_zero_qadd`.
Same `.arguments` patterns.

- [ ] **Step 6: Fix `test/indexing_test.jl`** (26 occurrences)

Same patterns. Also convert the 3 `@test_broken` for QAdd order-independence to `@test`.

- [ ] **Step 7: Fix `test/integration_test.jl`** (9 occurrences), `test/phase_space_test.jl` (3), `test/cluster_test.jl` (4), `test/average_test.jl` (4), `test/operators_test.jl` (1)

Same patterns. The `test/average_test.jl` line 381 `terms.arguments` needs special attention — it's a local variable named `terms` (not the function).

- [ ] **Step 8: Convert `@test_broken` to `@test`**

In `test/indexing_test.jl`:
- `@test_broken isequal(a + σi, σi + a)` → `@test isequal(a + σi, σi + a)`
- `@test_broken isequal(σj + σi, σi + σj)` → `@test isequal(σj + σi, σi + σj)`
- `@test_broken isequal(s1 + qadd, qadd + s1)` → `@test isequal(s1 + qadd, qadd + s1)`

In `test/average_test.jl`:
- `@test_broken isequal(average(2 * σi), 2 * average(σi))` → `@test isequal(average(2 * σi), 2 * average(σi))` (if auto-collection normalizes this)

---

### Task 9: Run tests and fix

- [ ] **Step 1: Run `make test`**

Run the full test suite. Fix any remaining failures from missed `.arguments` patterns or semantic changes.

- [ ] **Step 2: Verify quality gates**

Check that `aqua_test`, `concrete_test`, `explicit_imports_test` all pass. The `QTermDict = Dict{Vector{QSym}, CNum}` has concrete key and value types.

- [ ] **Step 3: Update allocation bounds**

Dict operations have different allocation profiles than Vector. Update `@allocations` thresholds in `qadd_test.jl`, `commutator_test.jl`, etc. if needed.

- [ ] **Step 4: Update TODOs.md**

Mark the Dict-based QAdd as done. Mark the resolved broken tests. Update the `@test_broken` count.

---

## Risk Mitigation

1. **Mutable key aliasing**: `terms()` yields `copy(ops)` to avoid mutating dict keys. `_simplify_qmul` already uses `copy(ops)` before mutation.

2. **Empty QAdd**: `isempty(q.dict)` is the zero check. `_zero_qadd()` returns fresh empty instances (no shared mutable singleton).

3. **Printing determinism**: `sorted_terms()` uses `_sort_key` (already defined) which returns `(Int, Int, Symbol)` tuples — well-defined total order.

4. **Constructor backward compatibility**: The `QAdd(Vector{QMul})` constructor still works — it auto-collects into the Dict. All existing code that builds `QMul` vectors and passes to `QAdd(...)` works unchanged.
