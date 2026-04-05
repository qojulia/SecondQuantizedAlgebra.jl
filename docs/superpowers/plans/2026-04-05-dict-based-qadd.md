# Dict-Based QAdd Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace `QAdd`'s internal `Vector{QMul}` storage with `Dict{Vector{QSym}, CNum}` so that like-term collection is automatic, equality is order-independent, and `_collect_like_terms` becomes a no-op.

**Architecture:** `QAdd` stores a `Dict{Vector{QSym}, CNum}` mapping operator content to prefactor. Construction auto-collects like terms and drops zeros. A `terms(::QAdd)` iterator yields `QMul` views for backward-compatible iteration. Printing and TermInterface use `sorted_terms()` for deterministic output. The `indices` and `non_equal` fields for symbolic sums are unchanged.

**Tech Stack:** Julia, SymbolicUtils, Symbolics, TermInterface

---

## File Structure

| File | Responsibility | Change |
|------|---------------|--------|
| `src/qadd.jl` | QAdd struct, arithmetic, iteration | **Major rewrite** — new struct, all +/* methods, helpers |
| `src/simplify.jl` | Worklist simplifier | Moderate — remove `_collect_like_terms`, adapt `_qsimplify` |
| `src/commutator.jl` | Commutator computation | Moderate — iteration pattern changes |
| `src/normal_order.jl` | Ground state rewriting | Moderate — iteration pattern changes |
| `src/interface.jl` | TermInterface integration | Small — `arguments()` and `maketerm()` |
| `src/printing.jl` | Unicode display | Small — use `sorted_terms()` |
| `src/latexify_recipes.jl` | LaTeX rendering | Small — use `sorted_terms()` |
| `src/numeric.jl` | Numeric conversion | Small — iteration pattern |
| `src/change_index.jl` | Index substitution | Small — iteration pattern |
| `src/expand_sums.jl` | Diagonal splitting | Small — iteration pattern |
| `src/qsum.jl` | Symbolic sums (Σ) | Small — construction pattern |
| `test/qadd_test.jl` | QAdd unit tests | **Major rewrite** — new semantics |
| `test/simplify_test.jl` | Simplify tests | Moderate — adapt assertions |
| `test/normal_order_test.jl` | Normal order tests | Moderate — adapt assertions |
| `test/commutator_test.jl` | Commutator tests | Small — adapt assertions |
| Other test files | Various | Small — adapt `.arguments[i]` patterns |

---

## Design Decisions

### Key type alias

```julia
const QTermDict = Dict{Vector{QSym}, CNum}
```

### Constructor contract

```julia
struct QAdd <: QTerm
    dict::QTermDict
    indices::Vector{Index}
    non_equal::Vector{Tuple{Index, Index}}
end
```

The inner constructor:
1. Merges like terms (same `Vector{QSym}` key → sum prefactors)
2. Drops zero-prefactor entries
3. Guarantees: empty dict = scalar zero

### Iteration helpers

```julia
# Iterate over QMul views (unordered) — for computation
terms(q::QAdd) = (QMul(c, ops) for (ops, c) in q.dict)

# Sorted iteration — for printing and TermInterface
function sorted_terms(q::QAdd)
    pairs = collect(q.dict)
    sort!(pairs; by = first, lt = _ops_lt)
    return [QMul(c, ops) for (ops, c) in pairs]
end
```

### Test helpers

Tests currently do `result.arguments[1].arg_c`. With Dict storage, term order is non-deterministic. Tests should use:
- `length(result)` for term count
- `result[ops_key]` for looking up a specific term's prefactor
- `has_term(result, ops, coeff)` helper for asserting specific terms exist
- `scalar_part(result)` for the coefficient of the identity (empty ops) term

---

### Task 1: Core QAdd struct and constructors

**Files:**
- Modify: `src/qadd.jl` (full rewrite of struct + constructors)
- Test: `test/qadd_test.jl`

- [ ] **Step 1: Write failing tests for new QAdd semantics**

Add to `test/qadd_test.jl`:
```julia
@testset "Dict-based construction" begin
    h = FockSpace(:c)
    a = Destroy(h, :a)
    ad = a'

    # Auto like-term collection
    s = QAdd(QMul[QMul(2, QSym[a]), QMul(3, QSym[a])])
    @test length(s) == 1  # merged into 5*a

    # Zero terms dropped
    s2 = QAdd(QMul[QMul(0, QSym[a]), QMul(3, QSym[ad])])
    @test length(s2) == 1

    # Order-independent equality
    s3 = QAdd(QMul[QMul(1, QSym[a]), QMul(2, QSym[ad])])
    s4 = QAdd(QMul[QMul(2, QSym[ad]), QMul(1, QSym[a])])
    @test isequal(s3, s4)
    @test hash(s3) == hash(s4)

    # Empty dict = zero
    s5 = QAdd(QMul[QMul(0, QSym[a])])
    @test iszero(s5)
end
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `make test` (or `julia --project -e 'using Pkg; Pkg.test()'`)
Expected: FAIL — current QAdd doesn't auto-collect or have order-independent equality.

- [ ] **Step 3: Rewrite QAdd struct and constructors**

Replace the entire struct and constructor section in `src/qadd.jl`:

```julia
const QTermDict = Dict{Vector{QSym}, CNum}

struct QAdd <: QTerm
    dict::QTermDict
    indices::Vector{Index}
    non_equal::Vector{Tuple{Index, Index}}
    function QAdd(dict::QTermDict, indices::Vector{Index},
        non_equal::Vector{Tuple{Index, Index}})
        return new(dict, indices, non_equal)
    end
end

# Primary constructor from QMul vector — auto-collects like terms
function QAdd(terms::Vector{QMul}, indices::Vector{Index},
    non_equal::Vector{Tuple{Index, Index}})
    d = QTermDict()
    for t in terms
        key = t.args_nc
        d[key] = get(d, key, _to_cnum(0)) + t.arg_c
    end
    # Drop zero entries
    filter!(p -> !iszero(p.second), d)
    return QAdd(d, indices, non_equal)
end
QAdd(terms::Vector{QMul}) = QAdd(terms, Index[], Tuple{Index, Index}[])

# Convenience: construct from a single dict entry
function QAdd(ops::Vector{QSym}, c::CNum,
    indices::Vector{Index} = Index[],
    non_equal::Vector{Tuple{Index, Index}} = Tuple{Index, Index}[])
    d = QTermDict(ops => c)
    return QAdd(d, indices, non_equal)
end
```

- [ ] **Step 4: Add length, iszero, zero, iteration helpers**

```julia
Base.length(a::QAdd) = length(a.dict)
Base.iszero(a::QAdd) = isempty(a.dict)
Base.zero(::QAdd) = QAdd(QTermDict(), Index[], Tuple{Index, Index}[])
Base.zero(::Type{QAdd}) = QAdd(QTermDict(), Index[], Tuple{Index, Index}[])

# Look up prefactor for a given operator sequence
function Base.getindex(a::QAdd, key::Vector{QSym})
    return get(a.dict, key, _to_cnum(0))
end

# Iterate over (ops, coeff) pairs
Base.iterate(q::QAdd, state...) = iterate(q.dict, state...)

# Yield QMul views (unordered) — main iteration for computation
terms(q::QAdd) = (QMul(c, ops) for (ops, c) in q.dict)

# Deterministic sort key for operator sequences
function _ops_lt(a::Pair{Vector{QSym}, CNum}, b::Pair{Vector{QSym}, CNum})
    ka = map(_sort_key, a.first)
    kb = map(_sort_key, b.first)
    length(ka) != length(kb) && return length(ka) < length(kb)
    for (x, y) in zip(ka, kb)
        x != y && return x < y
    end
    return false
end

# Sorted QMul vector — for printing and TermInterface
function sorted_terms(q::QAdd)
    pairs = sort!(collect(q.dict); lt = _ops_lt)
    return QMul[QMul(c, ops) for (ops, c) in pairs]
end
```

- [ ] **Step 5: Add equality and hashing**

```julia
function Base.isequal(a::QAdd, b::QAdd)
    isequal(a.dict, b.dict) || return false
    a.indices == b.indices || return false
    a.non_equal == b.non_equal || return false
    return true
end
Base.:(==)(a::QAdd, b::QAdd) = isequal(a, b)
Base.hash(q::QAdd, h::UInt) = hash(:QAdd, hash(q.dict, hash(q.indices, hash(q.non_equal, h))))
```

- [ ] **Step 6: Add adjoint**

```julia
function Base.adjoint(q::QAdd)
    d = QTermDict()
    for (ops, c) in q.dict
        adj_ops = QSym[adjoint(op) for op in reverse(ops)]
        _site_sort!(adj_ops)
        d[adj_ops] = get(d, adj_ops, _to_cnum(0)) + conj(c)
    end
    filter!(p -> !iszero(p.second), d)
    return QAdd(d, q.indices, q.non_equal)
end
```

- [ ] **Step 7: Run tests to verify new construction tests pass**

Run: `make test`
Expected: New construction tests PASS. Existing tests will FAIL (addressed in later tasks).

---

### Task 2: Addition operators

**Files:**
- Modify: `src/qadd.jl` (replace all `+` methods)
- Test: `test/qadd_test.jl`

- [ ] **Step 1: Write failing tests for auto-collection on addition**

```julia
@testset "Addition auto-collects like terms" begin
    h = FockSpace(:c)
    a = Destroy(h, :a)
    ad = a'

    # a + a = 2a (auto-collected)
    s = a + a
    @test length(s) == 1
    @test s[QSym[a]] == 2

    # a + ad + a = 2a + ad
    s2 = a + ad + a
    @test length(s2) == 2

    # Cancellation: a + (-a) = 0
    s3 = a + (-a)
    @test iszero(s3)

    # Commutativity: a + ad == ad + a
    @test isequal(a + ad, ad + a)
end
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `make test`
Expected: FAIL

- [ ] **Step 3: Implement addition operators**

Replace all `+` methods in `src/qadd.jl`. The key insight: every `+` builds a `QMul` vector and passes to the constructor (which auto-collects).

```julia
# Helpers: wrap QSym/scalar as QMul
_to_qmul(a::QSym) = QMul(_to_cnum(1), QSym[a])
_to_qmul(a::QMul) = a
_scalar_qmul(x::Number) = QMul(_to_cnum(x), QSym[])

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
    d[b.args_nc] = get(d, b.args_nc, _to_cnum(0)) + b.arg_c
    iszero(d[b.args_nc]) && delete!(d, b.args_nc)
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
        d[ops] = get(d, ops, _to_cnum(0)) + c
    end
    filter!(p -> !iszero(p.second), d)
    indices = vcat(a.indices, b.indices) |> unique
    non_equal = vcat(a.non_equal, b.non_equal) |> unique
    return QAdd(d, indices, non_equal)
end

# QField + Number
Base.:+(a::QSym, b::Number) = QAdd(QMul[_to_qmul(a), _scalar_qmul(b)])
Base.:+(a::Number, b::QSym) = b + a

Base.:+(a::QMul, b::Number) = QAdd(QMul[a, _scalar_qmul(b)])
Base.:+(a::Number, b::QMul) = b + a

function Base.:+(a::QAdd, b::Number)
    return a + _scalar_qmul(b)
end
Base.:+(a::Number, b::QAdd) = b + a

# Subtraction
Base.:-(a::QAdd) = QAdd(QTermDict(ops => -c for (ops, c) in a.dict), a.indices, a.non_equal)
Base.:-(a::QField, b::QField) = a + (-b)
Base.:-(a::QField, b::Number) = a + (-b)
Base.:-(a::Number, b::QField) = a + (-b)
```

- [ ] **Step 4: Run tests to verify addition tests pass**

Run: `make test`
Expected: Addition tests PASS.

---

### Task 3: Multiplication operators (distributive)

**Files:**
- Modify: `src/qadd.jl` (replace all `*` methods involving QAdd)
- Test: `test/qadd_test.jl`

- [ ] **Step 1: Write failing tests for distributive multiplication**

```julia
@testset "Distributive multiplication with auto-collection" begin
    h = FockSpace(:c)
    a = Destroy(h, :a)
    ad = a'

    # (a + ad) * 3 — scale each term
    s = (a + ad) * 3
    @test length(s) == 2
    @test s[QSym[a]] == 3
    @test s[QSym[ad]] == 3

    # (a + a) * ad — auto-collected before distribute
    s2 = (a + a) * ad
    @test length(s2) == 1  # 2*a*ad, single term
end
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `make test`
Expected: FAIL

- [ ] **Step 3: Implement multiplication operators**

```julia
## QAdd * ... (distributive)

# QAdd * Number
function Base.:*(a::QAdd, b::Number)
    cb = _to_cnum(b)
    d = QTermDict(ops => c * cb for (ops, c) in a.dict)
    return QAdd(d, a.indices, a.non_equal)
end
Base.:*(a::Number, b::QAdd) = b * a

# QAdd * QSym
function Base.:*(a::QAdd, b::QSym)
    result = QMul[]
    for (ops, c) in a.dict
        push!(result, QMul(c, ops) * b)
    end
    return QAdd(result, a.indices, a.non_equal)
end
function Base.:*(a::QSym, b::QAdd)
    result = QMul[]
    for (ops, c) in b.dict
        push!(result, a * QMul(c, ops))
    end
    return QAdd(result, b.indices, b.non_equal)
end

# QAdd * QMul
function Base.:*(a::QAdd, b::QMul)
    result = QMul[]
    for (ops, c) in a.dict
        push!(result, QMul(c, ops) * b)
    end
    return QAdd(result, a.indices, a.non_equal)
end
function Base.:*(a::QMul, b::QAdd)
    result = QMul[]
    for (ops, c) in b.dict
        push!(result, a * QMul(c, ops))
    end
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

- [ ] **Step 4: Run tests to verify multiplication tests pass**

Run: `make test`
Expected: Multiplication tests PASS.

---

### Task 4: Simplify integration

**Files:**
- Modify: `src/simplify.jl`
- Test: `test/simplify_test.jl`

- [ ] **Step 1: Adapt `_qsimplify` and remove `_collect_like_terms`**

The `_collect_like_terms` function is now redundant — the QAdd constructor auto-collects. Replace the simplify functions:

```julia
function _qsimplify(op::QSym, ::OrderingConvention)
    return QAdd(QMul[QMul(_to_cnum(1), QSym[op])])
end

function _qsimplify(m::QMul, ord::OrderingConvention)
    raw_terms = _simplify_qmul(m.arg_c, m.args_nc, ord)
    return QAdd(raw_terms)  # constructor auto-collects
end

function _qsimplify(s::QAdd, ord::OrderingConvention)
    all_terms = QMul[]
    for (ops, c) in s.dict
        append!(all_terms, _simplify_qmul(c, ops, ord))
    end
    return QAdd(all_terms, s.indices, s.non_equal)  # constructor auto-collects
end
```

Do the same for the `_qsimplify` methods with `HilbertSpace` argument — remove the `_collect_like_terms` wrapper:

```julia
function _qsimplify(op::QSym, ord::OrderingConvention, h::HilbertSpace)
    return _apply_ground_state(_qsimplify(op, ord), h)
end
function _qsimplify(m::QMul, ord::OrderingConvention, h::HilbertSpace)
    raw_terms = _simplify_qmul(m.arg_c, m.args_nc, ord)
    return _apply_ground_state(QAdd(raw_terms), h)
end
function _qsimplify(s::QAdd, ord::OrderingConvention, h::HilbertSpace)
    all_terms = QMul[]
    for (ops, c) in s.dict
        append!(all_terms, _simplify_qmul(c, ops, ord))
    end
    return _apply_ground_state(QAdd(all_terms, s.indices, s.non_equal), h)
end
```

Delete the `_collect_like_terms` function entirely (lines 191-199 of the current file).

- [ ] **Step 2: Adapt `Symbolics.expand`**

```julia
function Symbolics.expand(s::QAdd; kwargs...)
    result = QMul[]
    for (ops, c) in s.dict
        push!(result, QMul(_expand_prefactor(c; kwargs...), ops))
    end
    return QAdd(result, s.indices, s.non_equal)
end
function Symbolics.expand(m::QMul; kwargs...)
    return Symbolics.expand(QAdd(QMul[m]); kwargs...)
end
function Symbolics.expand(op::QSym; kwargs...)
    return QAdd(QMul[QMul(_to_cnum(1), QSym[op])])
end
```

- [ ] **Step 3: Update simplify tests**

Tests that check `result.arguments[1]` need to use the new helpers. Key pattern change:

```julia
# OLD:
@test result.arguments[1].arg_c == 5

# NEW (option A — use terms() and find):
t = only(collect(terms(result)))  # when expecting exactly 1 term
@test t.arg_c == 5

# NEW (option B — use Dict lookup):
@test result[QSym[ad, a]] == 5
```

Update each test in `test/simplify_test.jl` following this pattern. For tests that check `length(result.arguments)`, replace with `length(result)`.

- [ ] **Step 4: Run simplify tests**

Run: `make test`
Expected: PASS

---

### Task 5: Commutator adaptation

**Files:**
- Modify: `src/commutator.jl`
- Test: `test/commutator_test.jl`

- [ ] **Step 1: Replace `_ZERO_QADD` with a function**

```julia
_zero_qadd() = QAdd(QTermDict(), Index[], Tuple{Index, Index}[])
```

Replace all references to `_ZERO_QADD` with `_zero_qadd()`. This avoids the shared-mutable-singleton issue (old issue #5).

- [ ] **Step 2: Adapt iteration in commutator methods**

Every `for a_ in a.arguments` becomes `for a_ in terms(a)`. Every `QAdd(collapsed.arguments, ...)` becomes `QAdd(collapsed.dict, ...)` (or reconstruct from the dict).

Key pattern for the `_append_terms!` helper:

```julia
function _append_terms!(all_terms::Vector{QMul}, result::QAdd)
    for (ops, c) in result.dict
        push!(all_terms, QMul(c, ops))
    end
    return all_terms
end
```

The diagonal collapse code at lines 62-66, 82-86, 124-129, 147-152 reconstructs a QAdd with new indices. Adapt to use the dict:

```julia
collapsed = change_index(a, idx, b.index)
collapsed_qadd = QAdd(collapsed.dict, new_indices, new_ne)
```

- [ ] **Step 3: Update commutator tests**

Replace `result.arguments[i]` access patterns with `terms()`-based or Dict-lookup-based assertions.

- [ ] **Step 4: Run commutator tests**

Run: `make test`
Expected: PASS

---

### Task 6: Normal order, ground state, and remaining src/ files

**Files:**
- Modify: `src/normal_order.jl`, `src/change_index.jl`, `src/expand_sums.jl`, `src/qsum.jl`, `src/numeric.jl`
- Test: `test/normal_order_test.jl`

- [ ] **Step 1: Adapt normal_order.jl**

`_apply_ground_state` and `_expand_ground_state` iterate over terms:

```julia
function _apply_ground_state(expr::QAdd, h::NLevelSpace)
    all_terms = QMul[]
    for (ops, c) in expr.dict
        expanded = _expand_ground_state(QMul(c, ops), h)
        append!(all_terms, collect(terms(expanded)))
    end
    return QAdd(all_terms, expr.indices, expr.non_equal)
end
```

`_expand_ground_state` returns a `QAdd`, so its callers collecting into `all_terms` need:

```julia
for t in terms(expanded)
    push!(all_terms, t)
end
```

- [ ] **Step 2: Adapt change_index.jl**

```julia
function change_index(s::QAdd, from::Index, to::Index)
    new_terms = QMul[change_index(QMul(c, ops), from, to) for (ops, c) in s.dict]
    return QAdd(new_terms, s.indices, s.non_equal)
end
```

Also adapt `change_index(s::QAdd, ...)` for the version iterating `s.arguments` at line 79.

- [ ] **Step 3: Adapt expand_sums.jl**

```julia
function expand_sums(s::QAdd)
    isempty(s.indices) && return s
    result_terms = QMul[]
    result_ne = copy(s.non_equal)
    for (ops, c) in s.dict
        term = QMul(c, ops)
        # ... same logic, using term instead of iterating arguments
        push!(result_terms, term)
    end
    return QAdd(result_terms, copy(s.indices), result_ne)
end
```

- [ ] **Step 4: Adapt qsum.jl**

```julia
function Σ(expr::QAdd, i::Index, non_equal::Vector{Index} = Index[])
    ne_pairs = Tuple{Index, Index}[(i, j) for j in non_equal]
    all_indices = vcat(expr.indices, [i])
    all_ne = vcat(expr.non_equal, ne_pairs)
    return QAdd(expr.dict, all_indices, all_ne)  # pass dict directly
end
```

- [ ] **Step 5: Adapt numeric.jl**

```julia
function to_numeric(s::QAdd, b::QuantumOpticsBase.Basis; kwargs...)
    terms_num = [to_numeric(QMul(c, ops), b; kwargs...) for (ops, c) in s.dict]
    return sum(terms_num)
end
```

- [ ] **Step 6: Update normal_order tests**

Replace `result.arguments[1].arg_c` patterns. For Pauli tests that check a specific term's prefactor and axis, use:

```julia
# OLD: result.arguments[1].arg_c == im
# NEW:
t = only(collect(terms(result)))
@test t.arg_c == im
@test t.args_nc[1].axis == 3
```

For tests checking `length(result.arguments)`, replace with `length(result)`.

- [ ] **Step 7: Run all affected tests**

Run: `make test`
Expected: PASS

---

### Task 7: Printing and LaTeX

**Files:**
- Modify: `src/printing.jl`, `src/latexify_recipes.jl`
- Test: `test/printing_test.jl`, `test/latexify_test.jl`

- [ ] **Step 1: Adapt printing.jl**

Replace the `show(io, x::QAdd)` method (around line 194):

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
    return
```

- [ ] **Step 2: Adapt latexify_recipes.jl**

Replace `x.arguments` with `sorted_terms(x)`:

```julia
    if !isempty(x.indices)
        # ... index prefix logic unchanged ...
        st = sorted_terms(x)
        return Expr(:latexifymerge, prefix, " ", Expr(:call, :+, st...))
    end
    return Expr(:call, :+, sorted_terms(x)...)
```

- [ ] **Step 3: Run printing and LaTeX tests**

Run: `make test`
Expected: PASS

---

### Task 8: TermInterface integration

**Files:**
- Modify: `src/interface.jl`
- Test: `test/interface_test.jl`

- [ ] **Step 1: Adapt arguments() and maketerm()**

```julia
# QAdd — sums
SymbolicUtils.iscall(::QAdd) = true
SymbolicUtils.iscall(::Type{QAdd}) = true
SymbolicUtils.operation(::QAdd) = (+)
SymbolicUtils.arguments(a::QAdd) = sorted_terms(a)
TermInterface.metadata(::QAdd) = nothing

function TermInterface.maketerm(::Type{QAdd}, ::typeof(+), args, metadata)
    muls = QMul[x isa QMul ? x : QMul(_to_cnum(1), QSym[x]) for x in args]
    return QAdd(muls)
end
```

- [ ] **Step 2: Run interface tests**

Run: `make test`
Expected: PASS

---

### Task 9: Update remaining test files

**Files:**
- Modify: `test/integration_test.jl`, `test/phase_space_test.jl`, `test/indexing_test.jl`, `test/cluster_test.jl`

- [ ] **Step 1: Find and replace `.arguments` patterns**

In each test file, replace:
- `result.arguments` → `collect(terms(result))` or `sorted_terms(result)` (when order matters for indexing)
- `length(result.arguments)` → `length(result)`
- `result.arguments[1].arg_c` → use `only(collect(terms(result))).arg_c` (for single-term results) or Dict lookup
- `all(iszero, result.arguments)` → `iszero(result)`

- [ ] **Step 2: Run full test suite**

Run: `make test`
Expected: ALL PASS

---

### Task 10: Quality gates and performance verification

**Files:**
- Test: all test files

- [ ] **Step 1: Run full quality gate suite**

Run: `make test`
Verify:
- `aqua_test.jl` passes (no stale deps, undefined exports)
- `jet_test.jl` passes (no type errors)
- `concrete_test.jl` passes (QAdd fields are concretely typed: `Dict{Vector{QSym}, CNum}` is concrete)
- `explicit_imports_test.jl` passes

- [ ] **Step 2: Verify type stability**

All `@inferred` tests in `qadd_test.jl`, `simplify_test.jl`, `normal_order_test.jl` must pass. The Dict-based constructor should be inferrable since all types are concrete.

- [ ] **Step 3: Verify allocation bounds**

Run allocation tests. Expected changes:
- `QAdd + QMul` allocations may increase slightly (Dict copy + insert vs Vector push)
- `simplify` allocations should decrease (no `_collect_like_terms` pass)
- `QAdd * QAdd` allocations should be comparable

If allocation bounds are exceeded, adjust the thresholds in tests — the absolute counts may shift but should remain in the same order of magnitude.

- [ ] **Step 4: Run benchmarks**

Run: `make bench`
Compare against baseline. Key metrics:
- `simplify((a * a')^n)` for n=2,3,4 — should improve (no separate collect pass)
- `commutator` of multi-mode Hamiltonians — should improve (intermediate sums stay compact)
- Simple `a + b` — may show slight overhead (Dict vs Vector construction)

---

## Risk Mitigation

1. **Mutable key aliasing**: The worklist in `_simplify_qmul` uses `copy(ops)` before mutation. This pattern must be preserved — never pass a Dict key reference to code that mutates vectors. Add a comment at `_simplify_qmul` documenting this.

2. **Empty QAdd**: Ensure all code paths handle `isempty(q.dict)` correctly. The zero representation is an empty Dict, not a Dict with a single zero entry.

3. **Indexed sums**: The `indices` and `non_equal` fields are orthogonal to the Dict change. Verify with `test/indexing_test.jl` that all indexed-sum operations still work.

4. **Printing determinism**: `sorted_terms()` must produce identical output across Julia sessions. The sort key `_ops_lt` uses `_sort_key` which returns tuples of `(Int, Int, Symbol)` — these have a well-defined total order.
