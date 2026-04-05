# Indexed Operators: Ground-Up Redesign

## Problem

Many-body quantum systems have N identical subsystems. Users write `Σ_i g_i σ_i^{21} a` symbolically. The legacy implementation required ~5000 lines of parallel code in QC because indexed types didn't fit the standard algebra.

## Goal

Redesign SQA's type system with indexing built in from the ground up. Every operator carries an `Index` field. `QMul`/`QAdd` use `Num` prefactors (no type parameter). `QAdd` doubles as a symbolic sum when it has summation indices. The user API matches the legacy interface.

## Design Principles

1. **Every operator is indexed.** `index::Index` field on every `QSym`, sentinel `NO_INDEX` when not indexed.
2. **No wrapper types.** `IndexedOperator(op, i)` is a convenience function returning the same struct with `index` set — not a new type.
3. **Unified QAdd.** Empty `indices` = regular sum. Non-empty = `Σ_i`.
4. **`Num` everywhere.** Drop parametric `QMul{T}`/`QAdd{T}`. `Symbolics.Num` as prefactor type.
5. **Symbolics.jl for indexed variables.** `IndexedVariable(:g, i)` creates a `Num` via Symbolics, not a custom type.
6. **Lazy sums, explicit expansion.** No diagonal splitting on construction. `expand_sums` called explicitly.

---

## User-facing API

Matches legacy interface:

```julia
# 1. Spaces and operators
hc = FockSpace(:cavity)
ha = NLevelSpace(:atom, 2, 1)
h = hc ⊗ ha
@qnumbers a::Destroy(h)

# 2. Indices
@variables N Δ κ Γ R ν                 # symbolic parameters (Num from Symbolics)
i = Index(h, :i, N, ha)               # index i ranges 1:N on atom space
j = Index(h, :j, N, ha)

# 3. Indexed operators (same struct, index field set)
σ(α, β, i) = IndexedOperator(Transition(h, :σ, α, β), i)

# 4. Indexed variables (returns Num via Symbolics)
g(i) = IndexedVariable(:g, i)
Γ(i, j) = DoubleIndexedVariable(:Γ, i, j)

# 5. Symbolic sums
H = -Δ*a'a + Σ(g(i)*(a'*σ(1,2,i) + a*σ(2,1,i)), i)

# 6. Derive equations (QC)
ops = [a'*a, σ(2,2,j)]
eqs = meanfield(ops, H, J; rates=rates, order=2)
```

---

## Types

### `Index`

```julia
struct Index
    name::Symbol
    range::Num          # upper bound (Int or symbolic via Num)
    space_index::Int    # which space in ProductSpace
    sym::Num            # Symbolics representation (@syms i::Int wrapped in Num)
end

const NO_INDEX = Index(:_, Num(0), 0, Num(0))
has_index(idx::Index) = idx.space_index != 0  # cheap Int check, no Num comparison
```

The `sym` field is the Symbolics symbolic integer used in prefactor expressions. Created automatically:

```julia
function Index(h::HilbertSpace, name::Symbol, range, space::HilbertSpace)
    si = _find_space_index(h, space)
    sym_var = (@syms $name::Int)[1]
    return Index(name, Num(range), si, Num(sym_var))
end
function Index(h::HilbertSpace, name::Symbol, range, si::Int)
    sym_var = (@syms $name::Int)[1]
    return Index(name, Num(range), si, Num(sym_var))
end
```

`change_index` uses `Symbolics.substitute` on prefactors via the `sym` field.

### QSym subtypes — all have `index::Index`

Every operator struct has `index::Index` as its last field:

```julia
struct Destroy <: QSym
    name::Symbol
    space_index::Int
    copy_index::Int
    index::Index
end
Destroy(name::Symbol, si::Int, ci::Int) = Destroy(name, si, ci, NO_INDEX)
Destroy(name::Symbol, si::Int) = Destroy(name, si, 1, NO_INDEX)

struct Create <: QSym
    name::Symbol
    space_index::Int
    copy_index::Int
    index::Index
end
Create(name::Symbol, si::Int, ci::Int) = Create(name, si, ci, NO_INDEX)
Create(name::Symbol, si::Int) = Create(name, si, 1, NO_INDEX)

struct Transition <: QSym
    name::Symbol
    i::Int
    j::Int
    space_index::Int
    copy_index::Int
    index::Index
end
Transition(name::Symbol, i::Int, j::Int, si::Int, ci::Int) = Transition(name, i, j, si, ci, NO_INDEX)
Transition(name::Symbol, i::Int, j::Int, si::Int) = Transition(name, i, j, si, 1, NO_INDEX)

struct Pauli <: QSym
    name::Symbol
    axis::Int
    space_index::Int
    copy_index::Int
    index::Index
end
Pauli(name::Symbol, ax::Int, si::Int, ci::Int) = Pauli(name, ax, si, ci, NO_INDEX)
Pauli(name::Symbol, ax::Int, si::Int) = Pauli(name, ax, si, 1, NO_INDEX)

struct Spin <: QSym
    name::Symbol
    axis::Int
    space_index::Int
    copy_index::Int
    index::Index
end
Spin(name::Symbol, ax::Int, si::Int, ci::Int) = Spin(name, ax, si, ci, NO_INDEX)
Spin(name::Symbol, ax::Int, si::Int) = Spin(name, ax, si, 1, NO_INDEX)

struct Position <: QSym
    name::Symbol
    space_index::Int
    copy_index::Int
    index::Index
end
Position(name::Symbol, si::Int, ci::Int) = Position(name, si, ci, NO_INDEX)
Position(name::Symbol, si::Int) = Position(name, si, 1, NO_INDEX)

struct Momentum <: QSym
    name::Symbol
    space_index::Int
    copy_index::Int
    index::Index
end
Momentum(name::Symbol, si::Int, ci::Int) = Momentum(name, si, ci, NO_INDEX)
Momentum(name::Symbol, si::Int) = Momentum(name, si, 1, NO_INDEX)
```

All existing HilbertSpace convenience constructors remain unchanged — they call the inner constructors which default `index = NO_INDEX`.

Equality, hashing, adjoint all include `index`.

### `IndexedOperator` — convenience function, NOT a type

```julia
IndexedOperator(op::Destroy, i::Index) = Destroy(op.name, op.space_index, op.copy_index, i)
IndexedOperator(op::Create, i::Index) = Create(op.name, op.space_index, op.copy_index, i)
IndexedOperator(op::Transition, i::Index) = Transition(op.name, op.i, op.j, op.space_index, op.copy_index, i)
IndexedOperator(op::Pauli, i::Index) = Pauli(op.name, op.axis, op.space_index, op.copy_index, i)
IndexedOperator(op::Spin, i::Index) = Spin(op.name, op.axis, op.space_index, op.copy_index, i)
IndexedOperator(op::Position, i::Index) = Position(op.name, op.space_index, op.copy_index, i)
IndexedOperator(op::Momentum, i::Index) = Momentum(op.name, op.space_index, op.copy_index, i)
```

Returns the same struct type with `index` field set. No wrapper.

### `IndexedVariable` — convenience function returning `Num`

```julia
function IndexedVariable(name::Symbol, i::Index)
    f = SymbolicUtils.Sym{SymbolicUtils.SymReal}(name;
        type=SymbolicUtils.FnType{Tuple{Int}, Real, Nothing},
        shape=UnitRange{Int}[])
    return Num(f(Symbolics.unwrap(i.sym)))
end
```

Creates a symbolic function `g` and calls it with the index's Symbolics symbol. `IndexedVariable(:g, i)` returns a `Num` displaying as `g(i)`. `Symbolics.substitute(expr, unwrap(i.sym) => 3)` gives `g(3)`.

### `DoubleIndexedVariable` — convenience function returning `Num`

```julia
function DoubleIndexedVariable(name::Symbol, i::Index, j::Index; identical::Bool=true)
    if !identical && i == j
        return Num(0)
    end
    f = SymbolicUtils.Sym{SymbolicUtils.SymReal}(name;
        type=SymbolicUtils.FnType{Tuple{Int, Int}, Real, Nothing},
        shape=UnitRange{Int}[])
    return Num(f(Symbolics.unwrap(i.sym), Symbolics.unwrap(j.sym)))
end
```

Returns `Num` displaying as `Γ(i, j)`.

### `QMul` — drops type parameter

```julia
struct QMul <: QTerm
    arg_c::Num
    args_nc::Vector{QSym}
end
```

Construction auto-wraps to `Num`:
```julia
QMul(c::Number, ops::Vector{QSym}) = QMul(Num(c), ops)
```

All arithmetic, equality, hashing work via `Num` operations.

### `QAdd` — drops type parameter, gains summation

```julia
struct QAdd <: QTerm
    arguments::Vector{QMul}
    indices::Vector{Index}
    non_equal::Vector{Tuple{Index,Index}}
end
QAdd(args::Vector{QMul}) = QAdd(args, Index[], Tuple{Index,Index}[])
```

- `indices = []` → regular sum (backward compatible)
- `indices = [i]` → `Σ_i`
- `indices = [i, j]` → `Σ_i Σ_j`
- `non_equal = [(i,j)]` → `i ≠ j` constraint

### `Σ` — sum constructor

```julia
function Σ(expr, i::Index, non_equal::Vector{Index}=Index[])
    qadd = _to_qadd(expr)  # convert expr to QAdd if needed
    ne_pairs = [(i, j) for j in non_equal]
    return QAdd(qadd.arguments, [i], ne_pairs)
end
function Σ(expr, i::Index, j::Index, args...)
    # Multi-index sum: accumulate all indices on one QAdd
    all_indices = [i, j, args...]
    qadd = _to_qadd(expr)
    return QAdd(qadd.arguments, collect(Index, all_indices), Tuple{Index,Index}[])
end
const ∑ = Σ
```

Smart constructor: if `i` doesn't appear in any term, returns `(range - |non_equal|) * expr`.

---

## Simplification

### `_same_site` predicate

Replaces all scattered `a.space_index == b.space_index && a.copy_index == b.copy_index` checks:

```julia
function _same_site(a::QSym, b::QSym)
    a.space_index == b.space_index &&
    a.copy_index == b.copy_index &&
    a.index == b.index
end
```

Two operators interact (apply commutation/composition rules) only when on the same site: same space, same copy, same index.

Used in `_simplify_product!` for all rules (Transition composition, Pauli product, Fock swap, Spin swap, PhaseSpace swap).

---

## Commutator

### Same index or both NO_INDEX → normal rules

```julia
commutator(σ_i^{12}, σ_i^{21})  # same index → [σ^{12}, σ^{21}] = σ^{11} - σ^{22}
commutator(a, a')                # both NO_INDEX → [a, a†] = 1
```

### Different index, same space → zero

```julia
commutator(σ_i, σ_j)  # different index, same space → 0 (independent copies)
```

### QAdd with indices (sum) with indexed operator → diagonal collapse

```julia
[Σ_i f(i), B_j]  where i,j on same space
= [f(j), B_j]    (only i=j term survives, off-diagonal commutes to 0)
```

Implementation:
```julia
function commutator(a::QAdd, b::QSym)
    if !isempty(a.indices) && has_index(b.index)
        for idx in a.indices
            if idx.space_index == b.index.space_index
                # Diagonal collapse: substitute sum index → operator's index
                collapsed = change_index(a, idx, b.index)
                # Remove this index from summation (it's been substituted)
                return commutator(collapsed, b)
            end
        end
    end
    # Regular distribution over terms
    ...
end
```

### QAdd(sum) with QAdd(sum) → compose collapses

```julia
[Σ_i f(i), Σ_j g(j)] = Σ_j [f(j), g(j)]
```

---

## `change_index`

```julia
change_index(expr, from::Index, to::Index)
```

Recursively substitutes throughout the expression tree:

| Input | Action |
|---|---|
| `QSym` with `index == from` | Return copy with `index = to` |
| `QSym` with `index != from` | Unchanged |
| `QMul(c, ops)` | `QMul(substitute(c, from.sym => to.sym), [change_index(op, from, to) for op in ops])` |
| `QAdd(terms, indices, ne)` | Apply to terms, replace `from` in indices and non_equal |
| Scalar / `Num` | `Symbolics.substitute(expr, Dict(unwrap(from.sym) => unwrap(to.sym)))` |

The key: `change_index` on prefactors delegates to `Symbolics.substitute` using the `sym` fields. This is why `Index` stores a Symbolics symbol — so `g(i)` automatically becomes `g(j)` when changing `i → j`.

---

## `expand_sums`

```julia
expand_sums(expr) → expanded expression with diagonal splitting
```

Called explicitly by QC's `scale()`. Applies the splitting rules:

**Core rule:**
```
Σ_i(A_i * B_j) where i,j same space, j not in non_equal
→ Σ_{i≠j}(A_i * B_j) + A_j * B_j
```

All 15 legacy splitting rules live in this one function. Not called during construction or arithmetic.

---

## `get_indices`

```julia
get_indices(expr) → Vector{Index}
```

Collects all non-`NO_INDEX` indices in an expression. Used by QC to detect indexed systems.

---

## Printing

| Expression | Display |
|---|---|
| `Destroy(:a, 1, 1, i)` | `aᵢ` |
| `Transition(:σ, 1, 2, 2, 1, i)` | `σᵢ₁₂` |
| `IndexedVariable(:g, i)` | `gᵢ` (via Symbolics display) |
| `QAdd(terms, [i], [])` | `Σ(i=1:N) terms` |
| `QAdd(terms, [i], [(i,j)])` | `Σ(i=1:N,i≠j) terms` |

---

## Integration with QC

| QC step | What happens |
|---|---|
| `commutator(im*H, op)` | H is QAdd with indices. Distributes over terms. Same-space indexed ops trigger diagonal collapse. |
| `average(QAdd_with_indices)` | QC wraps: one new method `average(::QAdd)` that handles sums. |
| `cumulant_expansion` | Works on averaged terms inside sum. No special-casing. |
| `scale(me)` | Calls `expand_sums` → converts symbolic sums to N-dependent factors. |
| `complete(me)` | Expression tree traversal works unchanged. |

---

## File layout

| File | Contents |
|---|---|
| New: `src/index.jl` | `Index`, `NO_INDEX`, `has_index`, constructors, equality, hashing |
| Modify: `src/fock.jl` | Add `index::Index` field, `IndexedOperator` convenience |
| Modify: `src/nlevel.jl` | Same |
| Modify: `src/pauli.jl` | Same |
| Modify: `src/spin.jl` | Same |
| Modify: `src/phase_space.jl` | Same |
| Modify: `src/qmul.jl` | Drop `{T}`, use `Num` prefactor |
| Modify: `src/qadd.jl` | Drop `{T}`, add `indices`/`non_equal`, `Σ` constructor |
| New: `src/indexed_variables.jl` | `IndexedVariable`, `DoubleIndexedVariable` (return `Num`) |
| New: `src/change_index.jl` | `change_index`, `get_indices` |
| New: `src/expand_sums.jl` | `expand_sums` with diagonal splitting rules |
| Modify: `src/simplify.jl` | Use `_same_site` predicate |
| Modify: `src/commutator.jl` | Index checks, diagonal collapse for sums |
| Modify: `src/printing.jl` | Index subscripts, Σ notation |
| Modify: `src/SecondQuantizedAlgebra.jl` | Includes and exports |
| New: `test/index_test.jl` | Index, indexed operators, indexed variables |
| New: `test/qsum_test.jl` | QAdd with indices, Σ, commutator with sums |
| New: `test/expand_sums_test.jl` | Diagonal splitting rules |

## Exports

`Index`, `NO_INDEX`, `has_index`, `IndexedOperator`, `IndexedVariable`, `DoubleIndexedVariable`, `Σ`, `∑`, `change_index`, `expand_sums`, `get_indices`
