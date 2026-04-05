# Average System: Symbolic Operator Averaging

## Problem

QuantumCumulants.jl derives mean-field equations by taking quantum-mechanical averages of operator expressions: `⟨a†a⟩`, `⟨σ^{21}_i a⟩`, etc. The averaging step converts non-commutative operator expressions (`QField`) into commutative symbolic scalars (`BasicSymbolic`) that can participate in Symbolics/ModelingToolkit ODE systems.

The legacy implementation required three separate indexed-average types (`IndexedAverageSum`, `IndexedAverageDoubleSum`, `SpecialIndexedAverage`) totalling ~500 lines, because the old design had separate `SingleSum`/`DoubleSum` types. In the redesign, `QAdd` unifies regular sums and symbolic sums via its `indices`/`non_equal` fields, so all indexed-average types collapse into metadata on SymbolicUtils nodes.

## Goal

Implement `average()` and `undo_average()` in the redesign-v2 codebase. Averages are `BasicSymbolic{SymReal}` `Term` nodes whose `symtype` is `AvgSym`. Summation metadata from indexed `QAdd` is preserved via SymbolicUtils metadata on the resulting expression. No new indexed-average types.

## SymbolicUtils v4 Constraint

SymbolicUtils v4 uses Moshi's `@data` ADT: `BasicSymbolicImpl{T <: SymVariant}`. The type parameter `T` is the **variant** (`SymReal`, `SafeReal`, `TreeReal`) and controls algebraic behavior. It is **not** extensible — custom types cannot be used as `T`.

Instead, each `BasicSymbolic` node has a `type` field (accessed via `symtype()`) that stores the Julia type the expression represents. This is where `AvgSym` lives:

- **Variant** (`T`): always `SymReal` (default algebraic behavior)
- **Type** (`symtype`): `AvgSym` for average `Term` nodes, `Number` for sums/products of averages
- **Detection**: `symtype(x) === AvgSym` identifies average nodes (not `x isa BasicSymbolic{AvgSym}`)

## Design Principles

1. **Averages are SymbolicUtils scalars.** `average(op)` returns a `BasicSymbolic{SymReal}` `Term` node with `symtype = AvgSym`, which composes with `Num`/`CNum` via standard Symbolics arithmetic.
2. **Linearity.** `average` distributes over addition and pulls c-number prefactors out of products.
3. **No eager adjoint folding.** `average(op')` creates `avg(op')`, not `conj(avg(op))`. QC tracks `⟨a⟩` and `⟨a†⟩` as separate ODE variables.
4. **Summation metadata, not types.** Indexed `QAdd` summation info (`indices`, `non_equal`) is stored as SymbolicUtils metadata on the averaged result — no `IndexedAverageSum` wrapper types.
5. **Round-trip.** `undo_average(average(expr))` recovers the original operator expression, including summation indices.

---

## Types

### `AvgSym`

Marker type for the `symtype` of averaged expressions:

```julia
struct AvgSym end
```

Not a subtype of anything special — it's just a tag used in the `type` field of `BasicSymbolic` nodes. Detected via `symtype(x) === AvgSym`.

### Type alias

```julia
const Average = BasicSymbolic{SymReal}  # Note: same Julia type as other SymReal expressions;
                                         # distinguish via symtype(x) === AvgSym
```

### Helper predicate

```julia
is_average(x) = x isa BasicSymbolic && SymbolicUtils.iscall(x) && symtype(x) === AvgSym
```

### Symbolic function

```julia
const sym_average = SymbolicUtils.Sym{SymReal}(:avg;
    type = SymbolicUtils.FnType{Tuple{QField}, AvgSym, Nothing})
```

A symbolic function object. When called as `sym_average(op)`, SymbolicUtils creates `Term{SymReal}(sym_average, [op]; type = AvgSym)` via the built-in `(f::BasicSymbolic)(args...)` method.

### Metadata keys

```julia
struct SumIndices end     # key → Vector{Index}
struct SumNonEqual end    # key → Vector{Tuple{Index,Index}}
```

Attached to SymbolicUtils nodes via `setmetadata` when averaging an indexed `QAdd`.

---

## Functions

### `_average(op)` (internal)

Direct construction of a single average node:

```julia
function _average(op::QField)
    return sym_average(op)
end
```

Calls the symbolic function, which creates `Term{SymReal}(sym_average, [op]; type = AvgSym)`.

### `average(expr)` (public)

Linearly maps operator expressions to symbolic averages:

```julia
# Leaf operator: wrap in avg(...)
average(op::QSym) = _average(op)

# Product: pull out c-number prefactor
function average(op::QMul)
    isempty(op.args_nc) && return op.arg_c   # pure scalar
    avg = _average(QMul(1, op.args_nc))
    # CNum = Complex{Num}; multiply real/imag parts separately to stay
    # in Symbolics arithmetic (Complex{Num} * BasicSymbolic may not dispatch)
    r, i = real(op.arg_c), imag(op.arg_c)
    iszero(i) && return r * avg
    iszero(r) && return im * i * avg
    return (r + im * i) * avg
end

# Sum: distribute linearly, preserve summation metadata
function average(op::QAdd)
    result = sum(average(t) for t in op.arguments)
    # Attach summation metadata if this QAdd carries indices
    if !isempty(op.indices)
        result = SymbolicUtils.setmetadata(result, SumIndices, op.indices)
        result = SymbolicUtils.setmetadata(result, SumNonEqual, op.non_equal)
    end
    return result
end

# Scalar passthrough
average(x::Number) = x
```

### `undo_average(expr)` (public)

Recursively strips `avg(...)` to recover operator expressions:

```julia
function undo_average(x)
    if x isa BasicSymbolic && SymbolicUtils.iscall(x)
        f = SymbolicUtils.operation(x)
        if isequal(f, sym_average)
            return SymbolicUtils.arguments(x)[1]
        else
            args = map(undo_average, SymbolicUtils.arguments(x))
            result = f(args...)
            # Recover summation structure if the original node carried index metadata.
            if SymbolicUtils.hasmetadata(x, SumIndices)
                indices = SymbolicUtils.getmetadata(x, SumIndices)
                non_equal = SymbolicUtils.getmetadata(x, SumNonEqual)
                if result isa QAdd
                    result = QAdd(result.arguments, indices, non_equal)
                elseif result isa QMul
                    result = QAdd(QMul[result], indices, non_equal)
                elseif result isa QSym
                    result = QAdd(QMul[_to_qmul(result)], indices, non_equal)
                end
            end
            return result
        end
    else
        return x
    end
end
```

### `undo_average` for `Symbolics.Equation`

```julia
function undo_average(eq::Symbolics.Equation)
    return Symbolics.Equation(undo_average(eq.lhs), undo_average(eq.rhs))
end
```

### Metadata extraction helpers

For QC to recover summation info from averaged expressions:

```julia
function has_sum_metadata(x::BasicSymbolic)
    return SymbolicUtils.hasmetadata(x, SumIndices)
end

function get_sum_indices(x::BasicSymbolic)
    return SymbolicUtils.getmetadata(x, SumIndices)
end

function get_sum_non_equal(x::BasicSymbolic)
    return SymbolicUtils.getmetadata(x, SumNonEqual)
end
```

---

## Type promotion

No `Add`/`Mul` type redirects needed. In SymbolicUtils v4, `Add{T}` and `Mul{T}` are parameterized by the variant `T <: SymVariant`, not by `symtype`. Since all our expressions use `SymReal`, sums and products of averages naturally create `AddMul{SymReal}` nodes. The `symtype` of these combined expressions will be computed by `promote_symtype` (typically `Number`).

### `promote_symtype` for `sym_average`

This is handled automatically by SymbolicUtils' `FnType` machinery: when `sym_average(op)` is called, `promote_symtype` extracts `Y = AvgSym` from `FnType{Tuple{QField}, AvgSym, Nothing}` and sets it as the `type` of the resulting `Term`.

---

## `acts_on` helper

Extracts Hilbert space indices from averaged expressions (needed by QC for space-aware operations):

```julia
function acts_on(s::BasicSymbolic)
    if SymbolicUtils.iscall(s)
        f = SymbolicUtils.operation(s)
        if isequal(f, sym_average)
            return acts_on(SymbolicUtils.arguments(s)[1])
        else
            aon = Int[]
            for arg in SymbolicUtils.arguments(s)
                append!(aon, acts_on(arg))
            end
            unique!(aon)
            sort!(aon)
            return aon
        end
    else
        return Int[]
    end
end

acts_on(op::QSym) = [op.space_index]
function acts_on(op::QMul)
    aon = [x.space_index for x in op.args_nc]
    unique!(aon)
    sort!(aon)
    return aon
end
function acts_on(op::QAdd)
    aon = Int[]
    for t in op.arguments
        append!(aon, acts_on(t))
    end
    unique!(aon)
    sort!(aon)
    return aon
end
```

---

## Printing

Since averages are `BasicSymbolic{SymReal}` (same type as all other SymReal expressions), we cannot dispatch `show` on the type parameter. Instead, check `symtype` or `operation`:

```julia
# Hook into SymbolicUtils printing by checking if it's an average Term
function Base.show(io::IO, x::BasicSymbolic{SymReal})
    if SymbolicUtils.iscall(x) && isequal(SymbolicUtils.operation(x), sym_average)
        print(io, "⟨")
        show(io, SymbolicUtils.arguments(x)[1])
        print(io, "⟩")
    else
        # Fallback to default SymbolicUtils printing
        invoke(Base.show, Tuple{IO, BasicSymbolic}, io, x)
    end
end
```

Note: this overrides `show` for ALL `BasicSymbolic{SymReal}`, so the fallback must call the parent method. During implementation, verify the exact method signature needed and test that non-average SymReal expressions still print correctly.

### LaTeX recipe

```julia
@latexrecipe function f(x::BasicSymbolic{SymReal})
    if SymbolicUtils.iscall(x) && isequal(SymbolicUtils.operation(x), sym_average)
        op = SymbolicUtils.arguments(x)[1]
        return Expr(:latexifymerge, "\\langle ", op, " \\rangle")
    end
end
```

Same caveat: this recipe applies to all `BasicSymbolic{SymReal}`. Verify it doesn't interfere with existing Symbolics LaTeX rendering.

---

## Integration points

### Module (`src/SecondQuantizedAlgebra.jl`)

- Add `include("average.jl")` after `normal_order.jl` (averages depend on all operator types)
- Export: `average`, `undo_average`, `acts_on`, `has_sum_metadata`, `get_sum_indices`, `get_sum_non_equal`

### Existing code — no changes needed

- `QMul`, `QAdd`, `QSym` are unchanged
- `simplify`, `normal_order`, `expand_sums` are unchanged
- `numeric.jl` already has `numeric_average` which is independent

---

## Out of scope (QC's responsibility)

These live in QuantumCumulants.jl, not SQA:

- **Cumulant expansion** — QC walks the averaged expression tree, finds `avg(op)` nodes, expands products of operators into products of lower-order averages
- **`insert_index`** — QC substitutes concrete integer values for indices in averaged expressions
- **`NumberedOperator`** — concrete integer-indexed operators for numeric evaluation
- **`scale()`** — QC calls `expand_sums()` then does post-processing
- **`value_map`** — QC substitutes concrete values for numeric solving

## Eliminated by the redesign

These legacy types are **permanently removed**, not deferred. The redesign makes them unnecessary:

- **`IndexedAverageSum`** — replaced by SymbolicUtils metadata (`SumIndices`, `SumNonEqual`) on the averaged expression. The unified `QAdd` with `indices`/`non_equal` fields means averaging an indexed sum is just `average(QAdd)` with metadata preservation.
- **`IndexedAverageDoubleSum`** — same reason. A double sum is a `QAdd` with two indices; its average carries two-index metadata.
- **`SpecialIndexedAverage`** — was needed because the legacy had separate sum types with separate inequality constraints. The redesign's `non_equal` field on `QAdd` handles this uniformly.

---

## Examples

### Basic averaging

```julia
h = FockSpace(:cavity)
@qnumbers a::Destroy(h)

avg_a = average(a)          # ⟨a⟩ — BasicSymbolic{SymReal} with symtype = AvgSym
avg_ada = average(a' * a)   # ⟨a†a⟩

# Prefactor extraction
average(3 * a' * a)         # 3⟨a†a⟩

# Linearity over addition
average(a + a')             # ⟨a⟩ + ⟨a†⟩

# Round-trip
undo_average(avg_a) == a    # true

# Detection
is_average(avg_a)           # true
symtype(avg_a) === AvgSym   # true
```

### Indexed averaging

```julia
ha = NLevelSpace(:atom, 2)
h = FockSpace(:cavity) ⊗ ha
@qnumbers a::Destroy(h)
@variables N
i = Index(h, :i, N, ha)
σ(α, β) = IndexedOperator(Transition(h, :σ, α, β, 2), i)

# Sum of operator products
H_int = Σ(a' * σ(1, 2), i)    # QAdd with indices=[i]

avg_H = average(H_int)         # sum of averages with SumIndices metadata

# Recover summation info
has_sum_metadata(avg_H)         # true
get_sum_indices(avg_H)          # [i]
```

---

## Test plan (`test/average_test.jl`)

1. **Basic averaging**: `average(::QSym)` returns `BasicSymbolic` with `symtype === AvgSym`, `average(::QMul)` extracts prefactor, `average(::QAdd)` distributes
2. **Scalar passthrough**: `average(3) === 3`, `average(Num(x)) === Num(x)`
3. **Round-trip**: `undo_average(average(op))` recovers original for QSym, QMul, QAdd
4. **Indexed QAdd**: metadata preserved through `average()`, recoverable via `get_sum_indices`
5. **Arithmetic**: `average(a) + average(b)` produces valid SymbolicUtils expression, `average(a) * cnumber` works
6. **No double averaging**: `average(average(a))` — passthrough (averages are `BasicSymbolic` scalars, handled by `average(::Number)`)
7. **Adjoint**: `average(a')` creates `avg(a†)`, distinct from `average(a)`
8. **Multi-space**: operators on different spaces average independently
9. **Printing**: `⟨a⟩` unicode display, LaTeX `\langle a \rangle`
10. **`is_average` / `symtype`**: detection predicate works correctly
11. **Type stability**: `@inferred` on core paths where possible
