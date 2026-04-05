# Indexed Operators: Ground-Up Rewrite Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Rewrite SQA's type system with indexing built into every operator, `Num` prefactors (no type parameters), and unified `QAdd` that doubles as symbolic sums.

**Architecture:** Every `QSym` subtype gains `index::Index` (sentinel `NO_INDEX`). `QMul`/`QAdd` drop their `{T}` type parameter, using `Symbolics.Num` for all prefactors. `QAdd` gains `indices`/`non_equal` fields for symbolic sums. `IndexedOperator` is a convenience function (not a type). `IndexedVariable`/`DoubleIndexedVariable` return `Num` via Symbolics. New `Σ` constructor, `change_index`, `expand_sums`, `get_indices` utilities.

**Tech Stack:** Julia, Symbolics.jl (`Num`, `@syms`, `substitute`), SymbolicUtils.jl, existing SQA infrastructure.

**Note:** Tests may break during intermediate tasks. Full test suite should pass only after Task 13.

---

### Task 1: Create `src/index.jl` — Index type and sentinel

**Files:**
- Create: `src/index.jl`

- [ ] **Step 1: Create `src/index.jl`**

```julia
using Symbolics: Symbolics, Num, @syms
using SymbolicUtils: SymbolicUtils

"""
    Index

Symbolic summation index for many-body systems.

Fields:
- `name::Symbol` — display name
- `range::Num` — upper bound (Int or symbolic)
- `space_index::Int` — which space in ProductSpace
- `sym::Num` — Symbolics symbolic integer for variable substitution
"""
struct Index
    name::Symbol
    range::Num
    space_index::Int
    sym::Num
end

const NO_INDEX = Index(:_, Num(0), 0, Num(0))
has_index(idx::Index) = idx.space_index != 0

function Base.:(==)(a::Index, b::Index)
    return a.name == b.name && isequal(a.range, b.range) && a.space_index == b.space_index
end
function Base.hash(a::Index, h::UInt)
    return hash(:Index, hash(a.name, hash(a.space_index, h)))
end

# Construction from HilbertSpace
function Index(h::HilbertSpace, name::Symbol, range, space::HilbertSpace)
    si = _find_space_index(h, space)
    sym_var = Num(SymbolicUtils.Sym{SymbolicUtils.SymReal}(name;
        type = SymbolicUtils.FnType{Tuple{}, Int, Nothing},
        shape = UnitRange{Int}[]))
    return Index(name, Num(range), si, sym_var)
end
function Index(h::HilbertSpace, name::Symbol, range, si::Int)
    sym_var = Num(SymbolicUtils.Sym{SymbolicUtils.SymReal}(name;
        type = SymbolicUtils.FnType{Tuple{}, Int, Nothing},
        shape = UnitRange{Int}[]))
    return Index(name, Num(range), si, sym_var)
end

# Find space index in ProductSpace
_find_space_index(h::HilbertSpace, space::HilbertSpace) = 1
function _find_space_index(h::ProductSpace, space::HilbertSpace)
    for (i, s) in enumerate(h.spaces)
        actual = _unwrap_space(s)
        actual == space && return i
    end
    throw(ArgumentError("Space $space not found in $h"))
end
```

### Task 2: Add `index::Index` to all QSym subtypes

**Files:**
- Modify: `src/fock.jl`
- Modify: `src/nlevel.jl`
- Modify: `src/pauli.jl`
- Modify: `src/spin.jl`
- Modify: `src/phase_space.jl`

- [ ] **Step 1: Rewrite `src/fock.jl`**

```julia
"""
    Destroy <: QSym

Bosonic annihilation operator.
"""
struct Destroy <: QSym
    name::Symbol
    space_index::Int
    copy_index::Int
    index::Index
end
Destroy(name::Symbol, si::Int, ci::Int) = Destroy(name, si, ci, NO_INDEX)
Destroy(name::Symbol, si::Int) = Destroy(name, si, 1, NO_INDEX)

"""
    Create <: QSym

Bosonic creation operator.
"""
struct Create <: QSym
    name::Symbol
    space_index::Int
    copy_index::Int
    index::Index
end
Create(name::Symbol, si::Int, ci::Int) = Create(name, si, ci, NO_INDEX)
Create(name::Symbol, si::Int) = Create(name, si, 1, NO_INDEX)

# Construction from Hilbert spaces
Destroy(h::FockSpace, name::Symbol) = Destroy(name, 1)
Create(h::FockSpace, name::Symbol) = Create(name, 1)

function Destroy(h::ProductSpace, name::Symbol, idx::Int)
    1 <= idx <= length(h.spaces) || throw(ArgumentError("Index $idx out of range for ProductSpace with $(length(h.spaces)) spaces"))
    _unwrap_space(h.spaces[idx]) isa FockSpace || throw(ArgumentError("Space at index $idx is not a FockSpace"))
    return Destroy(name, idx)
end
function Create(h::ProductSpace, name::Symbol, idx::Int)
    1 <= idx <= length(h.spaces) || throw(ArgumentError("Index $idx out of range for ProductSpace with $(length(h.spaces)) spaces"))
    _unwrap_space(h.spaces[idx]) isa FockSpace || throw(ArgumentError("Space at index $idx is not a FockSpace"))
    return Create(name, idx)
end

# IndexedOperator convenience
IndexedOperator(op::Destroy, i::Index) = Destroy(op.name, op.space_index, op.copy_index, i)
IndexedOperator(op::Create, i::Index) = Create(op.name, op.space_index, op.copy_index, i)

# Adjoint
Base.adjoint(op::Destroy) = Create(op.name, op.space_index, op.copy_index, op.index)
Base.adjoint(op::Create) = Destroy(op.name, op.space_index, op.copy_index, op.index)

# Equality
Base.isequal(a::Destroy, b::Destroy) = a.name == b.name && a.space_index == b.space_index && a.copy_index == b.copy_index && a.index == b.index
Base.isequal(a::Create, b::Create) = a.name == b.name && a.space_index == b.space_index && a.copy_index == b.copy_index && a.index == b.index
Base.:(==)(a::Destroy, b::Destroy) = isequal(a, b)
Base.:(==)(a::Create, b::Create) = isequal(a, b)

# Hashing
Base.hash(a::Destroy, h::UInt) = hash(:Destroy, hash(a.name, hash(a.space_index, hash(a.copy_index, hash(a.index, h)))))
Base.hash(a::Create, h::UInt) = hash(:Create, hash(a.name, hash(a.space_index, hash(a.copy_index, hash(a.index, h)))))

# Canonical ordering
ladder(::Create) = 0
ladder(::Destroy) = 1
```

- [ ] **Step 2: Rewrite `src/nlevel.jl`**

Same pattern: add `index::Index` as last field to `Transition`, add convenience constructors defaulting to `NO_INDEX`, add `IndexedOperator(op::Transition, i::Index)`, update equality/hash/adjoint to include `index`.

```julia
"""
    NLevelSpace <: HilbertSpace

Hilbert space for N-level systems (atoms, qubits, etc.).
"""
struct NLevelSpace <: HilbertSpace
    name::Symbol
    n::Int
    ground_state::Int
    function NLevelSpace(name::Symbol, n::Int, ground_state::Int)
        1 <= ground_state <= n || throw(ArgumentError("Ground state $ground_state out of range 1:$n"))
        return new(name, n, ground_state)
    end
end
Base.:(==)(a::NLevelSpace, b::NLevelSpace) = a.name == b.name && a.n == b.n && a.ground_state == b.ground_state
Base.hash(a::NLevelSpace, h::UInt) = hash(:NLevelSpace, hash(a.name, hash(a.n, hash(a.ground_state, h))))

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

function Transition(h::NLevelSpace, name::Symbol, i::Int, j::Int)
    1 <= i <= h.n || throw(ArgumentError("Level i=$i out of range 1:$(h.n)"))
    1 <= j <= h.n || throw(ArgumentError("Level j=$j out of range 1:$(h.n)"))
    return Transition(name, i, j, 1)
end
function Transition(h::ProductSpace, name::Symbol, i::Int, j::Int, idx::Int)
    1 <= idx <= length(h.spaces) || throw(ArgumentError("Index $idx out of range"))
    space = _unwrap_space(h.spaces[idx])
    space isa NLevelSpace || throw(ArgumentError("Space at index $idx is not an NLevelSpace"))
    1 <= i <= space.n || throw(ArgumentError("Level i=$i out of range 1:$(space.n)"))
    1 <= j <= space.n || throw(ArgumentError("Level j=$j out of range 1:$(space.n)"))
    return Transition(name, i, j, idx)
end

IndexedOperator(op::Transition, i::Index) = Transition(op.name, op.i, op.j, op.space_index, op.copy_index, i)

Base.adjoint(op::Transition) = Transition(op.name, op.j, op.i, op.space_index, op.copy_index, op.index)

Base.isequal(a::Transition, b::Transition) = a.name == b.name && a.i == b.i && a.j == b.j && a.space_index == b.space_index && a.copy_index == b.copy_index && a.index == b.index
Base.:(==)(a::Transition, b::Transition) = isequal(a, b)
Base.hash(a::Transition, h::UInt) = hash(:Transition, hash(a.name, hash(a.i, hash(a.j, hash(a.space_index, hash(a.copy_index, hash(a.index, h)))))))

ladder(::Transition) = 0
```

- [ ] **Step 3: Rewrite `src/pauli.jl`**

Same pattern for `Pauli`:

```julia
struct PauliSpace <: HilbertSpace
    name::Symbol
end
Base.:(==)(a::PauliSpace, b::PauliSpace) = a.name == b.name
Base.hash(a::PauliSpace, h::UInt) = hash(:PauliSpace, hash(a.name, h))

struct Pauli <: QSym
    name::Symbol
    axis::Int
    space_index::Int
    copy_index::Int
    index::Index
    function Pauli(name::Symbol, axis::Int, si::Int, ci::Int, idx::Index)
        1 <= axis <= 3 || throw(ArgumentError("Pauli axis must be 1, 2, or 3, got $axis"))
        return new(name, axis, si, ci, idx)
    end
end
Pauli(name::Symbol, axis::Int, si::Int, ci::Int) = Pauli(name, axis, si, ci, NO_INDEX)
Pauli(name::Symbol, axis::Int, si::Int) = Pauli(name, axis, si, 1, NO_INDEX)

Pauli(h::PauliSpace, name::Symbol, axis::Int) = Pauli(name, axis, 1)
function Pauli(h::ProductSpace, name::Symbol, axis::Int, idx::Int)
    1 <= idx <= length(h.spaces) || throw(ArgumentError("Index $idx out of range"))
    _unwrap_space(h.spaces[idx]) isa PauliSpace || throw(ArgumentError("Space at index $idx is not a PauliSpace"))
    return Pauli(name, axis, idx)
end

IndexedOperator(op::Pauli, i::Index) = Pauli(op.name, op.axis, op.space_index, op.copy_index, i)

Base.adjoint(op::Pauli) = op
Base.isequal(a::Pauli, b::Pauli) = a.name == b.name && a.axis == b.axis && a.space_index == b.space_index && a.copy_index == b.copy_index && a.index == b.index
Base.:(==)(a::Pauli, b::Pauli) = isequal(a, b)
Base.hash(a::Pauli, h::UInt) = hash(:Pauli, hash(a.name, hash(a.axis, hash(a.space_index, hash(a.copy_index, hash(a.index, h))))))
ladder(::Pauli) = 0
```

- [ ] **Step 4: Rewrite `src/spin.jl`**

Same pattern for `Spin`:

```julia
struct SpinSpace <: HilbertSpace
    name::Symbol
    spin::Rational{Int}
    function SpinSpace(name::Symbol, spin::Rational{Int})
        spin > 0 || throw(ArgumentError("Spin must be positive, got $spin"))
        denominator(spin) in (1, 2) || throw(ArgumentError("Spin must be integer or half-integer, got $spin"))
        return new(name, spin)
    end
    SpinSpace(name::Symbol, spin::Integer) = SpinSpace(name, spin // 1)
end
Base.:(==)(a::SpinSpace, b::SpinSpace) = a.name == b.name && a.spin == b.spin
Base.hash(a::SpinSpace, h::UInt) = hash(:SpinSpace, hash(a.name, hash(a.spin, h)))

struct Spin <: QSym
    name::Symbol
    axis::Int
    space_index::Int
    copy_index::Int
    index::Index
    function Spin(name::Symbol, axis::Int, si::Int, ci::Int, idx::Index)
        1 <= axis <= 3 || throw(ArgumentError("Spin axis must be 1, 2, or 3, got $axis"))
        return new(name, axis, si, ci, idx)
    end
end
Spin(name::Symbol, axis::Int, si::Int, ci::Int) = Spin(name, axis, si, ci, NO_INDEX)
Spin(name::Symbol, axis::Int, si::Int) = Spin(name, axis, si, 1, NO_INDEX)

Spin(h::SpinSpace, name::Symbol, axis::Int) = Spin(name, axis, 1)
function Spin(h::ProductSpace, name::Symbol, axis::Int, idx::Int)
    1 <= idx <= length(h.spaces) || throw(ArgumentError("Index $idx out of range"))
    _unwrap_space(h.spaces[idx]) isa SpinSpace || throw(ArgumentError("Space at index $idx is not a SpinSpace"))
    return Spin(name, axis, idx)
end

IndexedOperator(op::Spin, i::Index) = Spin(op.name, op.axis, op.space_index, op.copy_index, i)

Base.adjoint(op::Spin) = op
Base.isequal(a::Spin, b::Spin) = a.name == b.name && a.axis == b.axis && a.space_index == b.space_index && a.copy_index == b.copy_index && a.index == b.index
Base.:(==)(a::Spin, b::Spin) = isequal(a, b)
Base.hash(a::Spin, h::UInt) = hash(:Spin, hash(a.name, hash(a.axis, hash(a.space_index, hash(a.copy_index, hash(a.index, h))))))
ladder(::Spin) = 0
```

- [ ] **Step 5: Rewrite `src/phase_space.jl`**

Same pattern for `Position` and `Momentum`:

```julia
struct PhaseSpace <: HilbertSpace
    name::Symbol
end
Base.:(==)(a::PhaseSpace, b::PhaseSpace) = a.name == b.name
Base.hash(a::PhaseSpace, h::UInt) = hash(:PhaseSpace, hash(a.name, h))

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

Position(h::PhaseSpace, name::Symbol) = Position(name, 1)
Momentum(h::PhaseSpace, name::Symbol) = Momentum(name, 1)
function Position(h::ProductSpace, name::Symbol, idx::Int)
    1 <= idx <= length(h.spaces) || throw(ArgumentError("Index $idx out of range for ProductSpace with $(length(h.spaces)) spaces"))
    _unwrap_space(h.spaces[idx]) isa PhaseSpace || throw(ArgumentError("Space at index $idx is not a PhaseSpace"))
    return Position(name, idx)
end
function Momentum(h::ProductSpace, name::Symbol, idx::Int)
    1 <= idx <= length(h.spaces) || throw(ArgumentError("Index $idx out of range for ProductSpace with $(length(h.spaces)) spaces"))
    _unwrap_space(h.spaces[idx]) isa PhaseSpace || throw(ArgumentError("Space at index $idx is not a PhaseSpace"))
    return Momentum(name, idx)
end

IndexedOperator(op::Position, i::Index) = Position(op.name, op.space_index, op.copy_index, i)
IndexedOperator(op::Momentum, i::Index) = Momentum(op.name, op.space_index, op.copy_index, i)

Base.adjoint(op::Position) = op
Base.adjoint(op::Momentum) = op
Base.isequal(a::Position, b::Position) = a.name == b.name && a.space_index == b.space_index && a.copy_index == b.copy_index && a.index == b.index
Base.isequal(a::Momentum, b::Momentum) = a.name == b.name && a.space_index == b.space_index && a.copy_index == b.copy_index && a.index == b.index
Base.:(==)(a::Position, b::Position) = isequal(a, b)
Base.:(==)(a::Momentum, b::Momentum) = isequal(a, b)
Base.hash(a::Position, h::UInt) = hash(:Position, hash(a.name, hash(a.space_index, hash(a.copy_index, hash(a.index, h)))))
Base.hash(a::Momentum, h::UInt) = hash(:Momentum, hash(a.name, hash(a.space_index, hash(a.copy_index, hash(a.index, h)))))
ladder(::Position) = 0
ladder(::Momentum) = 1
```

### Task 3: Rewrite `QMul` to use `Num`

**Files:**
- Modify: `src/qmul.jl`

- [ ] **Step 1: Rewrite `src/qmul.jl` — drop `{T}`, use `Num`**

Key changes:
- `struct QMul <: QTerm` (no type parameter)
- `arg_c::Num` instead of `arg_c::T`
- Inner constructor wraps any `Number` to `Num`
- Remove all `promote_type`/`convert` logic
- Sort key: `op -> (op.space_index, op.copy_index)`
- All `Number` in multiplication methods gets wrapped to `Num`

The full file is large. The key structural changes:

```julia
struct QMul <: QTerm
    arg_c::Num
    args_nc::Vector{QSym}
    function QMul(arg_c::Num, args_nc::Vector{QSym})
        return new(arg_c, args_nc)
    end
end
QMul(arg_c::Number, args_nc::Vector{QSym}) = QMul(Num(arg_c), args_nc)
QMul(args_nc::Vector{QSym}) = QMul(Num(1), args_nc)
```

All arithmetic methods stay the same structurally but remove type parameters:
- `Base.:*(a::QSym, b::QSym)` → `QMul(Num(1), sort!([a, b]; by=...))`
- `Base.:*(a::QSym, b::Number)` → `QMul(Num(b), QSym[a])`
- `Base.:*(a::QMul, b::Number)` → `QMul(a.arg_c * Num(b), a.args_nc)`
- etc.

Remove `promote_rule` and `convert` methods — no longer needed.

### Task 4: Rewrite `QAdd` to use `Num` and gain summation fields

**Files:**
- Modify: `src/qadd.jl`

- [ ] **Step 1: Rewrite `src/qadd.jl`**

Key changes:
- `struct QAdd <: QTerm` (no type parameter)
- Add `indices::Vector{Index}` and `non_equal::Vector{Tuple{Index,Index}}`
- Backward-compatible constructor: `QAdd(args::Vector{QMul}) = QAdd(args, Index[], Tuple{Index,Index}[])`
- Remove all `promote_type`/`convert` logic
- All `Number` gets wrapped to `Num`

```julia
struct QAdd <: QTerm
    arguments::Vector{QMul}
    indices::Vector{Index}
    non_equal::Vector{Tuple{Index,Index}}
    function QAdd(args::Vector{QMul}, indices::Vector{Index}, non_equal::Vector{Tuple{Index,Index}})
        return new(args, indices, non_equal)
    end
end
QAdd(args::Vector{QMul}) = QAdd(args, Index[], Tuple{Index,Index}[])
```

Helper functions stay but drop type parameters:
```julia
_to_qmul(a::QSym) = QMul(Num(1), QSym[a])
_to_qmul(a::QMul) = a
_scalar_qmul(x::Number) = QMul(Num(x), QSym[])
```

All arithmetic: remove `{S}`, `{T}`, `promote_type(S, T)`, `convert(QMul{TT}, ...)`. Just use `Num` arithmetic directly.

### Task 5: Create `src/indexed_variables.jl`

**Files:**
- Create: `src/indexed_variables.jl`

- [ ] **Step 1: Create `src/indexed_variables.jl`**

```julia
using SymbolicUtils: SymbolicUtils, SymReal, FnType
using Symbolics: Symbolics, Num, unwrap

"""
    IndexedVariable(name::Symbol, i::Index) -> Num

Create a symbolic indexed variable `name(i)`.
"""
function IndexedVariable(name::Symbol, i::Index)
    f = SymbolicUtils.Sym{SymReal}(name;
        type = FnType{Tuple{Int}, Real, Nothing},
        shape = UnitRange{Int}[])
    return Num(f(unwrap(i.sym)))
end

"""
    DoubleIndexedVariable(name::Symbol, i::Index, j::Index; identical::Bool=true) -> Num

Create a symbolic double-indexed variable `name(i, j)`.
If `identical=false`, returns `Num(0)` when `i == j`.
"""
function DoubleIndexedVariable(name::Symbol, i::Index, j::Index;
    identical::Bool = true)
    if !identical && i == j
        return Num(0)
    end
    f = SymbolicUtils.Sym{SymReal}(name;
        type = FnType{Tuple{Int, Int}, Real, Nothing},
        shape = UnitRange{Int}[])
    return Num(f(unwrap(i.sym), unwrap(j.sym)))
end
```

### Task 6: Create `src/qsum.jl` — `Σ` constructor

**Files:**
- Create: `src/qsum.jl`

- [ ] **Step 1: Create `src/qsum.jl`**

```julia
"""
    Σ(expr, i::Index, non_equal::Vector{Index}=Index[])

Create a symbolic sum over index `i`. Returns a `QAdd` with `indices = [i]`.
"""
function Σ(expr::QMul, i::Index, non_equal::Vector{Index} = Index[])
    ne_pairs = Tuple{Index,Index}[(i, j) for j in non_equal]
    return QAdd(QMul[expr], [i], ne_pairs)
end
function Σ(expr::QAdd, i::Index, non_equal::Vector{Index} = Index[])
    ne_pairs = Tuple{Index,Index}[(i, j) for j in non_equal]
    # Merge: if expr already has indices, add i to them
    all_indices = vcat(expr.indices, [i])
    all_ne = vcat(expr.non_equal, ne_pairs)
    return QAdd(expr.arguments, all_indices, all_ne)
end
function Σ(expr::QSym, i::Index, non_equal::Vector{Index} = Index[])
    return Σ(_to_qmul(expr), i, non_equal)
end
function Σ(expr::Number, i::Index, non_equal::Vector{Index} = Index[])
    return Σ(_scalar_qmul(expr), i, non_equal)
end

# Multi-index: Σ(expr, i, j) = double sum
function Σ(expr, i::Index, j::Index, rest::Index...)
    inner = Σ(expr, i)
    return Σ(inner, j, rest...)
end

const ∑ = Σ
```

### Task 7: Update `src/simplify.jl` — use `_same_site`

**Files:**
- Modify: `src/simplify.jl`

- [ ] **Step 1: Add `_same_site` predicate and update all rules**

Add at the top of the file:
```julia
function _same_site(a::QSym, b::QSym)
    return a.space_index == b.space_index &&
           a.copy_index == b.copy_index &&
           a.index == b.index
end
```

Replace all adjacent-pair checks. In `_simplify_product!`:

Transition rule — change:
```julia
if a isa Transition && b isa Transition && a.space_index == b.space_index && a.copy_index == b.copy_index && a.name == b.name
```
to:
```julia
if a isa Transition && b isa Transition && _same_site(a, b) && a.name == b.name
```

And update constructor calls to preserve index:
```julia
ops[i] = Transition(a.name, a.i, b.j, a.space_index, a.copy_index, a.index)
```

Same pattern for Pauli, Fock, Spin, PhaseSpace rules — replace scattered checks with `_same_site(a, b)` and preserve `a.index` in constructed operators.

Drop `_complex_promote` — `Num` handles complex natively. Replace `Complex{Int}` type annotations with `Num`. The `_collect_like_terms` function operates on `QMul` (no type parameter) with `Num` prefactors.

Remove `SymbolicUtils.simplify` and `Symbolics.expand` overloads that use type parameters — rewrite to work with non-parametric `QMul`/`QAdd`.

### Task 8: Update `src/commutator.jl` — index checks + sum collapse

**Files:**
- Modify: `src/commutator.jl`

- [ ] **Step 1: Update short-circuits and add sum diagonal collapse**

QSym short-circuit — use `_same_site`:
```julia
function commutator(a::QSym, b::QSym)
    _same_site(a, b) || return _ZERO_QADD
    isequal(a, b) && return _ZERO_QADD
    return simplify(a * b - b * a)
end
```

Update `_has_space` to check `(space_index, copy_index, index)`:
```julia
function _has_space(ops::Vector{QSym}, si::Int, ci::Int, idx::Index)
    for op in ops
        op.space_index == si && op.copy_index == ci && op.index == idx && return true
    end
    return false
end

function _has_shared_space(a::Vector{QSym}, b::Vector{QSym})
    for x in a, y in b
        _same_site(x, y) && return true
    end
    return false
end
```

Update call sites:
```julia
function commutator(a::QMul, b::QSym)
    _has_space(a.args_nc, b.space_index, b.copy_index, b.index) || return _ZERO_QADD
    return simplify(a * b - b * a)
end
```

Update `_ZERO_QADD`:
```julia
const _ZERO_QADD = QAdd(QMul[QMul(Num(0), QSym[])])
```

Add QAdd-with-indices commutator (diagonal collapse):
```julia
function commutator(a::QAdd, b::QSym)
    if !isempty(a.indices) && has_index(b.index)
        for (k, idx) in enumerate(a.indices)
            if idx.space_index == b.index.space_index
                collapsed = change_index(a, idx, b.index)
                # Remove collapsed index from indices
                new_indices = [a.indices[j] for j in eachindex(a.indices) if j != k]
                new_ne = filter(p -> p[1] != idx && p[2] != idx, a.non_equal)
                collapsed_qadd = QAdd(collapsed.arguments, new_indices, new_ne)
                return commutator(collapsed_qadd, b)
            end
        end
    end
    # Regular distribution
    all_terms = QMul[]
    for a_ in a.arguments
        _append_terms!(all_terms, commutator(a_, b))
    end
    return QAdd(all_terms, a.indices, a.non_equal)
end
```

Update `_append_terms!` to work without type parameter:
```julia
function _append_terms!(all_terms::Vector{QMul}, result::QAdd)
    for t in result.arguments
        push!(all_terms, t)
    end
    return all_terms
end
```

### Task 9: Create `src/change_index.jl`

**Files:**
- Create: `src/change_index.jl`

- [ ] **Step 1: Create `src/change_index.jl`**

```julia
using Symbolics: Symbolics, Num, unwrap, substitute

"""
    change_index(expr, from::Index, to::Index)

Substitute index `from` with `to` throughout an expression tree.
"""
change_index(x::Number, ::Index, ::Index) = x
change_index(x::Num, from::Index, to::Index) = Num(substitute(unwrap(x), Dict(unwrap(from.sym) => unwrap(to.sym))))

function change_index(op::Destroy, from::Index, to::Index)
    idx = op.index == from ? to : op.index
    return Destroy(op.name, op.space_index, op.copy_index, idx)
end
function change_index(op::Create, from::Index, to::Index)
    idx = op.index == from ? to : op.index
    return Create(op.name, op.space_index, op.copy_index, idx)
end
function change_index(op::Transition, from::Index, to::Index)
    idx = op.index == from ? to : op.index
    return Transition(op.name, op.i, op.j, op.space_index, op.copy_index, idx)
end
function change_index(op::Pauli, from::Index, to::Index)
    idx = op.index == from ? to : op.index
    return Pauli(op.name, op.axis, op.space_index, op.copy_index, idx)
end
function change_index(op::Spin, from::Index, to::Index)
    idx = op.index == from ? to : op.index
    return Spin(op.name, op.axis, op.space_index, op.copy_index, idx)
end
function change_index(op::Position, from::Index, to::Index)
    idx = op.index == from ? to : op.index
    return Position(op.name, op.space_index, op.copy_index, idx)
end
function change_index(op::Momentum, from::Index, to::Index)
    idx = op.index == from ? to : op.index
    return Momentum(op.name, op.space_index, op.copy_index, idx)
end

function change_index(m::QMul, from::Index, to::Index)
    new_c = change_index(m.arg_c, from, to)
    new_ops = QSym[change_index(op, from, to) for op in m.args_nc]
    return QMul(new_c, new_ops)
end

function change_index(s::QAdd, from::Index, to::Index)
    new_terms = QMul[change_index(t, from, to) for t in s.arguments]
    new_indices = [idx == from ? to : idx for idx in s.indices]
    new_ne = Tuple{Index,Index}[(a == from ? to : a, b == from ? to : b) for (a, b) in s.non_equal]
    return QAdd(new_terms, new_indices, new_ne)
end

"""
    get_indices(expr) -> Vector{Index}

Collect all non-NO_INDEX indices in an expression.
"""
get_indices(::Number) = Index[]
get_indices(::Num) = Index[]
function get_indices(op::QSym)
    has_index(op.index) ? Index[op.index] : Index[]
end
function get_indices(m::QMul)
    inds = Index[]
    for op in m.args_nc
        has_index(op.index) && op.index ∉ inds && push!(inds, op.index)
    end
    return inds
end
function get_indices(s::QAdd)
    inds = copy(s.indices)
    for t in s.arguments
        for idx in get_indices(t)
            idx ∉ inds && push!(inds, idx)
        end
    end
    return inds
end
```

### Task 10: Create `src/expand_sums.jl`

**Files:**
- Create: `src/expand_sums.jl`

- [ ] **Step 1: Create `src/expand_sums.jl`**

```julia
"""
    expand_sums(expr::QAdd) -> QAdd

Explicit diagonal splitting for symbolic sums. Called by QC's `scale()`.

Core rule: Σ_i(A_i * B_j) where i,j same space, j not in non_equal
→ Σ_{i≠j}(A_i * B_j) + A_j * B_j
"""
function expand_sums(s::QAdd)
    isempty(s.indices) && return s  # Nothing to expand

    result_terms = QMul[]
    result_indices = copy(s.indices)
    result_ne = copy(s.non_equal)

    for term in s.arguments
        term_indices = get_indices(term)
        new_ne_indices = Index[]

        for idx in term_indices
            # Check if this index conflicts with a sum index
            for sum_idx in s.indices
                if idx != sum_idx &&
                    idx.space_index == sum_idx.space_index &&
                    !((sum_idx, idx) in s.non_equal) &&
                    !((idx, sum_idx) in s.non_equal)
                    # Diagonal split needed
                    push!(new_ne_indices, idx)
                    push!(result_ne, (sum_idx, idx))
                    # Generate diagonal term: substitute sum_idx -> idx
                    diag_term = change_index(term, sum_idx, idx)
                    push!(result_terms, diag_term)
                end
            end
        end

        push!(result_terms, term)
    end

    return QAdd(result_terms, result_indices, result_ne)
end

# Passthrough for non-sum types
expand_sums(m::QMul) = m
expand_sums(op::QSym) = op
```

### Task 11: Update `src/normal_order.jl`

**Files:**
- Modify: `src/normal_order.jl`

- [ ] **Step 1: Update Transition constructor and add ClusterSpace no-op**

Update the Transition constructor call to include `index`:
```julia
new_op = Transition(op.name, k, k, op.space_index, op.copy_index, op.index)
```

Drop type parameter from `_apply_ground_state` and `_expand_ground_state` — they work with `QMul`/`QAdd` (no `{CT}`).

### Task 12: Update `src/interface.jl`, `src/printing.jl`, `src/SecondQuantizedAlgebra.jl`

**Files:**
- Modify: `src/interface.jl`
- Modify: `src/printing.jl`
- Modify: `src/SecondQuantizedAlgebra.jl`

- [ ] **Step 1: Update `src/interface.jl`**

Remove type parameters from `maketerm`, `arguments`, etc. Key changes:
- `SymbolicUtils.arguments(a::QMul)` returns `[a.arg_c, a.args_nc...]` — `arg_c` is now `Num`
- `TermInterface.maketerm(::Type{QMul}, ...)` — wrap `arg_c` in `Num`
- `TermInterface.maketerm(::Type{QAdd}, ...)` — no type promotion needed

- [ ] **Step 2: Update `src/printing.jl`**

Add Index display to operator show methods. Add `Σ` display for `QAdd` with indices. Add `_show_copy_suffix` and `_show_index_suffix`:

```julia
function _show_index_suffix(io::IO, idx::Index)
    if has_index(idx)
        write(io, "_")
        print(io, idx.name)  # subscript style: σ_i
    end
    return
end
```

Update all operator `show` methods to call `_show_index_suffix(io, x.index)`.

Add QAdd show for sums:
```julia
function Base.show(io::IO, x::QAdd)
    if !isempty(x.indices)
        # Show Σ prefix
        for idx in x.indices
            write(io, "Σ($(idx.name)=1:$(idx.range))")
        end
        if !isempty(x.non_equal)
            write(io, "(")
            # Show constraints
            for (a, b) in x.non_equal
                write(io, "$(a.name)≠$(b.name),")
            end
            write(io, ")")
        end
        write(io, " ")
    end
    # Show terms (existing logic)
    ...
end
```

- [ ] **Step 3: Update `src/SecondQuantizedAlgebra.jl`**

Add includes and exports:
```julia
# After types.jl, before fock.jl
include("index.jl")

# After cluster.jl, before qmul.jl
include("indexed_variables.jl")

# After qadd.jl
include("qsum.jl")

# After commutator.jl
include("change_index.jl")
include("expand_sums.jl")
```

Add to exports:
```julia
    Index, NO_INDEX, has_index, IndexedOperator,
    IndexedVariable, DoubleIndexedVariable,
    Σ, ∑, change_index, expand_sums, get_indices,
```

Add Symbolics import:
```julia
using Symbolics: Symbolics, Num, @syms, unwrap, substitute
```

### Task 13: Update tests

**Files:**
- Modify: all test files to work with `Num` prefactors and `index` field
- Create: `test/index_test.jl`
- Create: `test/qsum_test.jl`

- [ ] **Step 1: Update existing tests for Num prefactors**

In all test files, `arg_c` comparisons like `result.arguments[1].arg_c == 1` need to work with `Num`. `Num(1) == 1` is `true`, so most tests should work unchanged. The `@inferred` tests may need updating since return types change from `QAdd{Complex{Int}}` to `QAdd`.

Key changes across test files:
- Remove all `QMul{Int}`, `QMul{Complex{Int}}`, `QAdd{Int}`, `QAdd{Complex{Int}}` type annotations
- Replace `QMul{CT}[...]` with `QMul[...]`
- `_ZERO_QADD` type changes from `QAdd{Complex{Int}}` to `QAdd`

- [ ] **Step 2: Create `test/index_test.jl`**

```julia
using SecondQuantizedAlgebra
using Symbolics: @variables
using Test

@testset "Index" begin
    @testset "Construction" begin
        hc = FockSpace(:c)
        ha = NLevelSpace(:atom, 2, 1)
        h = hc ⊗ ha
        @variables N
        i = Index(h, :i, N, ha)
        @test i.name == :i
        @test i.space_index == 2
        @test has_index(i)
        @test !has_index(NO_INDEX)
    end

    @testset "IndexedOperator" begin
        ha = NLevelSpace(:atom, 2, 1)
        hc = FockSpace(:c)
        h = hc ⊗ ha
        @variables N
        i = Index(h, :i, N, ha)

        σ12 = Transition(h, :σ, 1, 2, 2)
        σ_i = IndexedOperator(σ12, i)
        @test σ_i isa Transition
        @test σ_i.index == i
        @test σ_i.space_index == 2
        @test σ_i.i == 1 && σ_i.j == 2

        a = Destroy(h, :a, 1)
        a_i = IndexedOperator(a, i)
        @test a_i isa Destroy
        @test a_i.index == i
    end

    @testset "Equality with index" begin
        @variables N
        ha = NLevelSpace(:atom, 2, 1)
        h = FockSpace(:c) ⊗ ha
        i = Index(h, :i, N, ha)
        j = Index(h, :j, N, ha)

        σ_i = IndexedOperator(Transition(h, :σ, 1, 2, 2), i)
        σ_j = IndexedOperator(Transition(h, :σ, 1, 2, 2), j)
        σ_i2 = IndexedOperator(Transition(h, :σ, 1, 2, 2), i)

        @test σ_i == σ_i2
        @test σ_i != σ_j
        @test hash(σ_i) == hash(σ_i2)
        @test hash(σ_i) != hash(σ_j)
    end

    @testset "IndexedVariable" begin
        @variables N
        ha = NLevelSpace(:atom, 2, 1)
        h = FockSpace(:c) ⊗ ha
        i = Index(h, :i, N, ha)

        g_i = IndexedVariable(:g, i)
        @test g_i isa Num
    end

    @testset "DoubleIndexedVariable" begin
        @variables N
        ha = NLevelSpace(:atom, 2, 1)
        h = FockSpace(:c) ⊗ ha
        i = Index(h, :i, N, ha)
        j = Index(h, :j, N, ha)

        Γ_ij = DoubleIndexedVariable(:Γ, i, j)
        @test Γ_ij isa Num

        Γ_diag = DoubleIndexedVariable(:Γ, i, i; identical = false)
        @test isequal(Γ_diag, Num(0))
    end
end
```

- [ ] **Step 3: Create `test/qsum_test.jl`**

```julia
using SecondQuantizedAlgebra
using Symbolics: @variables
using Test
import SecondQuantizedAlgebra: simplify, _ZERO_QADD

@testset "QSum" begin
    hc = FockSpace(:c)
    ha = NLevelSpace(:atom, 2, 1)
    h = hc ⊗ ha
    @variables N Δ
    @qnumbers a::Destroy(h)

    i = Index(h, :i, N, ha)
    j = Index(h, :j, N, ha)

    σ(α, β, idx) = IndexedOperator(Transition(h, :σ, α, β, 2), idx)
    g(idx) = IndexedVariable(:g, idx)

    @testset "Σ construction" begin
        s = Σ(g(i) * a' * σ(1, 2, i), i)
        @test s isa QAdd
        @test length(s.indices) == 1
        @test s.indices[1] == i
        @test isempty(s.non_equal)
    end

    @testset "Σ with non_equal" begin
        s = Σ(σ(2, 1, i) * σ(1, 2, j), i, [j])
        @test length(s.non_equal) == 1
        @test s.non_equal[1] == (i, j)
    end

    @testset "Same index operators simplify" begin
        result = simplify(σ(1, 2, i) * σ(2, 1, i))
        @test result isa QAdd
    end

    @testset "Different index operators don't simplify" begin
        m = σ(1, 2, i) * σ(2, 1, j)
        result = simplify(m)
        @test length(result.arguments) == 1
        @test length(result.arguments[1].args_nc) == 2
    end

    @testset "Commutator with different index = 0" begin
        result = commutator(σ(1, 2, i), σ(2, 1, j))
        @test all(t -> iszero(t.arg_c), result.arguments)
    end

    @testset "change_index" begin
        σ_i = σ(1, 2, i)
        σ_j = change_index(σ_i, i, j)
        @test σ_j.index == j
        @test σ_j.i == 1 && σ_j.j == 2

        m = g(i) * a' * σ(1, 2, i)
        m_j = change_index(m, i, j)
        @test m_j.args_nc[end].index == j
    end

    @testset "get_indices" begin
        @test isempty(get_indices(a))
        @test get_indices(σ(1, 2, i)) == [i]

        m = σ(1, 2, i) * σ(2, 1, j)
        inds = get_indices(m)
        @test i in inds && j in inds
    end
end
```

- [ ] **Step 4: Run full test suite**

Run: `julia --project -e 'using Pkg; Pkg.test()'`

Fix any remaining issues from the Num migration (type annotations in test files, `@inferred` expectations, etc.).

### Task 14: Update `test/concrete_test.jl`

**Files:**
- Modify: `test/concrete_test.jl`

- [ ] **Step 1: Skip `ClusterSpace` and `Index` in concrete check**

`Index` has `Num` fields which are concrete but `Num.val` wraps `BasicSymbolic`. CheckConcreteStructs may flag this.

```julia
const _CONCRETE_SKIP = Set{Symbol}([:ClusterSpace, :Index])
```
