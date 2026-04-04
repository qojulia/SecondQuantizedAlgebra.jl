# ClusterSpace Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `copy_index::Int` to all `QSym` subtypes and introduce `ClusterSpace{T}` with `cluster_expand` for representing N identical quantum subsystems.

**Architecture:** Every `QSym` subtype gains a `copy_index::Int` field (default `1`). Canonical ordering, simplification, and commutation use `(space_index, copy_index)` pairs. New `ClusterSpace{T}` wraps a HilbertSpace with `N::T` and `order::Int`. An explicit `cluster_expand` function produces `order` copies of an operator.

**Tech Stack:** Julia, SecondQuantizedAlgebra type system, existing simplification/commutation infrastructure.

---

### Task 1: Add `copy_index` to Fock operators

**Files:**
- Modify: `src/fock.jl`

- [ ] **Step 1: Update structs and all methods in `src/fock.jl`**

Replace the entire file with:

```julia
"""
    Destroy <: QSym

Bosonic annihilation operator.
"""
struct Destroy <: QSym
    name::Symbol
    space_index::Int
    copy_index::Int
end
Destroy(name::Symbol, space_index::Int) = Destroy(name, space_index, 1)

"""
    Create <: QSym

Bosonic creation operator.
"""
struct Create <: QSym
    name::Symbol
    space_index::Int
    copy_index::Int
end
Create(name::Symbol, space_index::Int) = Create(name, space_index, 1)

# Construction from Hilbert spaces (validation, then discard)
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

# Adjoint
Base.adjoint(op::Destroy) = Create(op.name, op.space_index, op.copy_index)
Base.adjoint(op::Create) = Destroy(op.name, op.space_index, op.copy_index)

# Equality
Base.isequal(a::Destroy, b::Destroy) = a.name == b.name && a.space_index == b.space_index && a.copy_index == b.copy_index
Base.isequal(a::Create, b::Create) = a.name == b.name && a.space_index == b.space_index && a.copy_index == b.copy_index
Base.:(==)(a::Destroy, b::Destroy) = isequal(a, b)
Base.:(==)(a::Create, b::Create) = isequal(a, b)

# Hashing
Base.hash(a::Destroy, h::UInt) = hash(:Destroy, hash(a.name, hash(a.space_index, hash(a.copy_index, h))))
Base.hash(a::Create, h::UInt) = hash(:Create, hash(a.name, hash(a.space_index, hash(a.copy_index, h))))

# Canonical ordering: Create (0) before Destroy (1), then by space_index, then name
"""
    ladder(op::QSym)

Returns 0 for creation operators, 1 for annihilation operators.
Used for canonical ordering within `QMul.args_nc`.
"""
ladder(::Create) = 0
ladder(::Destroy) = 1
```

Note: `_unwrap_space(h::HilbertSpace) = h` must be added to `src/hilbertspace.jl` (see Task 4b below) so it exists before these files are loaded. The `ClusterSpace` specialization is added in Task 7.

### Task 2: Add `copy_index` to Transition

**Files:**
- Modify: `src/nlevel.jl`

- [ ] **Step 1: Update structs and all methods in `src/nlevel.jl`**

Replace the entire file with:

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

"""
    Transition <: QSym

Transition operator |i⟩⟨j| on an [`NLevelSpace`](@ref).
"""
struct Transition <: QSym
    name::Symbol
    i::Int
    j::Int
    space_index::Int
    copy_index::Int
end
Transition(name::Symbol, i::Int, j::Int, space_index::Int) = Transition(name, i, j, space_index, 1)

# Construction from Hilbert spaces
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

# Adjoint: |i⟩⟨j|† = |j⟩⟨i|
Base.adjoint(op::Transition) = Transition(op.name, op.j, op.i, op.space_index, op.copy_index)

# Equality
Base.isequal(a::Transition, b::Transition) = a.name == b.name && a.i == b.i && a.j == b.j && a.space_index == b.space_index && a.copy_index == b.copy_index
Base.:(==)(a::Transition, b::Transition) = isequal(a, b)

# Hashing
Base.hash(a::Transition, h::UInt) = hash(:Transition, hash(a.name, hash(a.i, hash(a.j, hash(a.space_index, hash(a.copy_index, h))))))

# Ladder (not applicable to Transition)
ladder(::Transition) = 0
```

### Task 3: Add `copy_index` to Pauli, Spin, Position, Momentum

**Files:**
- Modify: `src/pauli.jl`
- Modify: `src/spin.jl`
- Modify: `src/phase_space.jl`

- [ ] **Step 1: Update `src/pauli.jl`**

Replace the entire file with:

```julia
"""
    PauliSpace <: HilbertSpace

Hilbert space for two-level Pauli operators.
"""
struct PauliSpace <: HilbertSpace
    name::Symbol
end
Base.:(==)(a::PauliSpace, b::PauliSpace) = a.name == b.name
Base.hash(a::PauliSpace, h::UInt) = hash(:PauliSpace, hash(a.name, h))

"""
    Pauli <: QSym

Pauli operator (σx, σy, σz) on a [`PauliSpace`](@ref).
Axis: 1=x, 2=y, 3=z.
"""
struct Pauli <: QSym
    name::Symbol
    axis::Int
    space_index::Int
    copy_index::Int
    function Pauli(name::Symbol, axis::Int, space_index::Int, copy_index::Int)
        1 <= axis <= 3 || throw(ArgumentError("Pauli axis must be 1, 2, or 3, got $axis"))
        return new(name, axis, space_index, copy_index)
    end
end
Pauli(name::Symbol, axis::Int, space_index::Int) = Pauli(name, axis, space_index, 1)

# Construction from Hilbert spaces
Pauli(h::PauliSpace, name::Symbol, axis::Int) = Pauli(name, axis, 1)
function Pauli(h::ProductSpace, name::Symbol, axis::Int, idx::Int)
    1 <= idx <= length(h.spaces) || throw(ArgumentError("Index $idx out of range"))
    _unwrap_space(h.spaces[idx]) isa PauliSpace || throw(ArgumentError("Space at index $idx is not a PauliSpace"))
    return Pauli(name, axis, idx)
end

# Adjoint — Hermitian
Base.adjoint(op::Pauli) = op

# Equality
Base.isequal(a::Pauli, b::Pauli) = a.name == b.name && a.axis == b.axis && a.space_index == b.space_index && a.copy_index == b.copy_index
Base.:(==)(a::Pauli, b::Pauli) = isequal(a, b)

# Hashing
Base.hash(a::Pauli, h::UInt) = hash(:Pauli, hash(a.name, hash(a.axis, hash(a.space_index, hash(a.copy_index, h)))))

# Ladder (not applicable)
ladder(::Pauli) = 0
```

- [ ] **Step 2: Update `src/spin.jl`**

Replace the entire file with:

```julia
"""
    SpinSpace <: HilbertSpace

Hilbert space for collective spin operators (Spin-S).
"""
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

"""
    Spin <: QSym

Angular momentum operator (Sx, Sy, Sz) on a [`SpinSpace`](@ref).
Axis: 1=x, 2=y, 3=z.
"""
struct Spin <: QSym
    name::Symbol
    axis::Int
    space_index::Int
    copy_index::Int
    function Spin(name::Symbol, axis::Int, space_index::Int, copy_index::Int)
        1 <= axis <= 3 || throw(ArgumentError("Spin axis must be 1, 2, or 3, got $axis"))
        return new(name, axis, space_index, copy_index)
    end
end
Spin(name::Symbol, axis::Int, space_index::Int) = Spin(name, axis, space_index, 1)

# Construction from Hilbert spaces
Spin(h::SpinSpace, name::Symbol, axis::Int) = Spin(name, axis, 1)
function Spin(h::ProductSpace, name::Symbol, axis::Int, idx::Int)
    1 <= idx <= length(h.spaces) || throw(ArgumentError("Index $idx out of range"))
    _unwrap_space(h.spaces[idx]) isa SpinSpace || throw(ArgumentError("Space at index $idx is not a SpinSpace"))
    return Spin(name, axis, idx)
end

# Adjoint — Hermitian
Base.adjoint(op::Spin) = op

# Equality
Base.isequal(a::Spin, b::Spin) = a.name == b.name && a.axis == b.axis && a.space_index == b.space_index && a.copy_index == b.copy_index
Base.:(==)(a::Spin, b::Spin) = isequal(a, b)

# Hashing
Base.hash(a::Spin, h::UInt) = hash(:Spin, hash(a.name, hash(a.axis, hash(a.space_index, hash(a.copy_index, h)))))

# Ladder (not applicable)
ladder(::Spin) = 0
```

- [ ] **Step 3: Update `src/phase_space.jl`**

Replace the entire file with:

```julia
"""
    PhaseSpace <: HilbertSpace

Hilbert space for position and momentum (quadrature) operators.
"""
struct PhaseSpace <: HilbertSpace
    name::Symbol
end
Base.:(==)(a::PhaseSpace, b::PhaseSpace) = a.name == b.name
Base.hash(a::PhaseSpace, h::UInt) = hash(:PhaseSpace, hash(a.name, h))

"""
    Position <: QSym

Position (quadrature) operator on a [`PhaseSpace`](@ref).
"""
struct Position <: QSym
    name::Symbol
    space_index::Int
    copy_index::Int
end
Position(name::Symbol, space_index::Int) = Position(name, space_index, 1)

"""
    Momentum <: QSym

Momentum (quadrature) operator on a [`PhaseSpace`](@ref).
"""
struct Momentum <: QSym
    name::Symbol
    space_index::Int
    copy_index::Int
end
Momentum(name::Symbol, space_index::Int) = Momentum(name, space_index, 1)

# Construction from Hilbert spaces
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

# Adjoint — Hermitian
Base.adjoint(op::Position) = op
Base.adjoint(op::Momentum) = op

# Equality
Base.isequal(a::Position, b::Position) = a.name == b.name && a.space_index == b.space_index && a.copy_index == b.copy_index
Base.isequal(a::Momentum, b::Momentum) = a.name == b.name && a.space_index == b.space_index && a.copy_index == b.copy_index
Base.:(==)(a::Position, b::Position) = isequal(a, b)
Base.:(==)(a::Momentum, b::Momentum) = isequal(a, b)

# Hashing
Base.hash(a::Position, h::UInt) = hash(:Position, hash(a.name, hash(a.space_index, hash(a.copy_index, h))))
Base.hash(a::Momentum, h::UInt) = hash(:Momentum, hash(a.name, hash(a.space_index, hash(a.copy_index, h))))

# Ladder (not applicable — Hermitian operators, no creation/annihilation distinction)
ladder(::Position) = 0
ladder(::Momentum) = 1
```

### Task 4: Update canonical ordering in `qmul.jl` and add `_unwrap_space` base method

**Files:**
- Modify: `src/qmul.jl`
- Modify: `src/hilbertspace.jl`

- [ ] **Step 0: Add `_unwrap_space` base method to `src/hilbertspace.jl`**

At the end of `src/hilbertspace.jl` (after the `isless` methods), add:

```julia
# Unwrap ClusterSpace to get the original space for validation.
# Base method returns the space as-is; ClusterSpace overrides in cluster.jl.
_unwrap_space(h::HilbertSpace) = h
```

- [ ] **Step 1: Change all sort keys from `op -> op.space_index` to `op -> (op.space_index, op.copy_index)`**

There are 6 occurrences on lines 44, 59, 77, 85, 95, 116. Replace all:

```
sort!(args_nc; by = op -> op.space_index)
```
with:
```
sort!(args_nc; by = op -> (op.space_index, op.copy_index))
```

And the one on line 59:
```
sort!(args; by = op -> op.space_index)
```
with:
```
sort!(args; by = op -> (op.space_index, op.copy_index))
```

### Task 5: Update simplification rules to check `copy_index`

**Files:**
- Modify: `src/simplify.jl`

- [ ] **Step 1: Add `copy_index` checks to ordering-independent reductions**

On line 116, change:
```julia
if a isa Transition && b isa Transition && a.space_index == b.space_index && a.name == b.name
```
to:
```julia
if a isa Transition && b isa Transition && a.space_index == b.space_index && a.copy_index == b.copy_index && a.name == b.name
```

On line 118, change:
```julia
ops[i] = Transition(a.name, a.i, b.j, a.space_index)
```
to:
```julia
ops[i] = Transition(a.name, a.i, b.j, a.space_index, a.copy_index)
```

On line 128, change:
```julia
if a isa Pauli && b isa Pauli && a.space_index == b.space_index && a.name == b.name
```
to:
```julia
if a isa Pauli && b isa Pauli && a.space_index == b.space_index && a.copy_index == b.copy_index && a.name == b.name
```

On line 134, change:
```julia
ops[i] = Pauli(a.name, 6 - a.axis - b.axis, a.space_index)
```
to:
```julia
ops[i] = Pauli(a.name, 6 - a.axis - b.axis, a.space_index, a.copy_index)
```

- [ ] **Step 2: Add `copy_index` checks to ordering-dependent swaps**

On line 164, change:
```julia
if a isa Destroy && b isa Create && a.space_index == b.space_index && a.name == b.name
```
to:
```julia
if a isa Destroy && b isa Create && a.space_index == b.space_index && a.copy_index == b.copy_index && a.name == b.name
```

On line 174, change:
```julia
if a isa Spin && b isa Spin && a.space_index == b.space_index && a.name == b.name && a.axis > b.axis
```
to:
```julia
if a isa Spin && b isa Spin && a.space_index == b.space_index && a.copy_index == b.copy_index && a.name == b.name && a.axis > b.axis
```

On line 179, change:
```julia
ops[i] = Spin(a.name, 6 - a.axis - b.axis, a.space_index)
```
to:
```julia
ops[i] = Spin(a.name, 6 - a.axis - b.axis, a.space_index, a.copy_index)
```

On line 187, change:
```julia
if a isa Momentum && b isa Position && a.space_index == b.space_index
```
to:
```julia
if a isa Momentum && b isa Position && a.space_index == b.space_index && a.copy_index == b.copy_index
```

### Task 6: Update commutator short-circuits

**Files:**
- Modify: `src/commutator.jl`

- [ ] **Step 1: Update QSym short-circuit and helper functions**

On line 19, change:
```julia
a.space_index == b.space_index || return _ZERO_QADD
```
to:
```julia
(a.space_index == b.space_index && a.copy_index == b.copy_index) || return _ZERO_QADD
```

Replace `_has_space` (lines 41-46) with:
```julia
function _has_space(ops::Vector{QSym}, si::Int, ci::Int)
    for op in ops
        op.space_index == si && op.copy_index == ci && return true
    end
    return false
end
```

Replace `_has_shared_space` (lines 48-52) with:
```julia
function _has_shared_space(a::Vector{QSym}, b::Vector{QSym})
    for x in a, y in b
        x.space_index == y.space_index && x.copy_index == y.copy_index && return true
    end
    return false
end
```

Update call sites — line 26:
```julia
_has_space(a.args_nc, b.space_index) || return _ZERO_QADD
```
to:
```julia
_has_space(a.args_nc, b.space_index, b.copy_index) || return _ZERO_QADD
```

Line 30:
```julia
_has_space(b.args_nc, a.space_index) || return _ZERO_QADD
```
to:
```julia
_has_space(b.args_nc, a.space_index, a.copy_index) || return _ZERO_QADD
```

### Task 7: Create `src/cluster.jl`

**Files:**
- Create: `src/cluster.jl`

- [ ] **Step 1: Create `src/cluster.jl`**

```julia
"""
    ClusterSpace{T} <: HilbertSpace

Hilbert space representing N identical copies of another space, with
correlations tracked up to a specified order.

Fields:
- `original_space::HilbertSpace` — the space being replicated
- `N::T` — total number of copies (Int or symbolic, used by QC for scaling)
- `order::Int` — number of distinct copies to track in the algebra
"""
struct ClusterSpace{T} <: HilbertSpace
    original_space::HilbertSpace
    N::T
    order::Int
    function ClusterSpace(original_space::HilbertSpace, N::T, order::Int) where {T}
        order >= 1 || throw(ArgumentError("Order must be >= 1, got $order"))
        return new{T}(original_space, N, order)
    end
end

Base.:(==)(a::ClusterSpace, b::ClusterSpace) = a.original_space == b.original_space && a.N == b.N && a.order == b.order
Base.hash(a::ClusterSpace, h::UInt) = hash(:ClusterSpace, hash(a.original_space, hash(a.N, hash(a.order, h))))

# ClusterSpace specialization (base method in hilbertspace.jl)
_unwrap_space(h::ClusterSpace) = h.original_space

"""
    has_cluster(h::HilbertSpace) -> Bool

Check if a Hilbert space contains any `ClusterSpace`.
"""
has_cluster(::HilbertSpace) = false
has_cluster(::ClusterSpace) = true
function has_cluster(h::ProductSpace)
    for space in h.spaces
        space isa ClusterSpace && return true
    end
    return false
end

"""
    cluster_expand(op::QSym, order::Int) -> Vector{QSym}

Create `order` copies of `op` with `copy_index` set to `1, 2, ..., order`
and name suffixed `_1`, `_2`, etc.
"""
function cluster_expand(op::QSym, order::Int)
    return [_with_copy(op, Symbol(op.name, :_, i), i) for i in 1:order]
end

"""
    cluster_expand(op::QSym, h::ProductSpace) -> Vector{QSym}

Create copies of `op` using the order from the `ClusterSpace` at `op.space_index`.
"""
function cluster_expand(op::QSym, h::ProductSpace)
    space = h.spaces[op.space_index]
    space isa ClusterSpace || throw(ArgumentError("Space at index $(op.space_index) is not a ClusterSpace"))
    return cluster_expand(op, space.order)
end

# Per-type copy constructors
_with_copy(op::Destroy, name::Symbol, ci::Int) = Destroy(name, op.space_index, ci)
_with_copy(op::Create, name::Symbol, ci::Int) = Create(name, op.space_index, ci)
_with_copy(op::Transition, name::Symbol, ci::Int) = Transition(name, op.i, op.j, op.space_index, ci)
_with_copy(op::Pauli, name::Symbol, ci::Int) = Pauli(name, op.axis, op.space_index, ci)
_with_copy(op::Spin, name::Symbol, ci::Int) = Spin(name, op.axis, op.space_index, ci)
_with_copy(op::Position, name::Symbol, ci::Int) = Position(name, op.space_index, ci)
_with_copy(op::Momentum, name::Symbol, ci::Int) = Momentum(name, op.space_index, ci)
```

### Task 8: Update normal_order.jl

**Files:**
- Modify: `src/normal_order.jl`

- [ ] **Step 1: Update Transition constructor call and add ClusterSpace no-op**

On line 52, change:
```julia
new_op = Transition(op.name, k, k, op.space_index)
```
to:
```julia
new_op = Transition(op.name, k, k, op.space_index, op.copy_index)
```

After line 37 (`_apply_ground_state(expr::QAdd, ::PhaseSpace) = expr`), add:
```julia
_apply_ground_state(expr::QAdd, ::ClusterSpace) = expr
```

### Task 9: Update printing

**Files:**
- Modify: `src/printing.jl`

- [ ] **Step 1: Add copy_index suffix display and ClusterSpace show**

After line 20 (the PhaseSpace show), add:
```julia
function Base.show(io::IO, h::ClusterSpace)
    write(io, "Cluster(")
    show(io, h.original_space)
    write(io, ", N=")
    print(io, h.N)
    write(io, ", order=")
    print(io, h.order)
    return write(io, ")")
end
```

Add a helper function for copy_index suffix. After the `_write_subscript` function (line 36), add:
```julia
function _show_copy_suffix(io::IO, ci::Int)
    if ci > 1
        write(io, "_")
        print(io, ci)
    end
    return
end
```

Update operator show methods to append copy suffix:

Line 23, change:
```julia
Base.show(io::IO, x::Destroy) = print(io, x.name)
```
to:
```julia
Base.show(io::IO, x::Destroy) = (print(io, x.name); _show_copy_suffix(io, x.copy_index))
```

Line 24, change:
```julia
Base.show(io::IO, x::Create) = (print(io, x.name); write(io, "†"))
```
to:
```julia
Base.show(io::IO, x::Create) = (print(io, x.name); _show_copy_suffix(io, x.copy_index); write(io, "†"))
```

Line 37-41 (Transition show), change to:
```julia
function Base.show(io::IO, x::Transition)
    print(io, x.name)
    _show_copy_suffix(io, x.copy_index)
    _write_subscript(io, x.i)
    return _write_subscript(io, x.j)
end
```

Line 45, change:
```julia
Base.show(io::IO, x::Pauli) = (print(io, x.name); print(io, _xyz_sym[x.axis]))
```
to:
```julia
Base.show(io::IO, x::Pauli) = (print(io, x.name); _show_copy_suffix(io, x.copy_index); print(io, _xyz_sym[x.axis]))
```

Line 46, change:
```julia
Base.show(io::IO, x::Spin) = (print(io, x.name); print(io, _xyz_sym[x.axis]))
```
to:
```julia
Base.show(io::IO, x::Spin) = (print(io, x.name); _show_copy_suffix(io, x.copy_index); print(io, _xyz_sym[x.axis]))
```

Line 49-50, change:
```julia
Base.show(io::IO, x::Position) = print(io, x.name)
Base.show(io::IO, x::Momentum) = print(io, x.name)
```
to:
```julia
Base.show(io::IO, x::Position) = (print(io, x.name); _show_copy_suffix(io, x.copy_index))
Base.show(io::IO, x::Momentum) = (print(io, x.name); _show_copy_suffix(io, x.copy_index))
```

### Task 10: Wire into module and exports

**Files:**
- Modify: `src/SecondQuantizedAlgebra.jl`

- [ ] **Step 1: Add include and exports**

After line 18 (`include("phase_space.jl")`), add:
```julia
include("cluster.jl")
```

In the export block, add after `PhaseSpace, Position, Momentum,`:
```julia
    ClusterSpace, cluster_expand, has_cluster,
```

### Task 11: Update concrete_test.jl to skip ClusterSpace

**Files:**
- Modify: `test/concrete_test.jl`

- [ ] **Step 1: Add ClusterSpace to skip set**

Change line 5:
```julia
const _CONCRETE_SKIP = Set{Symbol}()
```
to:
```julia
const _CONCRETE_SKIP = Set{Symbol}([:ClusterSpace])
```

### Task 12: Write tests

**Files:**
- Create: `test/cluster_test.jl`

- [ ] **Step 1: Create `test/cluster_test.jl`**

```julia
using SecondQuantizedAlgebra
using Test
import SecondQuantizedAlgebra: _unwrap_space, has_cluster, _ZERO_QADD, simplify

@testset "ClusterSpace" begin
    @testset "Construction" begin
        ha = NLevelSpace(:atom, 2, 1)
        c = ClusterSpace(ha, 10, 2)
        @test c.original_space == ha
        @test c.N == 10
        @test c.order == 2
        @test _unwrap_space(c) === ha
        @test _unwrap_space(ha) === ha
    end

    @testset "Equality and hashing" begin
        ha = NLevelSpace(:atom, 2, 1)
        c1 = ClusterSpace(ha, 10, 2)
        c2 = ClusterSpace(ha, 10, 2)
        c3 = ClusterSpace(ha, 20, 2)
        @test c1 == c2
        @test c1 != c3
        @test hash(c1) == hash(c2)
    end

    @testset "has_cluster" begin
        ha = NLevelSpace(:atom, 2, 1)
        hf = FockSpace(:c)
        c = ClusterSpace(ha, 10, 2)
        @test has_cluster(c) == true
        @test has_cluster(ha) == false
        @test has_cluster(hf) == false

        h = hf ⊗ c
        @test has_cluster(h) == true
        h2 = hf ⊗ ha
        @test has_cluster(h2) == false
    end

    @testset "ProductSpace with ClusterSpace" begin
        ha = NLevelSpace(:atom, 2, 1)
        hf = FockSpace(:c)
        c = ClusterSpace(ha, 10, 2)
        h = hf ⊗ c

        # Operators on FockSpace position work normally
        a = Destroy(h, :a, 1)
        @test a isa Destroy
        @test a.space_index == 1
        @test a.copy_index == 1

        # Operators on ClusterSpace position validate against original_space
        σ = Transition(h, :σ, 1, 2, 2)
        @test σ isa Transition
        @test σ.space_index == 2
        @test σ.copy_index == 1
    end

    @testset "cluster_expand" begin
        ha = NLevelSpace(:atom, 2, 1)
        hf = FockSpace(:c)
        c = ClusterSpace(ha, 10, 2)
        h = hf ⊗ c

        σ = Transition(h, :σ, 1, 2, 2)
        copies = cluster_expand(σ, 2)
        @test length(copies) == 2
        @test copies[1].name == :σ_1
        @test copies[1].copy_index == 1
        @test copies[1].space_index == 2
        @test copies[2].name == :σ_2
        @test copies[2].copy_index == 2
        @test copies[2].space_index == 2

        # Expand via ProductSpace
        copies2 = cluster_expand(σ, h)
        @test length(copies2) == 2
        @test copies2[1].name == :σ_1
        @test copies2[2].name == :σ_2

        # Fock operator expand
        a = Destroy(h, :a, 1)
        @test_throws ArgumentError cluster_expand(a, h)  # FockSpace, not ClusterSpace
        fock_copies = cluster_expand(a, 3)
        @test length(fock_copies) == 3
        @test fock_copies[1].name == :a_1
        @test fock_copies[2].copy_index == 2
    end

    @testset "copy_index field on all QSym types" begin
        # Default copy_index = 1
        @test Destroy(:a, 1).copy_index == 1
        @test Create(:a, 1).copy_index == 1
        @test Transition(:σ, 1, 2, 1).copy_index == 1
        @test Pauli(:σ, 1, 1).copy_index == 1
        @test Spin(:S, 1, 1).copy_index == 1
        @test Position(:x, 1).copy_index == 1
        @test Momentum(:p, 1).copy_index == 1

        # Explicit copy_index
        @test Destroy(:a, 1, 2).copy_index == 2
        @test Create(:a, 1, 3).copy_index == 3
        @test Transition(:σ, 1, 2, 1, 2).copy_index == 2
    end

    @testset "Adjoint preserves copy_index" begin
        a = Destroy(:a, 1, 2)
        @test a'.copy_index == 2
        @test a'.name == :a

        σ = Transition(:σ, 1, 2, 1, 3)
        @test σ'.copy_index == 3
        @test σ'.i == 2 && σ'.j == 1
    end

    @testset "Equality includes copy_index" begin
        a1 = Destroy(:a, 1, 1)
        a2 = Destroy(:a, 1, 2)
        @test a1 != a2
        @test !isequal(a1, a2)
        @test hash(a1) != hash(a2)
    end

    @testset "Canonical ordering uses (space_index, copy_index)" begin
        a1 = Destroy(:a, 1, 1)
        a2 = Destroy(:a, 1, 2)
        m = a2 * a1  # should sort by copy_index
        @test m.args_nc[1].copy_index == 1
        @test m.args_nc[2].copy_index == 2
    end

    @testset "Different copy_index operators commute" begin
        a1 = Destroy(:a, 1, 1)
        a2 = Destroy(:a, 1, 2)

        # Different copy → short-circuit to zero
        @test commutator(a1, a2') === _ZERO_QADD

        # Same copy → non-trivial
        result = commutator(a1, a1')
        @test result isa QAdd
        @test !all(iszero, result.arguments)
    end

    @testset "Simplify respects copy_index" begin
        a1 = Destroy(:a, 1, 1)
        a2 = Destroy(:a, 1, 2)

        # Same copy: a·a† → a†·a + 1
        result = simplify(a1 * a1')
        @test length(result.arguments) == 2

        # Different copy: a₁·a₂† stays as-is (no commutation rule)
        result = simplify(a1 * a2')
        @test length(result.arguments) == 1
        @test length(result.arguments[1].args_nc) == 2
    end

    @testset "Symbolic N" begin
        using Symbolics: @variables
        @variables N_sym
        ha = NLevelSpace(:atom, 2, 1)
        c = ClusterSpace(ha, N_sym, 2)
        @test c.N === N_sym
        @test c.order == 2
    end

    @testset "Type stability" begin
        ha = NLevelSpace(:atom, 2, 1)
        hf = FockSpace(:c)
        c = ClusterSpace(ha, 10, 2)
        h = hf ⊗ c

        a = Destroy(h, :a, 1)
        σ = Transition(h, :σ, 1, 2, 2)

        @inferred commutator(a, σ)
        @inferred simplify(a * a')

        a1 = Destroy(:a, 1, 1)
        a2 = Destroy(:a, 1, 2)
        @inferred commutator(a1, a2')
        @inferred commutator(a1, a1')
        @inferred simplify(a1 * a2')
    end
end
```

### Task 13: Run tests

- [ ] **Step 1: Run cluster tests**

Run: `julia --project -e 'using Pkg; Pkg.test(test_args=["cluster_test"])'`

Expected: All tests PASS.

- [ ] **Step 2: Run full test suite**

Run: `julia --project -e 'using Pkg; Pkg.test()'`

Expected: All tests PASS.
