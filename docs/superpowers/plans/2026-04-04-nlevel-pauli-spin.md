# NLevel, Pauli, and Spin Operators Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add NLevelSpace/Transition, PauliSpace/Pauli, and SpinSpace/Spin operators to the type-stable core, with lazy algebra rules applied in `normal_order()`.

**Architecture:** Each operator family gets its own source file following the Fock pattern (struct + constructors + adjoint + equality + hashing). The `normal_order` function is extended to handle all new operator pair rules. Existing files (hilbertspace.jl, printing.jl, latexify_recipes.jl, numeric.jl) are extended with new methods.

**Tech Stack:** Julia, SymbolicUtils v4, Symbolics v7, QuantumOpticsBase, Latexify, Combinatorics (for `levicivita`)

---

## File Structure

**New files:**
- `src/nlevel.jl` — NLevelSpace, Transition struct, constructors, adjoint, equality, hashing
- `src/pauli.jl` — PauliSpace, Pauli struct, constructors, adjoint, equality, hashing
- `src/spin.jl` — SpinSpace, Spin struct, constructors, adjoint, equality, hashing
- `test/nlevel_test.jl` — Tests for NLevel
- `test/pauli_test.jl` — Tests for Pauli
- `test/spin_test.jl` — Tests for Spin

**Modified files:**
- `src/SecondQuantizedAlgebra.jl` — Add includes, imports (`levicivita`), exports
- `src/hilbertspace.jl` — Add generic `⊗` fallback for non-FockSpace types
- `src/fock.jl` — Add `ladder` defaults for new operator types
- `src/normal_order.jl` — Extend `_normal_order_qmul` for Transition, Pauli, Spin rules
- `src/printing.jl` — Add `show` methods for new types
- `src/latexify_recipes.jl` — Add `@latexrecipe` for new types, `transition_superscript` toggle
- `src/numeric.jl` — Add `to_numeric` methods for new types and bases
- `test/normal_order_test.jl` — Add tests for new algebra rules
- `test/integration_test.jl` — Add multi-type integration tests

---

### Task 1: Generalize `⊗` for Non-Fock Spaces

**Files:**
- Modify: `src/hilbertspace.jl`

Currently `⊗` has specialized methods for `FockSpace` only. We need a generic binary fallback so `NLevelSpace(:a, 2, 1) ⊗ FockSpace(:c)` works.

- [ ] **Step 1: Add generic `⊗` methods**

Add to `src/hilbertspace.jl`, replacing the FockSpace-specific binary methods:

```julia
⊗(a::HilbertSpace, b::HilbertSpace) = ProductSpace((a, b))
⊗(a::ProductSpace, b::HilbertSpace) = ProductSpace((a.spaces..., b))
⊗(a::HilbertSpace, b::ProductSpace) = ProductSpace((a, b.spaces...))
⊗(a::ProductSpace, b::ProductSpace) = ProductSpace((a.spaces..., b.spaces...))
⊗(a::HilbertSpace, b::HilbertSpace, c::HilbertSpace...) = ⊗(a ⊗ b, c...)
⊗(a::HilbertSpace) = a
```

This replaces the four `FockSpace`-specific methods with generic `HilbertSpace` methods.

- [ ] **Step 2: Verify existing tests still pass**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: All existing tests pass (the generic methods subsume the FockSpace-specific ones).

---

### Task 2: NLevelSpace and Transition

**Files:**
- Create: `src/nlevel.jl`
- Create: `test/nlevel_test.jl`
- Modify: `src/SecondQuantizedAlgebra.jl`
- Modify: `src/fock.jl` (add `ladder` default)

- [ ] **Step 1: Write the test file**

```julia
# test/nlevel_test.jl
using SecondQuantizedAlgebra
using Test

@testset "nlevel" begin
    @testset "NLevelSpace" begin
        h = NLevelSpace(:atom, 3, 1)
        @test h isa HilbertSpace
        @test h.name == :atom
        @test h.n == 3
        @test h.ground_state == 1
        @test NLevelSpace(:a, 3, 1) == NLevelSpace(:a, 3, 1)
        @test NLevelSpace(:a, 3, 1) != NLevelSpace(:b, 3, 1)
    end

    @testset "Transition construction — single space" begin
        h = NLevelSpace(:atom, 3, 1)
        σ = Transition(h, :σ, 1, 2)
        @test σ isa Transition
        @test σ isa QSym
        @test σ.name == :σ
        @test σ.i == 1
        @test σ.j == 2
        @test σ.space_index == 1
    end

    @testset "Transition construction — product space" begin
        h = FockSpace(:c) ⊗ NLevelSpace(:atom, 2, 1)
        σ = Transition(h, :σ, 1, 2, 2)
        @test σ.space_index == 2
        @test_throws AssertionError Transition(h, :σ, 1, 2, 1)  # FockSpace, not NLevel
    end

    @testset "Adjoint" begin
        h = NLevelSpace(:atom, 3, 1)
        σ12 = Transition(h, :σ, 1, 2)
        σ21 = σ12'
        @test σ21 isa Transition
        @test σ21.i == 2
        @test σ21.j == 1
        @test σ12'' == σ12
    end

    @testset "Equality and hashing" begin
        h = NLevelSpace(:atom, 3, 1)
        σ1 = Transition(h, :σ, 1, 2)
        σ2 = Transition(h, :σ, 1, 2)
        σ3 = Transition(h, :σ, 2, 1)
        @test isequal(σ1, σ2)
        @test !isequal(σ1, σ3)
        @test hash(σ1) == hash(σ2)
        @test hash(σ1) != hash(σ3)
    end

    @testset "Arithmetic" begin
        h = NLevelSpace(:atom, 3, 1)
        σ12 = Transition(h, :σ, 1, 2)
        σ21 = Transition(h, :σ, 2, 1)

        m = σ12 * σ21
        @test m isa QMul{Int}
        @test length(m.args_nc) == 2

        s = σ12 + σ21
        @test s isa QAdd{Int}
    end

    @testset "@qnumbers" begin
        h = NLevelSpace(:atom, 2, 1)
        @qnumbers σ::Transition(h, 1, 2)
        @test σ isa Transition
        @test σ.name == :σ
    end

    @testset "Concrete types" begin
        using CheckConcreteStructs
        all_concrete(NLevelSpace)
        all_concrete(Transition)
    end
end
```

- [ ] **Step 2: Write the implementation**

```julia
# src/nlevel.jl

"""
    NLevelSpace <: HilbertSpace

Hilbert space for N-level systems (atoms, qubits, etc.).
"""
struct NLevelSpace <: HilbertSpace
    name::Symbol
    n::Int
    ground_state::Int
    function NLevelSpace(name::Symbol, n::Int, ground_state::Int)
        @assert 1 <= ground_state <= n "Ground state $ground_state out of range 1:$n"
        new(name, n, ground_state)
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
end

# Construction from Hilbert spaces
function Transition(h::NLevelSpace, name::Symbol, i::Int, j::Int)
    @assert 1 <= i <= h.n "Level i=$i out of range 1:$(h.n)"
    @assert 1 <= j <= h.n "Level j=$j out of range 1:$(h.n)"
    return Transition(name, i, j, 1)
end
function Transition(h::ProductSpace, name::Symbol, i::Int, j::Int, idx::Int)
    @assert 1 <= idx <= length(h.spaces) "Index $idx out of range"
    space = h.spaces[idx]
    @assert space isa NLevelSpace "Space at index $idx is not an NLevelSpace"
    @assert 1 <= i <= space.n "Level i=$i out of range 1:$(space.n)"
    @assert 1 <= j <= space.n "Level j=$j out of range 1:$(space.n)"
    return Transition(name, i, j, idx)
end

# Adjoint: |i⟩⟨j|† = |j⟩⟨i|
Base.adjoint(op::Transition) = Transition(op.name, op.j, op.i, op.space_index)

# Equality
Base.isequal(a::Transition, b::Transition) = a.name == b.name && a.i == b.i && a.j == b.j && a.space_index == b.space_index
Base.:(==)(a::Transition, b::Transition) = isequal(a, b)

# Hashing
Base.hash(a::Transition, h::UInt) = hash(:Transition, hash(a.name, hash(a.i, hash(a.j, hash(a.space_index, h)))))

# Ladder (not applicable to Transition, default value)
ladder(::Transition) = 0
```

Update `src/SecondQuantizedAlgebra.jl` — add `include("nlevel.jl")` after `include("fock.jl")`, and add to exports: `NLevelSpace, Transition`.

- [ ] **Step 3: Run tests**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: All tests pass including new nlevel_test.

---

### Task 3: PauliSpace and Pauli

**Files:**
- Create: `src/pauli.jl`
- Create: `test/pauli_test.jl`
- Modify: `src/SecondQuantizedAlgebra.jl`

- [ ] **Step 1: Write the test file**

```julia
# test/pauli_test.jl
using SecondQuantizedAlgebra
using Test

@testset "pauli" begin
    @testset "PauliSpace" begin
        h = PauliSpace(:p)
        @test h isa HilbertSpace
        @test h.name == :p
        @test PauliSpace(:p) == PauliSpace(:p)
        @test PauliSpace(:p) != PauliSpace(:q)
    end

    @testset "Pauli construction — single space" begin
        h = PauliSpace(:p)
        σx = Pauli(h, :σ, 1)
        σy = Pauli(h, :σ, 2)
        σz = Pauli(h, :σ, 3)
        @test σx isa Pauli
        @test σx isa QSym
        @test σx.axis == 1
        @test σy.axis == 2
        @test σz.axis == 3
        @test σx.space_index == 1
        @test_throws AssertionError Pauli(h, :σ, 4)
    end

    @testset "Pauli construction — product space" begin
        h = FockSpace(:c) ⊗ PauliSpace(:p)
        σx = Pauli(h, :σ, 1, 2)
        @test σx.space_index == 2
        @test_throws AssertionError Pauli(h, :σ, 1, 1)  # FockSpace, not Pauli
    end

    @testset "Adjoint — Hermitian" begin
        h = PauliSpace(:p)
        σx = Pauli(h, :σ, 1)
        @test σx' == σx
    end

    @testset "Equality and hashing" begin
        h = PauliSpace(:p)
        σ1 = Pauli(h, :σ, 1)
        σ2 = Pauli(h, :σ, 1)
        σ3 = Pauli(h, :σ, 2)
        @test isequal(σ1, σ2)
        @test !isequal(σ1, σ3)
        @test hash(σ1) == hash(σ2)
        @test hash(σ1) != hash(σ3)
    end

    @testset "Arithmetic" begin
        h = PauliSpace(:p)
        σx = Pauli(h, :σ, 1)
        σy = Pauli(h, :σ, 2)

        m = σx * σy
        @test m isa QMul{Int}
        @test length(m.args_nc) == 2

        s = σx + σy
        @test s isa QAdd{Int}
    end

    @testset "@qnumbers" begin
        h = PauliSpace(:p)
        @qnumbers σx::Pauli(h, 1)
        @test σx isa Pauli
        @test σx.name == :σx
        @test σx.axis == 1
    end

    @testset "Concrete types" begin
        using CheckConcreteStructs
        all_concrete(PauliSpace)
        all_concrete(Pauli)
    end
end
```

- [ ] **Step 2: Write the implementation**

```julia
# src/pauli.jl

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
    function Pauli(name::Symbol, axis::Int, space_index::Int)
        @assert 1 <= axis <= 3 "Pauli axis must be 1, 2, or 3, got $axis"
        new(name, axis, space_index)
    end
end

# Construction from Hilbert spaces
function Pauli(h::PauliSpace, name::Symbol, axis::Int)
    return Pauli(name, axis, 1)
end
function Pauli(h::ProductSpace, name::Symbol, axis::Int, idx::Int)
    @assert 1 <= idx <= length(h.spaces) "Index $idx out of range"
    @assert h.spaces[idx] isa PauliSpace "Space at index $idx is not a PauliSpace"
    return Pauli(name, axis, idx)
end

# Adjoint — Hermitian
Base.adjoint(op::Pauli) = op

# Equality
Base.isequal(a::Pauli, b::Pauli) = a.name == b.name && a.axis == b.axis && a.space_index == b.space_index
Base.:(==)(a::Pauli, b::Pauli) = isequal(a, b)

# Hashing
Base.hash(a::Pauli, h::UInt) = hash(:Pauli, hash(a.name, hash(a.axis, hash(a.space_index, h))))

# Ladder (not applicable)
ladder(::Pauli) = 0
```

Update `src/SecondQuantizedAlgebra.jl` — add `include("pauli.jl")` after `include("nlevel.jl")`, and add to exports: `PauliSpace, Pauli`.

- [ ] **Step 3: Run tests**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: All tests pass including new pauli_test.

---

### Task 4: SpinSpace and Spin

**Files:**
- Create: `src/spin.jl`
- Create: `test/spin_test.jl`
- Modify: `src/SecondQuantizedAlgebra.jl`

- [ ] **Step 1: Write the test file**

```julia
# test/spin_test.jl
using SecondQuantizedAlgebra
using Test

@testset "spin" begin
    @testset "SpinSpace" begin
        h = SpinSpace(:s, 1 // 2)
        @test h isa HilbertSpace
        @test h.name == :s
        @test h.spin == 1 // 2
        @test SpinSpace(:s, 1 // 2) == SpinSpace(:s, 1 // 2)
        @test SpinSpace(:s, 1 // 2) != SpinSpace(:s, 1)
    end

    @testset "Spin construction — single space" begin
        h = SpinSpace(:s, 1 // 2)
        Sx = Spin(h, :S, 1)
        Sy = Spin(h, :S, 2)
        Sz = Spin(h, :S, 3)
        @test Sx isa Spin
        @test Sx isa QSym
        @test Sx.axis == 1
        @test Sy.axis == 2
        @test Sz.axis == 3
        @test Sx.space_index == 1
        @test_throws AssertionError Spin(h, :S, 4)
    end

    @testset "Spin construction — product space" begin
        h = FockSpace(:c) ⊗ SpinSpace(:s, 1)
        Sx = Spin(h, :S, 1, 2)
        @test Sx.space_index == 2
        @test_throws AssertionError Spin(h, :S, 1, 1)  # FockSpace, not Spin
    end

    @testset "Adjoint — Hermitian" begin
        h = SpinSpace(:s, 1 // 2)
        Sx = Spin(h, :S, 1)
        @test Sx' == Sx
    end

    @testset "Equality and hashing" begin
        h = SpinSpace(:s, 1 // 2)
        S1 = Spin(h, :S, 1)
        S2 = Spin(h, :S, 1)
        S3 = Spin(h, :S, 2)
        @test isequal(S1, S2)
        @test !isequal(S1, S3)
        @test hash(S1) == hash(S2)
        @test hash(S1) != hash(S3)
    end

    @testset "Arithmetic" begin
        h = SpinSpace(:s, 1 // 2)
        Sx = Spin(h, :S, 1)
        Sy = Spin(h, :S, 2)

        m = Sx * Sy
        @test m isa QMul{Int}
        @test length(m.args_nc) == 2

        s = Sx + Sy
        @test s isa QAdd{Int}
    end

    @testset "@qnumbers" begin
        h = SpinSpace(:s, 1 // 2)
        @qnumbers Sx::Spin(h, 1)
        @test Sx isa Spin
        @test Sx.name == :Sx
        @test Sx.axis == 1
    end

    @testset "Concrete types" begin
        using CheckConcreteStructs
        all_concrete(SpinSpace)
        all_concrete(Spin)
    end
end
```

- [ ] **Step 2: Write the implementation**

```julia
# src/spin.jl

"""
    SpinSpace <: HilbertSpace

Hilbert space for collective spin operators (Spin-S).
"""
struct SpinSpace <: HilbertSpace
    name::Symbol
    spin::Rational{Int}
    function SpinSpace(name::Symbol, spin::Rational{Int})
        @assert spin > 0 "Spin must be positive, got $spin"
        @assert denominator(spin) in (1, 2) "Spin must be integer or half-integer, got $spin"
        new(name, spin)
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
    function Spin(name::Symbol, axis::Int, space_index::Int)
        @assert 1 <= axis <= 3 "Spin axis must be 1, 2, or 3, got $axis"
        new(name, axis, space_index)
    end
end

# Construction from Hilbert spaces
function Spin(h::SpinSpace, name::Symbol, axis::Int)
    return Spin(name, axis, 1)
end
function Spin(h::ProductSpace, name::Symbol, axis::Int, idx::Int)
    @assert 1 <= idx <= length(h.spaces) "Index $idx out of range"
    @assert h.spaces[idx] isa SpinSpace "Space at index $idx is not a SpinSpace"
    return Spin(name, axis, idx)
end

# Adjoint — Hermitian
Base.adjoint(op::Spin) = op

# Equality
Base.isequal(a::Spin, b::Spin) = a.name == b.name && a.axis == b.axis && a.space_index == b.space_index
Base.:(==)(a::Spin, b::Spin) = isequal(a, b)

# Hashing
Base.hash(a::Spin, h::UInt) = hash(:Spin, hash(a.name, hash(a.axis, hash(a.space_index, h))))

# Ladder (not applicable)
ladder(::Spin) = 0
```

Update `src/SecondQuantizedAlgebra.jl` — add `include("spin.jl")` after `include("pauli.jl")`, add `using Combinatorics: levicivita` to imports, and add to exports: `SpinSpace, Spin`.

- [ ] **Step 3: Run tests**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: All tests pass including new spin_test.

---

### Task 5: Extend `normal_order` for Transition, Pauli, Spin

**Files:**
- Modify: `src/normal_order.jl`
- Modify: `test/normal_order_test.jl`

- [ ] **Step 1: Add tests for new algebra rules**

Append to `test/normal_order_test.jl`:

```julia
@testset "normal_order — Transition" begin
    h = NLevelSpace(:atom, 3, 1)

    @testset "Composition: |1⟩⟨2| · |2⟩⟨3| = |1⟩⟨3|" begin
        σ12 = Transition(h, :σ, 1, 2)
        σ23 = Transition(h, :σ, 2, 3)
        m = σ12 * σ23
        result = normal_order(m)
        @test result isa QAdd
        @test length(result.arguments) == 1
        op = result.arguments[1].args_nc[1]
        @test op isa Transition
        @test op.i == 1 && op.j == 3
    end

    @testset "Orthogonality: |1⟩⟨2| · |3⟩⟨1| = 0" begin
        σ12 = Transition(h, :σ, 1, 2)
        σ31 = Transition(h, :σ, 3, 1)
        m = σ12 * σ31
        result = normal_order(m)
        @test result isa QAdd
        @test all(iszero, result.arguments)
    end

    @testset "Ground state rewriting" begin
        σ11 = Transition(h, :σ, 1, 1)
        m = QMul(1, QSym[σ11])
        result = normal_order(m, h)
        # |1⟩⟨1| = 1 - |2⟩⟨2| - |3⟩⟨3| for ground_state=1, n=3
        @test result isa QAdd
        @test length(result.arguments) == 3  # identity + two diagonal transitions
    end
end

@testset "normal_order — Pauli" begin
    h = PauliSpace(:p)
    σx = Pauli(h, :σ, 1)
    σy = Pauli(h, :σ, 2)
    σz = Pauli(h, :σ, 3)

    @testset "σx·σx = 1" begin
        m = σx * σx
        result = normal_order(m)
        @test result isa QAdd
        @test length(result.arguments) == 1
        @test isempty(result.arguments[1].args_nc)  # scalar
        @test result.arguments[1].arg_c == 1
    end

    @testset "σx·σy = iσz" begin
        m = σx * σy
        result = normal_order(m)
        @test result isa QAdd
        @test length(result.arguments) == 1
        t = result.arguments[1]
        @test t.arg_c == im
        @test length(t.args_nc) == 1
        @test t.args_nc[1] isa Pauli
        @test t.args_nc[1].axis == 3
    end

    @testset "σy·σx = -iσz" begin
        m = σy * σx
        result = normal_order(m)
        @test result isa QAdd
        @test length(result.arguments) == 1
        t = result.arguments[1]
        @test t.arg_c == -im
        @test t.args_nc[1].axis == 3
    end
end

@testset "normal_order — Spin" begin
    h = SpinSpace(:s, 1 // 2)
    Sx = Spin(h, :S, 1)
    Sy = Spin(h, :S, 2)
    Sz = Spin(h, :S, 3)

    @testset "[Sx, Sy] = iSz (swap + commutator)" begin
        # Sx * Sy stays as-is (axes in order), Sy * Sx gets reordered
        m = Sy * Sx  # out of axis order
        result = normal_order(m)
        @test result isa QAdd
        # Should produce: Sx·Sy - i·Sz (swap + commutator term)
        @test length(result.arguments) >= 2
    end

    @testset "Same axis — already ordered" begin
        m = Sx * Sx
        result = normal_order(m)
        @test result isa QAdd
        @test length(result.arguments) == 1  # no reordering needed
    end
end
```

- [ ] **Step 2: Extend `_normal_order_qmul` in `src/normal_order.jl`**

Replace the function with an extended version that handles all operator types. The scanning loop checks for adjacent pairs in priority order:

```julia
function _normal_order_qmul(arg_c::T, ops::Vector{QSym}) where {T}
    isempty(ops) && return QAdd(QMul{T}[QMul(arg_c, QSym[])])
    length(ops) == 1 && return QAdd(QMul{T}[QMul(arg_c, copy(ops))])

    for i in 1:(length(ops) - 1)
        a, b = ops[i], ops[i + 1]
        same_space = a.space_index == b.space_index && a.name == b.name

        # Fock: [a, a†] = 1
        if a isa Destroy && b isa Create && same_space
            swapped = QSym[ops[1:(i - 1)]..., b, a, ops[(i + 2):end]...]
            sort!(swapped; lt=canonical_lt)
            term1 = _normal_order_qmul(arg_c, swapped)
            contracted = QSym[ops[1:(i - 1)]..., ops[(i + 2):end]...]
            sort!(contracted; lt=canonical_lt)
            term2 = _normal_order_qmul(arg_c, contracted)
            return QAdd(QMul{T}[term1.arguments..., term2.arguments...])
        end

        # Transition: |i⟩⟨j| · |k⟩⟨l|
        if a isa Transition && b isa Transition && same_space
            if a.j == b.i
                # Composition: |a.i⟩⟨b.j|
                composed = Transition(a.name, a.i, b.j, a.space_index)
                new_ops = QSym[ops[1:(i - 1)]..., composed, ops[(i + 2):end]...]
                return _normal_order_qmul(arg_c, new_ops)
            else
                # Orthogonality: zero
                return QAdd(QMul{T}[QMul(zero(arg_c), QSym[])])
            end
        end

        # Pauli: σⱼ·σₖ = δⱼₖ + iϵⱼₖₗσₗ
        if a isa Pauli && b isa Pauli && same_space
            if a.axis == b.axis
                # σⱼ² = I
                new_ops = QSym[ops[1:(i - 1)]..., ops[(i + 2):end]...]
                return _normal_order_qmul(arg_c, new_ops)
            else
                # σⱼσₖ = iϵⱼₖₗσₗ
                eps = levicivita([a.axis, b.axis, 6 - a.axis - b.axis])
                l = 6 - a.axis - b.axis  # remaining axis: 1+2+3=6
                new_op = Pauli(a.name, l, a.space_index)
                new_ops = QSym[ops[1:(i - 1)]..., new_op, ops[(i + 2):end]...]
                new_c = arg_c * im * eps
                return _normal_order_qmul(new_c, new_ops)
            end
        end

        # Spin: [Sⱼ, Sₖ] = iϵⱼₖₗSₗ (swap out-of-order axes)
        if a isa Spin && b isa Spin && same_space && a.axis > b.axis
            # Swap: SⱼSₖ = SₖSⱼ + iϵⱼₖₗSₗ
            swapped = QSym[ops[1:(i - 1)]..., b, a, ops[(i + 2):end]...]
            sort!(swapped; lt=canonical_lt)
            term1 = _normal_order_qmul(arg_c, swapped)
            # Commutator term
            l = 6 - a.axis - b.axis
            eps = levicivita([a.axis, b.axis, l])
            comm_op = Spin(a.name, l, a.space_index)
            comm_ops = QSym[ops[1:(i - 1)]..., comm_op, ops[(i + 2):end]...]
            sort!(comm_ops; lt=canonical_lt)
            term2 = _normal_order_qmul(arg_c * im * eps, comm_ops)
            return QAdd(QMul{T}[term1.arguments..., term2.arguments...])
        end
    end

    # Already in normal order
    return QAdd(QMul{T}[QMul(arg_c, ops)])
end
```

Also update the docstring for `normal_order` and add the overload with spaces argument:

```julia
"""
    normal_order(expr)
    normal_order(expr, spaces)

Apply algebra rules to rewrite the expression in normal order.
The `spaces` argument (a `HilbertSpace` or mapping from space_index to
HilbertSpace) enables ground state rewriting for Transitions.
Always returns `QAdd`.
"""
function normal_order(m::QMul{T}, h::HilbertSpace) where {T}
    result = _normal_order_qmul(m.arg_c, m.args_nc)
    return _apply_ground_state(result, h)
end
function normal_order(s::QAdd{T}, h::HilbertSpace) where {T}
    result = normal_order(s)
    return _apply_ground_state(result, h)
end

function _apply_ground_state(expr::QAdd{T}, h::NLevelSpace) where {T}
    result = QMul{T}[]
    for term in expr.arguments
        expanded = _expand_ground_state(term, h)
        append!(result, expanded.arguments)
    end
    return QAdd(result)
end
_apply_ground_state(expr::QAdd, h::HilbertSpace) = expr  # no-op for non-NLevel

function _expand_ground_state(term::QMul{T}, h::NLevelSpace) where {T}
    g = h.ground_state
    n = h.n
    for (idx, op) in enumerate(term.args_nc)
        if op isa Transition && op.i == g && op.j == g
            # |g⟩⟨g| = 1 - Σ_{k≠g} |k⟩⟨k|
            terms = QMul{T}[]
            # Identity term (remove the transition)
            id_ops = QSym[term.args_nc[1:(idx - 1)]..., term.args_nc[(idx + 1):end]...]
            push!(terms, QMul(term.arg_c, id_ops))
            # Subtraction terms
            for k in 1:n
                k == g && continue
                new_op = Transition(op.name, k, k, op.space_index)
                sub_ops = QSym[term.args_nc[1:(idx - 1)]..., new_op, term.args_nc[(idx + 1):end]...]
                push!(terms, QMul(-term.arg_c, sub_ops))
            end
            return QAdd(terms)
        end
    end
    return QAdd(QMul{T}[term])
end
```

- [ ] **Step 3: Run tests**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: All tests pass including new normal_order tests.

---

### Task 6: Printing for New Types

**Files:**
- Modify: `src/printing.jl`
- Modify: `test/printing_test.jl`

- [ ] **Step 1: Add show methods to `src/printing.jl`**

Append to `src/printing.jl`:

```julia
# NLevel
Base.show(io::IO, h::NLevelSpace) = write(io, "ℋ(", string(h.name), ")")

# Pauli
Base.show(io::IO, h::PauliSpace) = write(io, "ℋ(", string(h.name), ")")

# Spin
Base.show(io::IO, h::SpinSpace) = write(io, "ℋ(", string(h.name), ")")

# Transition: σ₁₂
const _subscript_digits = ['₀', '₁', '₂', '₃', '₄', '₅', '₆', '₇', '₈', '₉']
function _to_subscript(n::Int)
    return join(_subscript_digits[d + 1] for d in reverse(digits(n)))
end

function Base.show(io::IO, x::Transition)
    write(io, string(x.name), _to_subscript(x.i), _to_subscript(x.j))
end

# Pauli: σx, σy, σz
const _xyz_sym = [:x, :y, :z]
function Base.show(io::IO, x::Pauli)
    write(io, string(x.name), string(_xyz_sym[x.axis]))
end

# Spin: Sx, Sy, Sz
function Base.show(io::IO, x::Spin)
    write(io, string(x.name), string(_xyz_sym[x.axis]))
end
```

- [ ] **Step 2: Add printing tests**

Append to `test/printing_test.jl`:

```julia
@testset "NLevel printing" begin
    h = NLevelSpace(:atom, 3, 1)
    σ12 = Transition(h, :σ, 1, 2)
    @test repr(σ12) == "σ₁₂"
    @test repr(h) == "ℋ(atom)"
end

@testset "Pauli printing" begin
    h = PauliSpace(:p)
    σx = Pauli(h, :σ, 1)
    σy = Pauli(h, :σ, 2)
    σz = Pauli(h, :σ, 3)
    @test repr(σx) == "σx"
    @test repr(σy) == "σy"
    @test repr(σz) == "σz"
    @test repr(h) == "ℋ(p)"
end

@testset "Spin printing" begin
    h = SpinSpace(:s, 1 // 2)
    Sx = Spin(h, :S, 1)
    Sy = Spin(h, :S, 2)
    Sz = Spin(h, :S, 3)
    @test repr(Sx) == "Sx"
    @test repr(Sy) == "Sy"
    @test repr(Sz) == "Sz"
    @test repr(h) == "ℋ(s)"
end
```

- [ ] **Step 3: Run tests**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: All tests pass.

---

### Task 7: LaTeX Recipes for New Types

**Files:**
- Modify: `src/latexify_recipes.jl`
- Modify: `test/latexify_test.jl`

- [ ] **Step 1: Add LaTeX recipes and transition_superscript toggle**

Add to `src/latexify_recipes.jl` (before the `QLaTeX` const):

```julia
# Transition superscript toggle
const transition_idx_script = Ref(:^)

"""
    transition_superscript(x::Bool)

Toggle whether Transition level indices are printed as superscript (true)
or subscript (false) in LaTeX. Default is superscript.
"""
function transition_superscript(x::Bool)
    transition_idx_script[] = x ? :^ : :_
    return x
end

@latexrecipe function f(x::Transition)
    return Expr(:latexifymerge, "{$(x.name)}$(transition_idx_script[]){{$(x.i)$(x.j)}}")
end

@latexrecipe function f(x::Pauli)
    ax = _xyz_sym[x.axis]
    return Expr(:latexifymerge, "{$(x.name)}_{{$ax}}")
end

@latexrecipe function f(x::Spin)
    ax = _xyz_sym[x.axis]
    return Expr(:latexifymerge, "{$(x.name)}_{{$ax}}")
end
```

Note: `_xyz_sym` is defined in `printing.jl` which is included before `latexify_recipes.jl`.

- [ ] **Step 2: Add LaTeX tests**

Append to `test/latexify_test.jl`:

```julia
@testset "Transition LaTeX" begin
    h = NLevelSpace(:atom, 3, 1)
    σ12 = Transition(h, :σ, 1, 2)
    s = latexify(σ12)
    @test occursin("σ", s)
    @test occursin("12", s)

    # Toggle superscript/subscript
    transition_superscript(false)
    s2 = latexify(σ12)
    @test occursin("_", s2)
    transition_superscript(true)  # reset
end

@testset "Pauli LaTeX" begin
    h = PauliSpace(:p)
    σx = Pauli(h, :σ, 1)
    s = latexify(σx)
    @test occursin("σ", s)
    @test occursin("x", s)
end

@testset "Spin LaTeX" begin
    h = SpinSpace(:s, 1 // 2)
    Sz = Spin(h, :S, 3)
    s = latexify(Sz)
    @test occursin("S", s)
    @test occursin("z", s)
end
```

- [ ] **Step 3: Run tests**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: All tests pass.

---

### Task 8: Numerical Conversion for New Types

**Files:**
- Modify: `src/numeric.jl`
- Modify: `test/numeric_test.jl`

- [ ] **Step 1: Add `to_numeric` methods to `src/numeric.jl`**

Add after the existing Fock methods:

```julia
# Transition
function to_numeric(op::Transition, b::QuantumOpticsBase.NLevelBasis; kwargs...)
    return QuantumOpticsBase.transition(b, op.i, op.j)
end

# Pauli
function to_numeric(op::Pauli, b::QuantumOpticsBase.SpinBasis; kwargs...)
    if op.axis == 1
        return QuantumOpticsBase.sigmax(b)
    elseif op.axis == 2
        return QuantumOpticsBase.sigmay(b)
    else
        return QuantumOpticsBase.sigmaz(b)
    end
end

# Spin
function to_numeric(op::Spin, b::QuantumOpticsBase.SpinBasis; kwargs...)
    if op.axis == 1
        return 0.5 * QuantumOpticsBase.sigmax(b)
    elseif op.axis == 2
        return 0.5 * QuantumOpticsBase.sigmay(b)
    else
        return 0.5 * QuantumOpticsBase.sigmaz(b)
    end
end
```

- [ ] **Step 2: Add numeric tests**

Append to `test/numeric_test.jl`:

```julia
@testset "NLevel numeric" begin
    h = NLevelSpace(:atom, 3, 1)
    σ12 = Transition(h, :σ, 1, 2)
    b = NLevelBasis(3)
    @test to_numeric(σ12, b) == transition(b, 1, 2)
    @test to_numeric(σ12', b) == transition(b, 2, 1)
end

@testset "Pauli numeric" begin
    h = PauliSpace(:p)
    σx = Pauli(h, :σ, 1)
    σy = Pauli(h, :σ, 2)
    σz = Pauli(h, :σ, 3)
    b = SpinBasis(1 // 2)
    @test to_numeric(σx, b) == sigmax(b)
    @test to_numeric(σy, b) == sigmay(b)
    @test to_numeric(σz, b) == sigmaz(b)
end

@testset "Spin numeric" begin
    h = SpinSpace(:s, 5 // 2)
    Sx = Spin(h, :S, 1)
    Sy = Spin(h, :S, 2)
    Sz = Spin(h, :S, 3)
    b = SpinBasis(5 // 2)
    @test to_numeric(Sx, b) == 0.5 * sigmax(b)
    @test to_numeric(Sy, b) == 0.5 * sigmay(b)
    @test to_numeric(Sz, b) == 0.5 * sigmaz(b)
end

@testset "Composite NLevel + Fock" begin
    h = FockSpace(:c) ⊗ NLevelSpace(:atom, 3, 1)
    @qnumbers a::Destroy(h, 1)
    σ12 = Transition(h, :σ, 1, 2, 2)
    bf = FockBasis(3)
    bn = NLevelBasis(3)
    bc = bf ⊗ bn
    @test to_numeric(σ12, bc) isa LazyTensor
    @test to_numeric(a' * σ12, bc) isa LazyTensor
end
```

- [ ] **Step 3: Run tests**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: All tests pass.

---

### Task 9: Update Module File and Exports

**Files:**
- Modify: `src/SecondQuantizedAlgebra.jl`
- Modify: `Project.toml` (re-add Combinatorics to deps if removed)

- [ ] **Step 1: Update `src/SecondQuantizedAlgebra.jl`**

The final module file should look like:

```julia
module SecondQuantizedAlgebra

using SymbolicUtils: SymbolicUtils
using Symbolics: Symbolics
using TermInterface: TermInterface

using Combinatorics: levicivita

using QuantumOpticsBase: QuantumOpticsBase
import QuantumOpticsBase: ⊗, tensor

using Latexify: Latexify, latexify, @latexrecipe

include("types.jl")
include("hilbertspace.jl")
include("fock.jl")
include("nlevel.jl")
include("pauli.jl")
include("spin.jl")
include("qmul.jl")
include("qadd.jl")
include("macros.jl")
include("interface.jl")
include("normal_order.jl")
include("simplify.jl")
include("printing.jl")
include("latexify_recipes.jl")
include("numeric.jl")

export QField, QSym, QTerm,
    HilbertSpace, FockSpace, ProductSpace,
    NLevelSpace, Transition,
    PauliSpace, Pauli,
    SpinSpace, Spin,
    ⊗, tensor,
    Destroy, Create,
    QMul, QAdd,
    @qnumbers,
    normal_order, simplify,
    to_numeric, numeric_average,
    transition_superscript

end
```

- [ ] **Step 2: Ensure Combinatorics is in `[deps]` of Project.toml**

Add if missing:
```toml
Combinatorics = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
```

- [ ] **Step 3: Run full test suite**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: ALL PASS

---

### Task 10: Integration Tests for Multi-Type Systems

**Files:**
- Modify: `test/integration_test.jl`

- [ ] **Step 1: Add multi-type integration tests**

Append to `test/integration_test.jl`:

```julia
@testset "Jaynes-Cummings model" begin
    # Fock + NLevel composite system
    hc = FockSpace(:cavity)
    ha = NLevelSpace(:atom, 2, 1)
    h = hc ⊗ ha
    @qnumbers a::Destroy(h, 1)
    σ12 = Transition(h, :σ, 1, 2, 2)
    σ21 = Transition(h, :σ, 2, 1, 2)

    # Interaction Hamiltonian: g(a†σ₁₂ + aσ₂₁)
    @variables g_jc
    H_int = g_jc * (a' * σ12 + a * σ21)
    @test H_int isa QAdd

    # Numeric check
    bf = FockBasis(5)
    bn = NLevelBasis(2)
    bc = bf ⊗ bn
    H_num = to_numeric(H_int, bc)
    @test H_num !== nothing  # just check it doesn't error
end

@testset "Pauli algebra via normal_order" begin
    h = PauliSpace(:p)
    σx = Pauli(h, :σ, 1)
    σy = Pauli(h, :σ, 2)
    σz = Pauli(h, :σ, 3)

    import SecondQuantizedAlgebra: simplify

    # σx·σy·σz = i·I (since σx·σy = iσz, then iσz·σz = i·I)
    m = σx * σy * σz
    result = simplify(normal_order(m))
    @test result isa QAdd
    @test length(result.arguments) == 1
    @test isempty(result.arguments[1].args_nc)
    @test result.arguments[1].arg_c == im
end

@testset "Mixed Fock + Pauli" begin
    h = FockSpace(:c) ⊗ PauliSpace(:p)
    @qnumbers a::Destroy(h, 1)
    σz = Pauli(h, :σ, 3, 2)

    # a†a ⊗ σz — operators on different spaces
    m = a' * a * σz
    @test m isa QMul
    @test length(m.args_nc) == 3

    # Normal order should commute a,a† on space 1, leave σz on space 2
    result = normal_order(a * a' * σz)
    @test result isa QAdd
end
```

- [ ] **Step 2: Run full test suite**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: ALL PASS
