# PhaseSpace / Position / Momentum Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add Position and Momentum quadrature operators on a PhaseSpace Hilbert space, with `[X, P] = i` commutation in simplify and numeric conversion via FockBasis.

**Architecture:** Two new concrete types (Position, Momentum) following the Destroy/Create pattern. PhaseSpace is a minimal HilbertSpace like PauliSpace. The commutation rule is an ordering-dependent swap in `_apply_ordering_rule` for NormalOrder: `P·X → X·P - i`.

**Tech Stack:** Julia, QuantumOpticsBase (FockBasis, destroy, create), Latexify (@latexrecipe), Combinatorics (not needed for this feature)

---

## File Structure

| File | Action | Responsibility |
|---|---|---|
| `src/phase_space.jl` | Create | PhaseSpace struct, Position/Momentum structs, constructors, adjoint, equality, hashing, ladder |
| `src/simplify.jl` | Modify | Add P·X → X·P - i ordering rule |
| `src/normal_order.jl` | Modify | Add `_apply_ground_state` no-op for PhaseSpace |
| `src/numeric.jl` | Modify | Add `to_numeric` for Position/Momentum on FockBasis |
| `src/printing.jl` | Modify | Add `show` methods for PhaseSpace, Position, Momentum |
| `src/latexify_recipes.jl` | Modify | Add `@latexrecipe` for Position, Momentum |
| `src/SecondQuantizedAlgebra.jl` | Modify | Add `include("phase_space.jl")` and exports |
| `test/phase_space_test.jl` | Create | All tests for the new types |

---

### Task 1: PhaseSpace, Position, Momentum structs + module wiring

**Files:**
- Create: `src/phase_space.jl`
- Modify: `src/SecondQuantizedAlgebra.jl`

- [ ] **Step 1: Create `src/phase_space.jl` with all type definitions**

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
end

"""
    Momentum <: QSym

Momentum (quadrature) operator on a [`PhaseSpace`](@ref).
"""
struct Momentum <: QSym
    name::Symbol
    space_index::Int
end

# Construction from Hilbert spaces
Position(h::PhaseSpace, name::Symbol) = Position(name, 1)
Momentum(h::PhaseSpace, name::Symbol) = Momentum(name, 1)

function Position(h::ProductSpace, name::Symbol, idx::Int)
    1 <= idx <= length(h.spaces) || throw(ArgumentError("Index $idx out of range for ProductSpace with $(length(h.spaces)) spaces"))
    h.spaces[idx] isa PhaseSpace || throw(ArgumentError("Space at index $idx is not a PhaseSpace"))
    return Position(name, idx)
end
function Momentum(h::ProductSpace, name::Symbol, idx::Int)
    1 <= idx <= length(h.spaces) || throw(ArgumentError("Index $idx out of range for ProductSpace with $(length(h.spaces)) spaces"))
    h.spaces[idx] isa PhaseSpace || throw(ArgumentError("Space at index $idx is not a PhaseSpace"))
    return Momentum(name, idx)
end

# Adjoint — Hermitian
Base.adjoint(op::Position) = op
Base.adjoint(op::Momentum) = op

# Equality
Base.isequal(a::Position, b::Position) = a.name == b.name && a.space_index == b.space_index
Base.isequal(a::Momentum, b::Momentum) = a.name == b.name && a.space_index == b.space_index
Base.:(==)(a::Position, b::Position) = isequal(a, b)
Base.:(==)(a::Momentum, b::Momentum) = isequal(a, b)

# Hashing
Base.hash(a::Position, h::UInt) = hash(:Position, hash(a.name, hash(a.space_index, h)))
Base.hash(a::Momentum, h::UInt) = hash(:Momentum, hash(a.name, hash(a.space_index, h)))

# Ladder (not applicable — Hermitian operators, no creation/annihilation distinction)
ladder(::Position) = 0
ladder(::Momentum) = 1
```

- [ ] **Step 2: Wire into module — add include and exports**

In `src/SecondQuantizedAlgebra.jl`, add `include("phase_space.jl")` after the `include("spin.jl")` line (line 19). The include block should read:

```julia
include("types.jl")
include("hilbertspace.jl")
include("fock.jl")
include("nlevel.jl")
include("pauli.jl")
include("spin.jl")
include("phase_space.jl")
include("qmul.jl")
include("qadd.jl")
include("macros.jl")
include("interface.jl")
include("simplify.jl")
include("normal_order.jl")
include("printing.jl")
include("latexify_recipes.jl")
include("numeric.jl")
```

Add to the export list (after `Spin`):

```julia
export QField, QSym, QTerm,
    HilbertSpace, FockSpace, ProductSpace,
    NLevelSpace, Transition,
    PauliSpace, Pauli,
    SpinSpace, Spin,
    PhaseSpace, Position, Momentum,
    ⊗, tensor,
    Destroy, Create,
    QMul, QAdd,
    @qnumbers,
    OrderingConvention, NormalOrder,
    normal_order, simplify,
    to_numeric, numeric_average,
    transition_superscript
```

---

### Task 2: Printing and LaTeX

**Files:**
- Modify: `src/printing.jl`
- Modify: `src/latexify_recipes.jl`

- [ ] **Step 1: Add REPL printing**

In `src/printing.jl`, add after the `Base.show(io::IO, h::SpinSpace)` line (line 15):

```julia
Base.show(io::IO, h::PhaseSpace) = write(io, "ℋ(", string(h.name), ")")
```

Add after the Pauli/Spin show methods (after line 33):

```julia
# Operators — Position/Momentum
Base.show(io::IO, x::Position) = write(io, string(x.name))
Base.show(io::IO, x::Momentum) = write(io, string(x.name))
```

- [ ] **Step 2: Add LaTeX recipes**

In `src/latexify_recipes.jl`, add after the Spin recipe (after line 35):

```julia
@latexrecipe function f(x::Position)
    return Expr(:latexifymerge, "\\hat{$(x.name)}")
end

@latexrecipe function f(x::Momentum)
    return Expr(:latexifymerge, "\\hat{$(x.name)}")
end
```

---

### Task 3: Simplify ordering rule

**Files:**
- Modify: `src/simplify.jl`

- [ ] **Step 1: Add the `[X, P] = i` ordering rule**

In `src/simplify.jl`, in the `_apply_ordering_rule` function for `NormalOrder` (around line 117-146), add a new block **after** the Spin block and **before** the final `return nothing`:

```julia
# PhaseSpace: [X, P] = i (swap P·X → X·P - i)
if a isa Momentum && b isa Position && same_space
    swapped = QSym[ops[1:(i - 1)]..., b, a, ops[(i + 2):end]...]
    sort!(swapped; lt=canonical_lt)
    term1 = _simplify_qmul(arg_c, swapped, NormalOrder())
    contracted = QSym[ops[1:(i - 1)]..., ops[(i + 2):end]...]
    sort!(contracted; lt=canonical_lt)
    term2 = _simplify_qmul(-im * arg_c, contracted, NormalOrder())
    all_terms = QMul[term1.arguments..., term2.arguments...]
    return _collect_qadd(all_terms)
end
```

---

### Task 4: Ground state no-op

**Files:**
- Modify: `src/normal_order.jl`

- [ ] **Step 1: Add PhaseSpace no-op**

In `src/normal_order.jl`, add after the `_apply_ground_state(expr::QAdd, ::SpinSpace) = expr` line (line 37):

```julia
_apply_ground_state(expr::QAdd, ::PhaseSpace) = expr
```

---

### Task 5: Numeric conversion

**Files:**
- Modify: `src/numeric.jl`

- [ ] **Step 1: Add `to_numeric` for Position and Momentum**

In `src/numeric.jl`, add after the Spin `to_numeric` block (after line 39):

```julia
# Position: X = (a + a†) / √2
function to_numeric(op::Position, b::QuantumOpticsBase.FockBasis; kwargs...)
    return (QuantumOpticsBase.destroy(b) + QuantumOpticsBase.create(b)) / sqrt(2)
end

# Momentum: P = i(a† - a) / √2
function to_numeric(op::Momentum, b::QuantumOpticsBase.FockBasis; kwargs...)
    return im * (QuantumOpticsBase.create(b) - QuantumOpticsBase.destroy(b)) / sqrt(2)
end
```

---

### Task 6: Tests

**Files:**
- Create: `test/phase_space_test.jl`

- [ ] **Step 1: Create the test file**

```julia
using SecondQuantizedAlgebra
using QuantumOpticsBase
using Latexify
using Test
import SecondQuantizedAlgebra: simplify

@testset "phase_space" begin
    @testset "PhaseSpace" begin
        h = PhaseSpace(:q)
        @test h isa HilbertSpace
        @test h.name == :q
        @test PhaseSpace(:q) == PhaseSpace(:q)
        @test PhaseSpace(:q) != PhaseSpace(:r)
        @test hash(PhaseSpace(:q)) == hash(PhaseSpace(:q))
    end

    @testset "Position construction — single space" begin
        h = PhaseSpace(:q)
        x = Position(h, :x)
        @test x isa Position
        @test x isa QSym
        @test x.name == :x
        @test x.space_index == 1
    end

    @testset "Momentum construction — single space" begin
        h = PhaseSpace(:q)
        p = Momentum(h, :p)
        @test p isa Momentum
        @test p isa QSym
        @test p.name == :p
        @test p.space_index == 1
    end

    @testset "Construction — product space" begin
        h = FockSpace(:c) ⊗ PhaseSpace(:q)
        x = Position(h, :x, 2)
        p = Momentum(h, :p, 2)
        @test x.space_index == 2
        @test p.space_index == 2
        @test_throws ArgumentError Position(h, :x, 1)  # FockSpace, not PhaseSpace
        @test_throws ArgumentError Momentum(h, :p, 1)
        @test_throws ArgumentError Position(h, :x, 3)  # out of range
    end

    @testset "Adjoint — Hermitian" begin
        h = PhaseSpace(:q)
        x = Position(h, :x)
        p = Momentum(h, :p)
        @test x' == x
        @test p' == p
    end

    @testset "Equality and hashing" begin
        h = PhaseSpace(:q)
        x1 = Position(h, :x)
        x2 = Position(h, :x)
        p1 = Momentum(h, :p)
        @test isequal(x1, x2)
        @test !isequal(x1, Position(h, :y))
        @test hash(x1) == hash(x2)
        @test hash(x1) != hash(p1)
    end

    @testset "Lazy multiplication" begin
        h = PhaseSpace(:q)
        x = Position(h, :x)
        p = Momentum(h, :p)

        m = x * p
        @test m isa QMul{Int}
        @test length(m.args_nc) == 2

        s = x + p
        @test s isa QAdd{Int}
    end

    @testset "Simplify: [X, P] = i" begin
        h = PhaseSpace(:q)
        x = Position(h, :x)
        p = Momentum(h, :p)

        # p*x out of order → x*p - i
        result = simplify(p * x)
        @test result isa QAdd
        @test length(result.arguments) == 2  # x·p + (-i)

        # x*p already ordered → x·p (1 term)
        result2 = simplify(x * p)
        @test result2 isa QAdd
        @test length(result2.arguments) == 1
    end

    @testset "Simplify: mixed spaces don't interact" begin
        h = PhaseSpace(:q1) ⊗ PhaseSpace(:q2)
        x1 = Position(h, :x, 1)
        p2 = Momentum(h, :p, 2)

        # Different spaces — no commutation
        result = simplify(p2 * x1)
        @test length(result.arguments) == 1
    end

    @testset "Numeric conversion" begin
        h = PhaseSpace(:q)
        x = Position(h, :x)
        p = Momentum(h, :p)
        b = FockBasis(10)

        x_num = to_numeric(x, b)
        p_num = to_numeric(p, b)
        @test x_num ≈ (destroy(b) + create(b)) / sqrt(2)
        @test p_num ≈ im * (create(b) - destroy(b)) / sqrt(2)
    end

    @testset "Numeric: commutator [X, P] = i" begin
        h = PhaseSpace(:q)
        x = Position(h, :x)
        p = Momentum(h, :p)
        b = FockBasis(30)

        x_num = to_numeric(x, b)
        p_num = to_numeric(p, b)
        comm = x_num * p_num - p_num * x_num
        # Should be approximately i*I (exact in infinite basis, approximate in truncated)
        @test comm.data[1, 1] ≈ im atol = 1e-10
    end

    @testset "Numeric: composite basis" begin
        h = FockSpace(:c) ⊗ PhaseSpace(:q)
        x = Position(h, :x, 2)
        bf = FockBasis(3)
        bq = FockBasis(5)
        bc = bf ⊗ bq
        @test to_numeric(x, bc) isa LazyTensor
    end

    @testset "Printing" begin
        h = PhaseSpace(:q)
        @test sprint(show, h) == "ℋ(q)"

        x = Position(h, :x)
        p = Momentum(h, :p)
        @test sprint(show, x) == "x"
        @test sprint(show, p) == "p"
    end

    @testset "LaTeX" begin
        h = PhaseSpace(:q)
        x = Position(h, :x)
        p = Momentum(h, :p)

        sx = latexify(x)
        sp = latexify(p)
        @test occursin("hat", sx)
        @test occursin("x", sx)
        @test occursin("hat", sp)
        @test occursin("p", sp)

        # MIME dispatch
        @test repr(MIME"text/latex"(), x) == sx
        @test repr(MIME"text/latex"(), p) == sp
    end

    @testset "@qnumbers" begin
        h = PhaseSpace(:q)
        @qnumbers x::Position(h)
        @qnumbers p::Momentum(h)
        @test x isa Position
        @test x.name == :x
        @test p isa Momentum
        @test p.name == :p
    end

    @testset "Concrete types" begin
        using CheckConcreteStructs
        all_concrete(PhaseSpace)
        all_concrete(Position)
        all_concrete(Momentum)
    end
end
```

- [ ] **Step 2: Run all tests**

Run: `make test`
Expected: All tests pass, including the new `phase_space_test.jl`.

---

## Self-Review

**Spec coverage:**
- PhaseSpace struct: Task 1, Step 1
- Position/Momentum structs: Task 1, Step 1
- Construction (single + product space): Task 1, Step 1
- Adjoint (Hermitian): Task 1, Step 1
- Equality/hashing: Task 1, Step 1
- Module wiring (include + exports): Task 1, Step 2
- Commutation `[X, P] = i`: Task 3, Step 1
- Ground state no-op: Task 4, Step 1
- Numeric conversion: Task 5, Step 1
- REPL printing: Task 2, Step 1
- LaTeX recipes: Task 2, Step 2
- Tests (construction, adjoint, equality, lazy mul, simplify, numeric, printing, LaTeX, mixed spaces): Task 6, Step 1

**Placeholder scan:** No TBD, TODO, or incomplete sections. All code blocks are complete.

**Type consistency:** Position/Momentum fields are `name::Symbol, space_index::Int` throughout. Constructor signatures match struct definitions. `_apply_ordering_rule` uses `a isa Momentum && b isa Position` consistently.
