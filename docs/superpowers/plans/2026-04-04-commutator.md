# Commutator Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement a type-stable `commutator(a, b)` function that computes `[a, b] = a*b - b*a`, simplifies, and always returns `QAdd{CT}`.

**Architecture:** Single file `src/commutator.jl` with 9 methods dispatching on operator types. Short-circuits for trivially commuting operands (scalars, different spaces, identical ops) return a pre-allocated `_ZERO_QADD` constant. Non-trivial cases delegate to `simplify(a*b - b*a)`. `QAdd` bilinearity distributes via `reduce(+, ...)` using existing promotion.

**Tech Stack:** Julia, SecondQuantizedAlgebra type system (`QField`, `QSym`, `QTerm`, `QMul`, `QAdd`), `simplify` from `src/simplify.jl`.

---

### Task 1: Create `src/commutator.jl` with core implementation

**Files:**
- Create: `src/commutator.jl`

- [ ] **Step 1: Create `src/commutator.jl` with all methods**

```julia
"""
    commutator(a, b)

Compute the commutator `[a, b] = a*b - b*a` and simplify the result.
Always returns `QAdd`. Short-circuits to zero when operands commute
(different Hilbert spaces, identical operators, or scalar arguments).
"""
function commutator end

const _ZERO_QADD = QAdd(QMul{Complex{Int}}[QMul(zero(Complex{Int}), QSym[])])

# Scalars commute with everything
commutator(::Number, ::Number) = _ZERO_QADD
commutator(::Number, ::QField) = _ZERO_QADD
commutator(::QField, ::Number) = _ZERO_QADD

# QSym, QSym: short-circuit when on different spaces or identical
function commutator(a::QSym, b::QSym)
    a.space_index == b.space_index || return _ZERO_QADD
    isequal(a, b) && return _ZERO_QADD
    return simplify(a * b - b * a)
end

# QMul, QSym: short-circuit when no shared space
function commutator(a::QMul, b::QSym)
    idx = findfirst(x -> x.space_index == b.space_index, a.args_nc)
    idx === nothing && return _ZERO_QADD
    return simplify(a * b - b * a)
end
function commutator(a::QSym, b::QMul)
    idx = findfirst(x -> x.space_index == a.space_index, b.args_nc)
    idx === nothing && return _ZERO_QADD
    return simplify(a * b - b * a)
end

# QMul, QMul: short-circuit when no shared spaces
function commutator(a::QMul, b::QMul)
    aon_a = map(x -> x.space_index, a.args_nc)
    aon_b = map(x -> x.space_index, b.args_nc)
    isempty(intersect(aon_a, aon_b)) && return _ZERO_QADD
    return simplify(a * b - b * a)
end

# QAdd: distribute (bilinearity)
function commutator(a::QAdd, b::QAdd)
    results = [commutator(a_, b_) for a_ in a.arguments, b_ in b.arguments]
    return reduce(+, results)
end
function commutator(a::QAdd, b::Union{QSym, QMul})
    results = [commutator(a_, b) for a_ in a.arguments]
    return reduce(+, results)
end
function commutator(a::Union{QSym, QMul}, b::QAdd)
    results = [commutator(a, b_) for b_ in b.arguments]
    return reduce(+, results)
end
```

### Task 2: Wire into module and export

**Files:**
- Modify: `src/SecondQuantizedAlgebra.jl`

- [ ] **Step 1: Add include after `simplify.jl`**

In `src/SecondQuantizedAlgebra.jl`, add the include line after line 23 (`include("simplify.jl")`):

```julia
include("commutator.jl")
```

- [ ] **Step 2: Add `commutator` to exports**

In the `export` block of `src/SecondQuantizedAlgebra.jl`, add `commutator` after `simplify`:

```julia
    normal_order, simplify, commutator,
```

### Task 3: Write tests

**Files:**
- Create: `test/commutator_test.jl`

- [ ] **Step 1: Create `test/commutator_test.jl`**

```julia
using SecondQuantizedAlgebra
using Test
import SecondQuantizedAlgebra: _ZERO_QADD, simplify

@testset "commutator" begin
    h = FockSpace(:c)
    a = Destroy(h, :a)
    ad = a'

    @testset "Scalar short-circuits" begin
        @test commutator(1, 2) ==== _ZERO_QADD
        @test commutator(3, a) ==== _ZERO_QADD
        @test commutator(a, 5) ==== _ZERO_QADD
        @test commutator(2.0, ad) ==== _ZERO_QADD
    end

    @testset "Fock: [a, a†] = 1" begin
        result = commutator(a, ad)
        @test result isa QAdd
        # simplify(a*a' - a'*a) = 1 (scalar identity)
        @test length(result.arguments) == 1
        @test isempty(result.arguments[1].args_nc)
        @test result.arguments[1].arg_c == 1
    end

    @testset "Fock: [a†, a] = -1" begin
        result = commutator(ad, a)
        @test result isa QAdd
        @test length(result.arguments) == 1
        @test result.arguments[1].arg_c == -1
    end

    @testset "Identical operators: [a, a] = 0" begin
        result = commutator(a, a)
        @test result === _ZERO_QADD
    end

    @testset "Different spaces: [a, b] = 0" begin
        h2 = FockSpace(:a) ⊗ FockSpace(:b)
        a1 = Destroy(h2, :a, 1)
        a2 = Destroy(h2, :b, 2)
        @test commutator(a1, a2) === _ZERO_QADD
        @test commutator(a1, a2') === _ZERO_QADD
    end

    @testset "QMul with QSym — shared space" begin
        result = commutator(ad * a, a)
        @test result isa QAdd
    end

    @testset "QMul with QSym — no shared space" begin
        h2 = FockSpace(:a) ⊗ FockSpace(:b)
        a1 = Destroy(h2, :a, 1)
        a2 = Destroy(h2, :b, 2)
        @test commutator(a1' * a1, a2) === _ZERO_QADD
        @test commutator(a2, a1' * a1) === _ZERO_QADD
    end

    @testset "QMul with QMul — no shared space" begin
        h2 = FockSpace(:a) ⊗ FockSpace(:b)
        a1 = Destroy(h2, :a, 1)
        a2 = Destroy(h2, :b, 2)
        @test commutator(a1' * a1, a2' * a2) === _ZERO_QADD
    end

    @testset "QAdd bilinearity: [a + a†, a] = [a,a] + [a†,a]" begin
        expr = a + ad
        result = commutator(expr, a)
        @test result isa QAdd
        # [a,a] = 0, [a†,a] = -1 → total = -1
        simplified = simplify(result)
        @test length(simplified.arguments) == 1
        @test simplified.arguments[1].arg_c == -1
    end

    @testset "QAdd, QAdd bilinearity" begin
        result = commutator(a + ad, a + ad)
        @test result isa QAdd
        # [a+a†, a+a†] = [a,a]+[a,a†]+[a†,a]+[a†,a†] = 0+1-1+0 = 0
        simplified = simplify(result)
        @test all(iszero, simplified.arguments)
    end

    @testset "Nested commutator" begin
        # [a, [a, a†]] = [a, 1] = 0
        inner = commutator(a, ad)
        # inner is QAdd representing 1
        # commutator with QAdd dispatches to bilinearity
        # each term is a scalar QMul (no ops) → commutator(a, scalar_qmul)
        # QSym with QMul that has no shared space → _ZERO_QADD
        result = commutator(a, inner)
        @test result isa QAdd
        @test all(iszero, simplify(result).arguments)
    end

    @testset "Spin: [Sx, Sy] = iSz" begin
        hs = SpinSpace(:s, 1 // 2)
        Sx = Spin(hs, :S, 1)
        Sy = Spin(hs, :S, 2)
        Sz = Spin(hs, :S, 3)
        result = commutator(Sx, Sy)
        @test result isa QAdd
        # [Sx, Sy] = iSz → simplify(SxSy - SySx)
        # SxSy is already ordered, SySx needs swap: SySx = SxSy - iSz
        # So SxSy - (SxSy - iSz) = iSz
        simplified = simplify(result)
        @test length(simplified.arguments) == 1
        @test simplified.arguments[1].arg_c == im
        @test simplified.arguments[1].args_nc[1] isa Spin
        @test simplified.arguments[1].args_nc[1].axis == 3
    end

    @testset "PhaseSpace: [X, P] = i" begin
        hps = PhaseSpace(:q)
        x = Position(hps, :x)
        p = Momentum(hps, :p)
        result = commutator(x, p)
        @test result isa QAdd
        simplified = simplify(result)
        @test length(simplified.arguments) == 1
        @test isempty(simplified.arguments[1].args_nc)
        @test simplified.arguments[1].arg_c == im
    end

    @testset "Return type is always QAdd" begin
        @test commutator(a, ad) isa QAdd
        @test commutator(a, a) isa QAdd
        @test commutator(1, a) isa QAdd
        @test commutator(a + ad, a) isa QAdd
        @test commutator(a + ad, a + ad) isa QAdd
        @test commutator(ad * a, a) isa QAdd
    end

    @testset "Type stability" begin
        @inferred commutator(a, ad)
        @inferred commutator(ad, a)
        @inferred commutator(a, a)
        @inferred commutator(1, a)
        @inferred commutator(a, 1)
        @inferred commutator(ad * a, a)
        @inferred commutator(a, ad * a)
        @inferred commutator(ad * a, a * ad)
        @inferred commutator(a + ad, a)
        @inferred commutator(a, a + ad)
        @inferred commutator(a + ad, a + ad)

        h2 = FockSpace(:a) ⊗ FockSpace(:b)
        a1 = Destroy(h2, :a, 1)
        a2 = Destroy(h2, :b, 2)
        @inferred commutator(a1, a2)
        @inferred commutator(a1, a2')
    end

    @testset "Allocations" begin
        # Warmup
        commutator(a, ad)
        commutator(a, a)
        commutator(ad * a, a)
        commutator(a + ad, a)

        # Short-circuit (zero) should be very cheap
        @test @allocations(commutator(a, a)) == 0
        @test @allocations(commutator(1, a)) == 0

        # Non-trivial but simple
        @test @allocations(commutator(a, ad)) < 100
        @test @allocations(commutator(ad * a, a)) < 200

        # Bilinearity over QAdd
        @test @allocations(commutator(a + ad, a)) < 200
    end
end
```

### Task 4: Verify tests pass

- [ ] **Step 1: Run commutator tests**

Run: `julia --project -e 'using Pkg; Pkg.test(test_args=["commutator_test"])'`

Expected: All tests PASS.

- [ ] **Step 2: Run full test suite**

Run: `julia --project -e 'using Pkg; Pkg.test()'`

Expected: All tests PASS (including existing tests unchanged).
