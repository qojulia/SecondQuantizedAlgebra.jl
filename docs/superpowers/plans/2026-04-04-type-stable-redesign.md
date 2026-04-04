# Type-Stable Redesign Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the type-unstable SecondQuantizedAlgebra.jl internals with a concrete, static, type-stable architecture for Fock space operators, following KeldyshContraction.jl patterns.

**Architecture:** Old code is preserved in `src/legacy/` and `test/legacy/` for reference. New code is built in `src/` with fresh test files in `test/`. The module entry point `src/SecondQuantizedAlgebra.jl` includes only new files. Each task builds one layer, tested before moving to the next.

**Tech Stack:** Julia, SymbolicUtils.jl, Symbolics.jl, TermInterface.jl, QuantumOpticsBase.jl, Latexify.jl, MacroTools.jl

**Reference files:** Old implementations in `src/legacy/` — especially `qnumber.jl`, `fock.jl`, `hilbertspace.jl`, `utils.jl`, `printing.jl`, `latexify_recipes.jl`. KeldyshContraction.jl at `/tmp/KeldyshContraction.jl/` for design patterns.

---

### Task 1: Abstract Type Hierarchy

**Files:**
- Create: `src/types.jl`
- Modify: `src/SecondQuantizedAlgebra.jl`
- Test: `test/types_test.jl`

- [ ] **Step 1: Write the test file**

```julia
# test/types_test.jl
using SecondQuantizedAlgebra
using Test

@testset "type hierarchy" begin
    @test QSym <: QField
    @test QTerm <: QField
    @test !(QSym <: QTerm)
    @test !(QTerm <: QSym)
    @test isabstracttype(QField)
    @test isabstracttype(QSym)
    @test isabstracttype(QTerm)
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project -e 'using Pkg; Pkg.test(test_args=["types_test"])'`
Expected: FAIL — `QField` not defined

- [ ] **Step 3: Write the implementation**

```julia
# src/types.jl
"""
    QField

Abstract type representing any expression involving operators.
"""
abstract type QField end

"""
    QSym <: QField

Abstract type representing fundamental operator types (leaves in the expression tree).
"""
abstract type QSym <: QField end

"""
    QTerm <: QField

Abstract type representing compound noncommutative expressions.
"""
abstract type QTerm <: QField end
```

Update `src/SecondQuantizedAlgebra.jl` to include and export:

```julia
module SecondQuantizedAlgebra

using SymbolicUtils: SymbolicUtils, BasicSymbolic, arguments, iscall, operation, substitute
using Symbolics: Symbolics
using TermInterface: TermInterface

using Combinatorics: combinations

using SciMLBase: SciMLBase
using QuantumOpticsBase: QuantumOpticsBase
import QuantumOpticsBase: ⊗, tensor

using LaTeXStrings: LaTeXStrings, @L_str, latexstring
using Latexify: Latexify, latexify, @latexrecipe
using MacroTools: MacroTools

include("types.jl")

export QField, QSym, QTerm

end
```

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project -e 'using Pkg; Pkg.test(test_args=["types_test"])'`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/types.jl src/SecondQuantizedAlgebra.jl test/types_test.jl
git commit -m "feat: add abstract type hierarchy (QField, QSym, QTerm)"
```

---

### Task 2: Hilbert Spaces

**Files:**
- Create: `src/hilbertspace.jl`
- Modify: `src/SecondQuantizedAlgebra.jl`
- Test: `test/hilbertspace_test.jl`
- Reference: `src/legacy/hilbertspace.jl`

- [ ] **Step 1: Write the test file**

```julia
# test/hilbertspace_test.jl
using SecondQuantizedAlgebra
using Test

@testset "hilbert spaces" begin
    @testset "FockSpace" begin
        h1 = FockSpace(:a)
        h2 = FockSpace(:b)
        h3 = FockSpace(:a)

        @test h1 == h3
        @test h1 != h2
        @test h1.name == :a
        @test h1 isa HilbertSpace
    end

    @testset "ProductSpace" begin
        h1 = FockSpace(:a)
        h2 = FockSpace(:b)
        h3 = FockSpace(:c)
        h4 = FockSpace(:d)

        h12 = h1 ⊗ h2
        @test h12 isa ProductSpace
        @test h12.spaces == (h1, h2)

        # Associativity
        h123_a = (h1 ⊗ h2) ⊗ h3
        h123_b = h1 ⊗ (h2 ⊗ h3)
        h123_c = h1 ⊗ h2 ⊗ h3
        @test h123_a == h123_b == h123_c
        @test h123_a.spaces == (h1, h2, h3)

        # 4 spaces
        h1234 = h1 ⊗ h2 ⊗ h3 ⊗ h4
        @test h1234.spaces == (h1, h2, h3, h4)
        @test (h1 ⊗ h2) ⊗ (h3 ⊗ h4) == h1234

        # tensor alias
        @test tensor(h1, h2, h3, h4) == h1234
    end

    @testset "Concrete types" begin
        using CheckConcreteStructs
        all_concrete(FockSpace)
        all_concrete(ProductSpace)
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project -e 'using Pkg; Pkg.test(test_args=["hilbertspace_test"])'`
Expected: FAIL — `FockSpace` not defined

- [ ] **Step 3: Write the implementation**

```julia
# src/hilbertspace.jl

"""
    HilbertSpace

Abstract type for representing Hilbert spaces. Used at construction time
for operator validation — not stored on operators at runtime.
"""
abstract type HilbertSpace end

"""
    FockSpace <: HilbertSpace

Hilbert space for bosonic operators (quantum harmonic oscillator).
"""
struct FockSpace <: HilbertSpace
    name::Symbol
end
Base.:(==)(a::FockSpace, b::FockSpace) = a.name == b.name
Base.hash(a::FockSpace, h::UInt) = hash(:FockSpace, hash(a.name, h))

"""
    ProductSpace{T} <: HilbertSpace

Composite Hilbert space consisting of multiple subspaces.
Uses a Tuple for fully concrete storage.
"""
struct ProductSpace{T<:Tuple{Vararg{HilbertSpace}}} <: HilbertSpace
    spaces::T
end
function Base.:(==)(a::ProductSpace, b::ProductSpace)
    return a.spaces == b.spaces
end
function Base.hash(a::ProductSpace, h::UInt)
    return hash(:ProductSpace, hash(a.spaces, h))
end

"""
    ⊗(spaces::HilbertSpace...)

Create a [`ProductSpace`](@ref) consisting of multiple subspaces.
Unicode `\\otimes<tab>`.
"""
⊗(a::FockSpace, b::FockSpace) = ProductSpace((a, b))
⊗(a::ProductSpace, b::FockSpace) = ProductSpace((a.spaces..., b))
⊗(a::FockSpace, b::ProductSpace) = ProductSpace((a, b.spaces...))
⊗(a::ProductSpace, b::ProductSpace) = ProductSpace((a.spaces..., b.spaces...))
⊗(a::HilbertSpace, b::HilbertSpace, c::HilbertSpace...) = ⊗(a ⊗ b, c...)
⊗(a::HilbertSpace) = a

"""
    tensor(spaces::HilbertSpace...)

Create a [`ProductSpace`](@ref). Alias for [`⊗`](@ref).
"""
tensor(args::Vararg{HilbertSpace}) = ⊗(args...)
```

Update `src/SecondQuantizedAlgebra.jl`:

```julia
include("types.jl")
include("hilbertspace.jl")

export QField, QSym, QTerm,
    HilbertSpace, FockSpace, ProductSpace,
    ⊗, tensor
```

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project -e 'using Pkg; Pkg.test(test_args=["hilbertspace_test"])'`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/hilbertspace.jl src/SecondQuantizedAlgebra.jl test/hilbertspace_test.jl
git commit -m "feat: add FockSpace and ProductSpace with concrete Tuple storage"
```

---

### Task 3: Destroy and Create Operators

**Files:**
- Create: `src/fock.jl`
- Modify: `src/SecondQuantizedAlgebra.jl`
- Test: `test/fock_test.jl`
- Reference: `src/legacy/fock.jl`, `/tmp/KeldyshContraction.jl/src/keldysh_algebra/keldysh_algebra.jl`

- [ ] **Step 1: Write the test file**

```julia
# test/fock_test.jl
using SecondQuantizedAlgebra
using Test

@testset "fock operators" begin
    @testset "Construction — single space" begin
        h = FockSpace(:c)
        a = Destroy(h, :a)
        @test a isa Destroy
        @test a isa QSym
        @test a.name == :a
        @test a.space_index == 1
    end

    @testset "Construction — product space" begin
        h = FockSpace(:a) ⊗ FockSpace(:b)
        a = Destroy(h, :a, 1)
        b = Destroy(h, :b, 2)
        @test a.space_index == 1
        @test b.space_index == 2
        @test_throws AssertionError Destroy(h, :c, 3)
    end

    @testset "Adjoint" begin
        h = FockSpace(:c)
        a = Destroy(h, :a)
        ad = a'
        @test ad isa Create
        @test ad.name == :a
        @test ad.space_index == 1
        @test ad' == a
        @test a'' == a
    end

    @testset "Equality and hashing" begin
        h = FockSpace(:c)
        a1 = Destroy(h, :a)
        a2 = Destroy(h, :a)
        b = Destroy(h, :b)
        @test isequal(a1, a2)
        @test !isequal(a1, b)
        @test hash(a1) == hash(a2)
        @test hash(a1) != hash(b)
        @test hash(a1) != hash(a1')
    end

    @testset "Canonical ordering helpers" begin
        h = FockSpace(:c)
        a = Destroy(h, :a)
        import SecondQuantizedAlgebra: ladder
        @test ladder(a) == 1
        @test ladder(a') == 0
    end

    @testset "Concrete types" begin
        using CheckConcreteStructs
        all_concrete(Destroy)
        all_concrete(Create)
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project -e 'using Pkg; Pkg.test(test_args=["fock_test"])'`
Expected: FAIL — `Destroy` not defined

- [ ] **Step 3: Write the implementation**

```julia
# src/fock.jl

"""
    Destroy <: QSym

Bosonic annihilation operator. Stores only name and space index.
"""
struct Destroy <: QSym
    name::Symbol
    space_index::Int
end

"""
    Create <: QSym

Bosonic creation operator. Stores only name and space index.
"""
struct Create <: QSym
    name::Symbol
    space_index::Int
end

# Construction from Hilbert spaces (validation, then discard)
Destroy(h::FockSpace, name::Symbol) = Destroy(name, 1)
Create(h::FockSpace, name::Symbol) = Create(name, 1)

function Destroy(h::ProductSpace, name::Symbol, idx::Int)
    @assert 1 <= idx <= length(h.spaces) "Index $idx out of range for ProductSpace with $(length(h.spaces)) spaces"
    @assert h.spaces[idx] isa FockSpace "Space at index $idx is not a FockSpace"
    return Destroy(name, idx)
end
function Create(h::ProductSpace, name::Symbol, idx::Int)
    @assert 1 <= idx <= length(h.spaces) "Index $idx out of range for ProductSpace with $(length(h.spaces)) spaces"
    @assert h.spaces[idx] isa FockSpace "Space at index $idx is not a FockSpace"
    return Create(name, idx)
end

# Adjoint
Base.adjoint(op::Destroy) = Create(op.name, op.space_index)
Base.adjoint(op::Create) = Destroy(op.name, op.space_index)

# Equality
Base.isequal(a::Destroy, b::Destroy) = a.name == b.name && a.space_index == b.space_index
Base.isequal(a::Create, b::Create) = a.name == b.name && a.space_index == b.space_index
Base.:(==)(a::Destroy, b::Destroy) = isequal(a, b)
Base.:(==)(a::Create, b::Create) = isequal(a, b)

# Hashing
Base.hash(a::Destroy, h::UInt) = hash(:Destroy, hash(a.name, hash(a.space_index, h)))
Base.hash(a::Create, h::UInt) = hash(:Create, hash(a.name, hash(a.space_index, h)))

# Canonical ordering: Create (0) before Destroy (1), then by space_index, then name
"""
    ladder(op::QSym)

Returns 0 for creation operators, 1 for annihilation operators.
Used for canonical ordering within `QMul.args_nc`.
"""
ladder(::Create) = 0
ladder(::Destroy) = 1

"""
    canonical_lt(a::QSym, b::QSym)

Canonical ordering comparator: sort by (ladder, space_index, name).
"""
function canonical_lt(a::QSym, b::QSym)
    la, lb = ladder(a), ladder(b)
    la != lb && return la < lb
    a.space_index != b.space_index && return a.space_index < b.space_index
    return a.name < b.name
end
```

Update `src/SecondQuantizedAlgebra.jl`:

```julia
include("types.jl")
include("hilbertspace.jl")
include("fock.jl")

export QField, QSym, QTerm,
    HilbertSpace, FockSpace, ProductSpace,
    ⊗, tensor,
    Destroy, Create
```

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project -e 'using Pkg; Pkg.test(test_args=["fock_test"])'`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/fock.jl src/SecondQuantizedAlgebra.jl test/fock_test.jl
git commit -m "feat: add Destroy and Create operators with concrete fields"
```

---

### Task 4: QMul — Lazy Multiplication

**Files:**
- Create: `src/qmul.jl`
- Modify: `src/SecondQuantizedAlgebra.jl`
- Test: `test/qmul_test.jl`
- Reference: `/tmp/KeldyshContraction.jl/src/keldysh_algebra/QTerm.jl`

- [ ] **Step 1: Write the test file**

```julia
# test/qmul_test.jl
using SecondQuantizedAlgebra
using Test

@testset "QMul" begin
    h = FockSpace(:c)
    a = Destroy(h, :a)
    ad = a'

    @testset "Construction" begin
        m = a * ad
        @test m isa QMul{Int}
        @test m.arg_c == 1
        @test length(m.args_nc) == 2
        # Canonical order: Create before Destroy
        @test m.args_nc[1] isa Create
        @test m.args_nc[2] isa Destroy
    end

    @testset "Scalar multiplication" begin
        m1 = 3 * a
        @test m1 isa QMul{Int}
        @test m1.arg_c == 3
        @test m1.args_nc == [a]

        m2 = a * 2.0
        @test m2 isa QMul{Float64}
        @test m2.arg_c == 2.0

        m3 = 0 * a
        @test m3 isa QMul{Int}
        @test m3.arg_c == 0
    end

    @testset "QSym * QSym" begin
        h2 = FockSpace(:a) ⊗ FockSpace(:b)
        a1 = Destroy(h2, :a, 1)
        a2 = Destroy(h2, :b, 2)
        m = a1 * a2
        @test m isa QMul{Int}
        # Both are Destroy (ladder=1), so sorted by space_index
        @test m.args_nc[1].space_index <= m.args_nc[2].space_index
    end

    @testset "QMul * QSym and QSym * QMul" begin
        m1 = (2 * a) * ad
        @test m1 isa QMul{Int}
        @test m1.arg_c == 2
        @test length(m1.args_nc) == 2

        m2 = ad * (3 * a)
        @test m2 isa QMul{Int}
        @test m2.arg_c == 3
        @test length(m2.args_nc) == 2
    end

    @testset "QMul * QMul" begin
        m1 = 2 * a
        m2 = 3 * ad
        m3 = m1 * m2
        @test m3 isa QMul{Int}
        @test m3.arg_c == 6
        @test length(m3.args_nc) == 2
    end

    @testset "QMul * Number" begin
        m = (2 * a) * 3
        @test m isa QMul{Int}
        @test m.arg_c == 6
    end

    @testset "Division" begin
        m = a / 2
        @test m isa QMul
        @test m.arg_c == 1 // 2
    end

    @testset "Power" begin
        m = a^3
        @test m isa QMul
        @test length(m.args_nc) == 3
    end

    @testset "Negation" begin
        m = -a
        @test m isa QMul{Int}
        @test m.arg_c == -1
    end

    @testset "Lazy — no commutation" begin
        # a * a' should NOT expand to a'a + 1
        m = a * ad
        @test m isa QMul  # not QAdd
        @test length(m.args_nc) == 2
    end

    @testset "Equality and hashing" begin
        m1 = 2 * a * ad
        m2 = 2 * a * ad
        @test isequal(m1, m2)
        @test hash(m1) == hash(m2)
    end

    @testset "iszero" begin
        m = 0 * a
        @test iszero(m)
        @test !iszero(2 * a)
    end

    @testset "Concrete types" begin
        using CheckConcreteStructs
        all_concrete(QMul)
    end

    @testset "Type stability" begin
        @inferred a * ad
        @inferred 3 * a
        @inferred a * 2.0
        @inferred (2 * a) * ad
        @inferred (2 * a) * (3 * ad)
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project -e 'using Pkg; Pkg.test(test_args=["qmul_test"])'`
Expected: FAIL — `QMul` not defined

- [ ] **Step 3: Write the implementation**

```julia
# src/qmul.jl

"""
    QMul{T<:Number} <: QTerm

Lazy product of quantum operators with a commutative prefactor.
`*` never applies commutation relations — it just collects operators
in canonical order and multiplies prefactors.

Fields:
- `arg_c::T` — commutative (c-number) prefactor
- `args_nc::Vector{QSym}` — non-commutative operator factors, canonically sorted
"""
struct QMul{T<:Number} <: QTerm
    arg_c::T
    args_nc::Vector{QSym}
    function QMul(arg_c::T, args_nc::Vector{QSym}) where {T<:Number}
        return new{T}(arg_c, args_nc)
    end
    QMul(args_nc::Vector{QSym}) = new{Int}(1, args_nc)
end

Base.length(a::QMul) = length(a.args_nc)
Base.iszero(a::QMul) = iszero(a.arg_c)
Base.zero(::QMul{T}) where {T} = QMul(zero(T), QSym[])
Base.zero(::Type{QMul{T}}) where {T} = QMul(zero(T), QSym[])

# Equality and hashing
function Base.isequal(a::QMul, b::QMul)
    isequal(a.arg_c, b.arg_c) || return false
    length(a.args_nc) == length(b.args_nc) || return false
    for (x, y) in zip(a.args_nc, b.args_nc)
        isequal(x, y) || return false
    end
    return true
end
Base.hash(q::QMul, h::UInt) = hash(:QMul, hash(q.arg_c, hash(q.args_nc, h)))

# Adjoint
function Base.adjoint(q::QMul)
    args_nc = QSym[adjoint(op) for op in reverse(q.args_nc)]
    sort!(args_nc; lt=canonical_lt)
    return QMul(conj(q.arg_c), args_nc)
end

# Promote rules
Base.promote_rule(::Type{QMul{S}}, ::Type{QMul{T}}) where {S,T} = QMul{promote_type(S, T)}
function Base.convert(::Type{QMul{T}}, x::QMul{S}) where {T<:Number,S<:Number}
    return QMul(convert(T, x.arg_c), x.args_nc)
end

## Multiplication — always returns QMul

# QSym * QSym → QMul{Int}
function Base.:*(a::QSym, b::QSym)
    args = QSym[a, b]
    sort!(args; lt=canonical_lt)
    return QMul(1, args)
end

# QSym * Number → QMul
Base.:*(a::QSym, b::Number) = QMul(b, QSym[a])
Base.:*(b::Number, a::QSym) = QMul(b, QSym[a])

# QMul * Number → QMul
function Base.:*(a::QMul, b::Number)
    return QMul(a.arg_c * b, a.args_nc)
end
Base.:*(b::Number, a::QMul) = a * b

# QSym * QMul → QMul
function Base.:*(a::QSym, b::QMul)
    args_nc = QSym[a, b.args_nc...]
    sort!(args_nc; lt=canonical_lt)
    return QMul(b.arg_c, args_nc)
end
function Base.:*(a::QMul, b::QSym)
    args_nc = QSym[a.args_nc..., b]
    sort!(args_nc; lt=canonical_lt)
    return QMul(a.arg_c, args_nc)
end

# QMul * QMul → QMul
function Base.:*(a::QMul, b::QMul)
    args_nc = QSym[a.args_nc..., b.args_nc...]
    sort!(args_nc; lt=canonical_lt)
    return QMul(a.arg_c * b.arg_c, args_nc)
end

# Division
Base.:/(a::QSym, b::Number) = QMul(1 // b, QSym[a])
Base.:/(a::QMul, b::Number) = QMul(a.arg_c // b, a.args_nc)

# Power
function Base.:^(a::QSym, n::Integer)
    @assert n >= 0 "Negative powers not supported"
    n == 0 && return QMul(1, QSym[])
    args_nc = QSym[a for _ in 1:n]
    return QMul(1, args_nc)
end
function Base.:^(a::QMul, n::Integer)
    @assert n >= 0 "Negative powers not supported"
    n == 0 && return QMul(1, QSym[])
    result = a
    for _ in 2:n
        result = result * a
    end
    return result
end

# Negation
Base.:-(a::QSym) = QMul(-1, QSym[a])
Base.:-(a::QMul) = QMul(-a.arg_c, a.args_nc)
```

Update `src/SecondQuantizedAlgebra.jl` — add `include("qmul.jl")` and `export QMul`.

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project -e 'using Pkg; Pkg.test(test_args=["qmul_test"])'`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/qmul.jl src/SecondQuantizedAlgebra.jl test/qmul_test.jl
git commit -m "feat: add QMul with lazy multiplication, canonical ordering"
```

---

### Task 5: QAdd — Addition

**Files:**
- Create: `src/qadd.jl`
- Modify: `src/SecondQuantizedAlgebra.jl`
- Test: `test/qadd_test.jl`
- Reference: `/tmp/KeldyshContraction.jl/src/keldysh_algebra/QTerm.jl`, `/tmp/KeldyshContraction.jl/src/keldysh_algebra/field_math.jl`

- [ ] **Step 1: Write the test file**

```julia
# test/qadd_test.jl
using SecondQuantizedAlgebra
using Test

@testset "QAdd" begin
    h = FockSpace(:c)
    a = Destroy(h, :a)
    ad = a'

    @testset "QSym + QSym" begin
        s = a + ad
        @test s isa QAdd{Int}
        @test length(s.arguments) == 2
        # Each QSym is wrapped in QMul
        @test all(x -> x isa QMul{Int}, s.arguments)
    end

    @testset "QMul + QMul" begin
        s = (2 * a) + (3 * ad)
        @test s isa QAdd{Int}
        @test length(s.arguments) == 2
    end

    @testset "QMul + QSym" begin
        s1 = (2 * a) + ad
        @test s1 isa QAdd{Int}
        @test length(s1.arguments) == 2

        s2 = ad + (2 * a)
        @test s2 isa QAdd{Int}
    end

    @testset "QAdd + QMul" begin
        s = (a + ad) + (3 * a)
        @test s isa QAdd{Int}
        @test length(s.arguments) == 3
    end

    @testset "QAdd + QAdd" begin
        s = (a + ad) + (a + ad)
        @test s isa QAdd{Int}
        @test length(s.arguments) == 4
    end

    @testset "QField + Number" begin
        s1 = a + 5
        @test s1 isa QAdd{Int}
        @test length(s1.arguments) == 2

        s2 = 5 + a
        @test s2 isa QAdd{Int}

        s3 = (a + ad) + 3
        @test s3 isa QAdd{Int}
        @test length(s3.arguments) == 3

        s4 = 3 + (a + ad)
        @test s4 isa QAdd{Int}
    end

    @testset "Subtraction" begin
        s = a - ad
        @test s isa QAdd
        @test length(s.arguments) == 2

        s2 = a - 3
        @test s2 isa QAdd
    end

    @testset "QAdd * Number" begin
        s = (a + ad) * 3
        @test s isa QAdd{Int}
        @test all(x -> x.arg_c == 3, s.arguments)
    end

    @testset "QAdd * QSym (distributes)" begin
        s = (a + ad) * a
        @test s isa QAdd{Int}
        @test length(s.arguments) == 2
        @test all(x -> length(x.args_nc) == 2, s.arguments)
    end

    @testset "QAdd * QMul (distributes)" begin
        s = (a + ad) * (2 * a)
        @test s isa QAdd{Int}
        @test length(s.arguments) == 2
    end

    @testset "QAdd * QAdd (distributes)" begin
        s = (a + ad) * (a + ad)
        @test s isa QAdd{Int}
        @test length(s.arguments) == 4
    end

    @testset "Equality and hashing" begin
        s1 = a + ad
        s2 = a + ad
        @test isequal(s1, s2)
        @test hash(s1) == hash(s2)
    end

    @testset "Adjoint" begin
        s = a + 2 * ad
        sd = s'
        @test sd isa QAdd
        @test length(sd.arguments) == 2
    end

    @testset "Concrete types" begin
        using CheckConcreteStructs
        all_concrete(QAdd)
    end

    @testset "Type stability" begin
        @inferred a + ad
        @inferred (2 * a) + (3 * ad)
        @inferred (a + ad) + a
        @inferred (a + ad) * 3
        @inferred (a + ad) * a
        @inferred (a + ad) * (a + ad)
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project -e 'using Pkg; Pkg.test(test_args=["qadd_test"])'`
Expected: FAIL — `QAdd` not defined

- [ ] **Step 3: Write the implementation**

```julia
# src/qadd.jl

"""
    QAdd{T<:Number} <: QTerm

Sum of [`QMul{T}`](@ref) terms.

Fields:
- `arguments::Vector{QMul{T}}` — terms of the sum
"""
struct QAdd{T<:Number} <: QTerm
    arguments::Vector{QMul{T}}
    function QAdd(args::Vector{QMul{T}}) where {T<:Number}
        return new{T}(args)
    end
end

Base.length(a::QAdd) = length(a.arguments)

# Equality and hashing
function Base.isequal(a::QAdd, b::QAdd)
    length(a.arguments) == length(b.arguments) || return false
    for (x, y) in zip(a.arguments, b.arguments)
        isequal(x, y) || return false
    end
    return true
end
Base.hash(q::QAdd, h::UInt) = hash(:QAdd, hash(q.arguments, h))

# Adjoint
Base.adjoint(q::QAdd) = QAdd(QMul[adjoint(t) for t in q.arguments])

# Promote
Base.promote_rule(::Type{QAdd{S}}, ::Type{QAdd{T}}) where {S,T} = QAdd{promote_type(S, T)}
function Base.convert(::Type{QAdd{T}}, x::QAdd{S}) where {T<:Number,S<:Number}
    return QAdd(QMul{T}[convert(QMul{T}, m) for m in x.arguments])
end

# Helper: wrap QSym as QMul{T}
_to_qmul(a::QSym, ::Type{T}) where {T} = QMul(one(T), QSym[a])
_to_qmul(a::QMul{T}) where {T} = a
_to_qmul(a::QMul{S}, ::Type{T}) where {S,T} = convert(QMul{T}, a)

# Helper: wrap scalar as QMul
_scalar_qmul(x::T) where {T<:Number} = QMul(x, QSym[])

## Addition — always returns QAdd

# QSym + QSym → QAdd{Int}
function Base.:+(a::QSym, b::QSym)
    return QAdd(QMul{Int}[_to_qmul(a, Int), _to_qmul(b, Int)])
end

# QMul + QMul → QAdd
function Base.:+(a::QMul{S}, b::QMul{T}) where {S,T}
    TT = promote_type(S, T)
    return QAdd(QMul{TT}[convert(QMul{TT}, a), convert(QMul{TT}, b)])
end

# QMul + QSym → QAdd
function Base.:+(a::QMul{T}, b::QSym) where {T}
    return QAdd(QMul{T}[a, _to_qmul(b, T)])
end
Base.:+(a::QSym, b::QMul{T}) where {T} = b + a

# QAdd + QMul → QAdd
function Base.:+(a::QAdd{S}, b::QMul{T}) where {S,T}
    TT = promote_type(S, T)
    args = QMul{TT}[convert(QMul{TT}, x) for x in a.arguments]
    push!(args, convert(QMul{TT}, b))
    return QAdd(args)
end
Base.:+(a::QMul, b::QAdd) = b + a

# QAdd + QSym → QAdd
function Base.:+(a::QAdd{T}, b::QSym) where {T}
    return a + _to_qmul(b, T)
end
Base.:+(a::QSym, b::QAdd) = b + a

# QAdd + QAdd → QAdd
function Base.:+(a::QAdd{S}, b::QAdd{T}) where {S,T}
    TT = promote_type(S, T)
    args = QMul{TT}[convert(QMul{TT}, x) for x in a.arguments]
    for x in b.arguments
        push!(args, convert(QMul{TT}, x))
    end
    return QAdd(args)
end

# QField + Number → QAdd
function Base.:+(a::QSym, b::T) where {T<:Number}
    return QAdd(QMul{T}[_to_qmul(a, T), _scalar_qmul(b)])
end
Base.:+(a::Number, b::QSym) = b + a

function Base.:+(a::QMul{S}, b::T) where {S,T<:Number}
    TT = promote_type(S, T)
    return QAdd(QMul{TT}[convert(QMul{TT}, a), QMul(convert(TT, b), QSym[])])
end
Base.:+(a::Number, b::QMul) = b + a

function Base.:+(a::QAdd{S}, b::T) where {S,T<:Number}
    TT = promote_type(S, T)
    args = QMul{TT}[convert(QMul{TT}, x) for x in a.arguments]
    push!(args, QMul(convert(TT, b), QSym[]))
    return QAdd(args)
end
Base.:+(a::Number, b::QAdd) = b + a

# Subtraction
Base.:-(a::QAdd) = QAdd(QMul[-t for t in a.arguments])
Base.:-(a::QField, b::QField) = a + (-b)
Base.:-(a::QField, b::Number) = a + (-b)
Base.:-(a::Number, b::QField) = a + (-b)

## QAdd * ... (distributive)

# QAdd * Number → QAdd
function Base.:*(a::QAdd{S}, b::T) where {S,T<:Number}
    TT = promote_type(S, T)
    return QAdd(QMul{TT}[convert(QMul{TT}, t) * b for t in a.arguments])
end
Base.:*(a::Number, b::QAdd) = b * a

# QAdd * QSym → QAdd
function Base.:*(a::QAdd{T}, b::QSym) where {T}
    return QAdd(QMul{T}[t * b for t in a.arguments])
end
Base.:*(a::QSym, b::QAdd) = QAdd(QMul[a * t for t in b.arguments])

# QAdd * QMul → QAdd
function Base.:*(a::QAdd{S}, b::QMul{T}) where {S,T}
    TT = promote_type(S, T)
    return QAdd(QMul{TT}[t * b for t in a.arguments])
end
function Base.:*(a::QMul{S}, b::QAdd{T}) where {S,T}
    TT = promote_type(S, T)
    return QAdd(QMul{TT}[a * t for t in b.arguments])
end

# QAdd * QAdd → QAdd
function Base.:*(a::QAdd{S}, b::QAdd{T}) where {S,T}
    TT = promote_type(S, T)
    args = QMul{TT}[]
    for ai in a.arguments, bi in b.arguments
        push!(args, ai * bi)
    end
    return QAdd(args)
end
```

Update `src/SecondQuantizedAlgebra.jl` — add `include("qadd.jl")` after `qmul.jl`, and `export QAdd`.

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project -e 'using Pkg; Pkg.test(test_args=["qadd_test"])'`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/qadd.jl src/SecondQuantizedAlgebra.jl test/qadd_test.jl
git commit -m "feat: add QAdd with type-stable addition and distribution"
```

---

### Task 6: @qnumbers Macro

**Files:**
- Create: `src/macros.jl`
- Modify: `src/SecondQuantizedAlgebra.jl`
- Test: `test/macros_test.jl`
- Reference: `src/legacy/qnumber.jl` (the `@qnumbers` macro section)

- [ ] **Step 1: Write the test file**

```julia
# test/macros_test.jl
using SecondQuantizedAlgebra
using Test

@testset "@qnumbers" begin
    @testset "Single space" begin
        h = FockSpace(:c)
        @qnumbers a::Destroy(h)
        @test a isa Destroy
        @test a.name == :a
        @test a.space_index == 1
    end

    @testset "Multiple operators" begin
        h = FockSpace(:c)
        @qnumbers a::Destroy(h) b::Destroy(h)
        @test a.name == :a
        @test b.name == :b
    end

    @testset "Product space" begin
        h = FockSpace(:a) ⊗ FockSpace(:b)
        @qnumbers a::Destroy(h, 1) b::Destroy(h, 2)
        @test a.space_index == 1
        @test b.space_index == 2
    end

    @testset "Create operators" begin
        h = FockSpace(:c)
        @qnumbers ad::Create(h)
        @test ad isa Create
        @test ad.name == :ad
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project -e 'using Pkg; Pkg.test(test_args=["macros_test"])'`
Expected: FAIL — `@qnumbers` not defined

- [ ] **Step 3: Write the implementation**

```julia
# src/macros.jl

"""
    @qnumbers

Convenience macro for the construction of operators.

Examples
========
```julia
h = FockSpace(:fock)
@qnumbers a::Destroy(h)

h = FockSpace(:one) ⊗ FockSpace(:two)
@qnumbers a::Destroy(h, 1) b::Destroy(h, 2)
```
"""
macro qnumbers(qs...)
    ex = Expr(:block)
    qnames = []
    for q in qs
        @assert q isa Expr && q.head == :(::)
        name = q.args[1]
        @assert name isa Symbol
        push!(qnames, name)
        f = q.args[2]
        @assert f isa Expr && f.head == :call
        op_type = f.args[1]
        op_args = f.args[2:end]
        name_quoted = Expr(:quote, name)
        construction = Expr(:call, esc(op_type), map(esc, op_args)..., name_quoted)
        push!(ex.args, :($(esc(name)) = $(construction)))
    end
    push!(ex.args, Expr(:tuple, map(esc, qnames)...))
    return ex
end
```

Update `src/SecondQuantizedAlgebra.jl` — add `include("macros.jl")` and `export @qnumbers`.

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project -e 'using Pkg; Pkg.test(test_args=["macros_test"])'`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/macros.jl src/SecondQuantizedAlgebra.jl test/macros_test.jl
git commit -m "feat: add @qnumbers macro for operator construction"
```

---

### Task 7: TermInterface Integration

**Files:**
- Create: `src/interface.jl`
- Modify: `src/SecondQuantizedAlgebra.jl`
- Test: `test/interface_test.jl`
- Reference: `/tmp/KeldyshContraction.jl/src/keldysh_algebra/interface.jl`

- [ ] **Step 1: Write the test file**

```julia
# test/interface_test.jl
using SecondQuantizedAlgebra
using SymbolicUtils
using TermInterface
using Test

@testset "TermInterface" begin
    h = FockSpace(:c)
    a = Destroy(h, :a)
    ad = a'

    @testset "QSym is not callable" begin
        @test SymbolicUtils.iscall(a) == false
        @test SymbolicUtils.iscall(ad) == false
    end

    @testset "QMul TermInterface" begin
        m = 3 * a * ad
        @test SymbolicUtils.iscall(m) == true
        @test SymbolicUtils.operation(m) == (*)
        args = SymbolicUtils.arguments(m)
        @test args[1] == 3
        @test length(args) == 3
        @test TermInterface.metadata(m) === nothing
    end

    @testset "QAdd TermInterface" begin
        s = a + ad
        @test SymbolicUtils.iscall(s) == true
        @test SymbolicUtils.operation(s) == (+)
        args = SymbolicUtils.arguments(s)
        @test length(args) == 2
        @test all(x -> x isa QMul, args)
        @test TermInterface.metadata(s) === nothing
    end

    @testset "maketerm QMul" begin
        m = 2 * a * ad
        args = SymbolicUtils.arguments(m)
        m2 = TermInterface.maketerm(typeof(m), *, args, nothing)
        @test isequal(m, m2)
    end

    @testset "maketerm QAdd" begin
        s = a + ad
        args = SymbolicUtils.arguments(s)
        s2 = TermInterface.maketerm(typeof(s), +, args, nothing)
        @test isequal(s, s2)
    end

    @testset "symtype" begin
        @test SymbolicUtils.symtype(a) == Destroy
        @test SymbolicUtils.symtype(ad) == Create
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project -e 'using Pkg; Pkg.test(test_args=["interface_test"])'`
Expected: FAIL — SymbolicUtils methods not defined

- [ ] **Step 3: Write the implementation**

```julia
# src/interface.jl

## TermInterface / SymbolicUtils integration

# Head
TermInterface.head(::QField) = :call

# QSym — leaves (not callable)
SymbolicUtils.iscall(::QSym) = false
TermInterface.metadata(::QSym) = nothing

# QMul — products
SymbolicUtils.iscall(::QMul) = true
SymbolicUtils.iscall(::Type{<:QMul}) = true
SymbolicUtils.operation(::QMul) = (*)
SymbolicUtils.arguments(a::QMul) = Any[a.arg_c, a.args_nc...]
TermInterface.metadata(::QMul) = nothing

function TermInterface.maketerm(::Type{<:QMul}, ::typeof(*), args, metadata)
    args_c = filter(x -> !(x isa QField), args)
    args_nc = filter(x -> x isa QField, args)
    arg_c = isempty(args_c) ? 1 : *(args_c...)
    return QMul(arg_c, QSym[args_nc...])
end

# QAdd — sums
SymbolicUtils.iscall(::QAdd) = true
SymbolicUtils.iscall(::Type{<:QAdd}) = true
SymbolicUtils.operation(::QAdd) = (+)
SymbolicUtils.arguments(a::QAdd) = a.arguments
TermInterface.metadata(::QAdd) = nothing
TermInterface.maketerm(::Type{<:QAdd}, ::typeof(+), args, metadata) = QAdd(args)

# Type promotion
SymbolicUtils.symtype(x::T) where {T<:QField} = T

for f in SymbolicUtils.basic_diadic
    @eval SymbolicUtils.promote_symtype(::$(typeof(f)), Ts::Type{<:QField}...) = promote_type(Ts...)
    @eval SymbolicUtils.promote_symtype(::$(typeof(f)), T::Type{<:QField}, Ts...) = T
    @eval SymbolicUtils.promote_symtype(::$(typeof(f)), T::Type{<:QField}, S::Type{<:Number}) = T
    @eval SymbolicUtils.promote_symtype(::$(typeof(f)), T::Type{<:Number}, S::Type{<:QField}) = S
    @eval SymbolicUtils.promote_symtype(::$(typeof(f)), T::Type{<:QField}, S::Type{<:QField}) = promote_type(T, S)
end

# one / zero / isone / iszero for QField
Base.one(::T) where {T<:QField} = one(T)
Base.one(::Type{<:QField}) = 1
Base.isone(::QField) = false
Base.zero(::T) where {T<:QField} = zero(T)
Base.zero(::Type{<:QField}) = 0
```

Update `src/SecondQuantizedAlgebra.jl` — add `include("interface.jl")` after `qadd.jl`.

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project -e 'using Pkg; Pkg.test(test_args=["interface_test"])'`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/interface.jl src/SecondQuantizedAlgebra.jl test/interface_test.jl
git commit -m "feat: add TermInterface/SymbolicUtils integration for QMul, QAdd"
```

---

### Task 8: Normal Ordering

**Files:**
- Create: `src/normal_order.jl`
- Modify: `src/SecondQuantizedAlgebra.jl`
- Test: `test/normal_order_test.jl`

- [ ] **Step 1: Write the test file**

```julia
# test/normal_order_test.jl
using SecondQuantizedAlgebra
using Test

@testset "normal_order" begin
    h = FockSpace(:c)
    a = Destroy(h, :a)
    ad = a'

    @testset "Already normal-ordered" begin
        # a†a is already normal-ordered
        m = ad * a
        result = normal_order(m)
        @test result isa QAdd
        @test length(result.arguments) == 1
        @test isequal(result.arguments[1], m)
    end

    @testset "Single commutation: a*a† → a†a + 1" begin
        m = a * ad  # stored as QMul(1, [Create, Destroy]) due to canonical sort
        result = normal_order(m)
        @test result isa QAdd
        # Should have two terms: a†a and 1
        @test length(result.arguments) == 2
    end

    @testset "QSym passthrough" begin
        result = normal_order(a)
        @test result isa QAdd
        @test length(result.arguments) == 1
    end

    @testset "QAdd — normal-orders each term" begin
        expr = a * ad + ad * a
        result = normal_order(expr)
        @test result isa QAdd
    end

    @testset "Multi-operator: a*a†*a† → a†*a†*a + 2*a†" begin
        m = a * ad * ad
        result = normal_order(m)
        @test result isa QAdd
        # a a† a† = a† a a† + a† = a† (a† a + 1) + a† = a†a†a + a† + a† = a†a†a + 2a†
        @test length(result.arguments) >= 2
    end

    @testset "Different spaces don't commute" begin
        h2 = FockSpace(:a) ⊗ FockSpace(:b)
        a1 = Destroy(h2, :a, 1)
        a2 = Destroy(h2, :b, 2)
        # a1 * a2† — different spaces, already normal-ordered within each space
        m = a1 * a2'
        result = normal_order(m)
        @test result isa QAdd
        @test length(result.arguments) == 1  # no commutation needed
    end

    @testset "Return type is always QAdd" begin
        @test normal_order(a) isa QAdd
        @test normal_order(2 * a) isa QAdd
        @test normal_order(a + ad) isa QAdd
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project -e 'using Pkg; Pkg.test(test_args=["normal_order_test"])'`
Expected: FAIL — `normal_order` not defined

- [ ] **Step 3: Write the implementation**

```julia
# src/normal_order.jl

"""
    normal_order(expr)

Apply bosonic commutation relation [a, a†] = 1 to rewrite the expression
in normal order (all creation operators to the left of annihilation operators).
Always returns `QAdd`.

Does NOT collect like terms — use `simplify()` for that.
"""
function normal_order(op::QSym)
    return QAdd(QMul{Int}[QMul(1, QSym[op])])
end

function normal_order(m::QMul{T}) where {T}
    return _normal_order_qmul(m.arg_c, m.args_nc)
end

function normal_order(s::QAdd{T}) where {T}
    result = QMul{T}[]
    for term in s.arguments
        ordered = normal_order(term)
        append!(result, ordered.arguments)
    end
    return QAdd(result)
end

"""
    _normal_order_qmul(prefactor, ops)

Recursively normal-order a product of operators. Uses bubble-sort approach:
scan for adjacent Destroy-Create pairs on the same space, swap and add identity.
"""
function _normal_order_qmul(arg_c::T, ops::Vector{QSym}) where {T}
    isempty(ops) && return QAdd(QMul{T}[QMul(arg_c, QSym[])])
    length(ops) == 1 && return QAdd(QMul{T}[QMul(arg_c, copy(ops))])

    # Find first out-of-normal-order pair on the same space
    for i in 1:(length(ops) - 1)
        a, b = ops[i], ops[i + 1]
        if a isa Destroy && b isa Create && a.space_index == b.space_index
            # Apply [a, a†] = 1: swap a,a† → a†,a and add identity term

            # Term 1: swapped (a† a) — continue normal ordering
            swapped = QSym[ops[1:i-1]..., b, a, ops[i+2:end]...]
            sort!(swapped; lt=canonical_lt)
            term1 = _normal_order_qmul(arg_c, swapped)

            # Term 2: identity (contraction) — remove both, multiply by δ
            contracted = QSym[ops[1:i-1]..., ops[i+2:end]...]
            sort!(contracted; lt=canonical_lt)
            term2 = _normal_order_qmul(arg_c, contracted)

            # Combine
            all_terms = QMul{T}[term1.arguments..., term2.arguments...]
            return QAdd(all_terms)
        end
    end

    # Already in normal order
    return QAdd(QMul{T}[QMul(arg_c, ops)])
end
```

Update `src/SecondQuantizedAlgebra.jl` — add `include("normal_order.jl")` and `export normal_order`.

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project -e 'using Pkg; Pkg.test(test_args=["normal_order_test"])'`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/normal_order.jl src/SecondQuantizedAlgebra.jl test/normal_order_test.jl
git commit -m "feat: add normal_order() with bosonic commutation relations"
```

---

### Task 9: Simplify

**Files:**
- Create: `src/simplify.jl`
- Modify: `src/SecondQuantizedAlgebra.jl`
- Test: `test/simplify_test.jl`

- [ ] **Step 1: Write the test file**

```julia
# test/simplify_test.jl
using SecondQuantizedAlgebra
using SymbolicUtils
using Symbolics
using Test

@testset "simplify" begin
    h = FockSpace(:c)
    a = Destroy(h, :a)
    ad = a'

    @testset "Collect like terms" begin
        s = QAdd(QMul{Int}[QMul(2, QSym[ad, a]), QMul(3, QSym[ad, a])])
        result = simplify(s)
        @test result isa QAdd
        @test length(result.arguments) == 1
        @test result.arguments[1].arg_c == 5
    end

    @testset "Remove zero terms" begin
        s = QAdd(QMul{Int}[QMul(0, QSym[ad, a]), QMul(3, QSym[ad, a])])
        result = simplify(s)
        @test length(result.arguments) == 1
        @test result.arguments[1].arg_c == 3
    end

    @testset "a + a = 2a" begin
        s = a + a
        result = simplify(s)
        @test length(result.arguments) == 1
        @test result.arguments[1].arg_c == 2
    end

    @testset "Symbolic prefactors" begin
        @variables g h_sym
        s = QAdd([QMul(g, QSym[ad, a]), QMul(h_sym, QSym[ad, a])])
        result = simplify(s)
        @test length(result.arguments) == 1
        # Prefactor should be g + h_sym (a symbolic expression)
    end

    @testset "SymbolicUtils.simplify on QField" begin
        @variables g
        s = QAdd([QMul(g, QSym[ad, a]), QMul(g, QSym[ad, a])])
        result = SymbolicUtils.simplify(s)
        @test result isa QAdd
    end

    @testset "Symbolics.expand on QField" begin
        expr = (a + ad) * (a + ad)
        result = Symbolics.expand(expr)
        @test result isa QAdd
        @test length(result.arguments) == 4
    end

    @testset "Return type is QAdd" begin
        @test simplify(a + ad) isa QAdd
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project -e 'using Pkg; Pkg.test(test_args=["simplify_test"])'`
Expected: FAIL — `simplify` not defined on new types

- [ ] **Step 3: Write the implementation**

```julia
# src/simplify.jl

"""
    simplify(s::QAdd)

Group terms with identical `args_nc`, sum their `arg_c` prefactors.
Remove zero-prefactor terms. Returns `QAdd`.
"""
function simplify(s::QAdd{T}) where {T}
    # Group by args_nc hash
    groups = Dict{UInt, Tuple{Vector{QSym}, Any}}()
    order = UInt[]  # preserve insertion order

    for term in s.arguments
        key = hash(term.args_nc)
        if haskey(groups, key)
            existing_ops, existing_c = groups[key]
            # Verify it's actually the same args_nc (not just hash collision)
            if isequal(existing_ops, term.args_nc)
                groups[key] = (existing_ops, existing_c + term.arg_c)
            else
                # Hash collision — use a different key
                key2 = hash(term.args_nc, key)
                if haskey(groups, key2)
                    _, ec = groups[key2]
                    groups[key2] = (term.args_nc, ec + term.arg_c)
                else
                    groups[key2] = (term.args_nc, term.arg_c)
                    push!(order, key2)
                end
            end
        else
            groups[key] = (term.args_nc, term.arg_c)
            push!(order, key)
        end
    end

    # Build result, skipping zeros
    result = QMul[]
    for key in order
        ops, c = groups[key]
        iszero(c) && continue
        push!(result, QMul(c, ops))
    end

    isempty(result) && return QAdd(QMul[QMul(zero(T), QSym[])])
    return QAdd(result)
end

simplify(m::QMul) = QAdd([m])
simplify(op::QSym) = QAdd(QMul{Int}[QMul(1, QSym[op])])

# SymbolicUtils.simplify — also simplify each prefactor
function SymbolicUtils.simplify(s::QAdd; kwargs...)
    simplified = simplify(s)
    result = QMul[
        QMul(_simplify_prefactor(t.arg_c; kwargs...), t.args_nc)
        for t in simplified.arguments
    ]
    return QAdd(result)
end
SymbolicUtils.simplify(m::QMul; kwargs...) = SymbolicUtils.simplify(QAdd([m]); kwargs...)
SymbolicUtils.simplify(op::QSym; kwargs...) = SymbolicUtils.simplify(simplify(op); kwargs...)

_simplify_prefactor(x::Number; kwargs...) = x
function _simplify_prefactor(x; kwargs...)
    return SymbolicUtils.simplify(x; kwargs...)
end

# Symbolics.expand — distribute products, then expand prefactors
function Symbolics.expand(s::QAdd; kwargs...)
    # QAdd is already expanded (sum of products)
    result = QMul[
        QMul(_expand_prefactor(t.arg_c; kwargs...), t.args_nc)
        for t in s.arguments
    ]
    return QAdd(result)
end
Symbolics.expand(m::QMul; kwargs...) = Symbolics.expand(QAdd([m]); kwargs...)
Symbolics.expand(op::QSym; kwargs...) = QAdd(QMul{Int}[QMul(1, QSym[op])])

_expand_prefactor(x::Number; kwargs...) = x
function _expand_prefactor(x; kwargs...)
    return Symbolics.expand(x; kwargs...)
end
```

Update `src/SecondQuantizedAlgebra.jl` — add `include("simplify.jl")` and `export simplify`.

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project -e 'using Pkg; Pkg.test(test_args=["simplify_test"])'`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/simplify.jl src/SecondQuantizedAlgebra.jl test/simplify_test.jl
git commit -m "feat: add simplify(), SymbolicUtils.simplify, Symbolics.expand for QField"
```

---

### Task 10: Printing (REPL Display)

**Files:**
- Create: `src/printing.jl`
- Modify: `src/SecondQuantizedAlgebra.jl`
- Test: `test/printing_test.jl`
- Reference: `src/legacy/printing.jl`

- [ ] **Step 1: Write the test file**

```julia
# test/printing_test.jl
using SecondQuantizedAlgebra
using Test

@testset "printing" begin
    h = FockSpace(:cavity)
    a = Destroy(h, :a)
    ad = a'

    @testset "HilbertSpace" begin
        @test repr(h) == "ℋ(cavity)"
        h2 = FockSpace(:a) ⊗ FockSpace(:b)
        @test repr(h2) == "ℋ(a) ⊗ ℋ(b)"
    end

    @testset "Operators" begin
        @test repr(a) == "a"
        @test repr(ad) == "a†"
    end

    @testset "QMul" begin
        @test repr(3 * ad * a) == "3 * a† * a"
        @test repr(1 * ad * a) == "a† * a"
        @test repr(-1 * a) == "-a"
        @test repr(QMul(5, QSym[])) == "5"
    end

    @testset "QAdd" begin
        @test repr(ad * a + 1) == "a† * a + 1"
        s = ad * a + (-1 * a)
        r = repr(s)
        @test occursin("a† * a", r)
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project -e 'using Pkg; Pkg.test(test_args=["printing_test"])'`
Expected: FAIL — show methods not defined

- [ ] **Step 3: Write the implementation**

```julia
# src/printing.jl

# HilbertSpace
Base.show(io::IO, h::FockSpace) = write(io, "ℋ(", string(h.name), ")")
function Base.show(io::IO, h::ProductSpace)
    show(io, h.spaces[1])
    for i in 2:length(h.spaces)
        write(io, " ⊗ ")
        show(io, h.spaces[i])
    end
end

# Operators
Base.show(io::IO, x::Destroy) = write(io, string(x.name))
Base.show(io::IO, x::Create) = write(io, string(x.name), "†")

# QMul
function Base.show(io::IO, x::QMul)
    if isempty(x.args_nc)
        # Pure scalar
        show(io, x.arg_c)
        return
    end
    if x.arg_c == -1
        write(io, "-")
    elseif !isone(x.arg_c)
        show(io, x.arg_c)
        write(io, " * ")
    end
    show(io, x.args_nc[1])
    for i in 2:length(x.args_nc)
        write(io, " * ")
        show(io, x.args_nc[i])
    end
end

# QAdd
function Base.show(io::IO, x::QAdd)
    isempty(x.arguments) && return write(io, "0")
    show(io, x.arguments[1])
    for i in 2:length(x.arguments)
        term = x.arguments[i]
        if term.arg_c isa Real && term.arg_c < 0
            write(io, " - ")
            show(io, QMul(-term.arg_c, term.args_nc))
        else
            write(io, " + ")
            show(io, term)
        end
    end
end
```

Update `src/SecondQuantizedAlgebra.jl` — add `include("printing.jl")`.

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project -e 'using Pkg; Pkg.test(test_args=["printing_test"])'`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/printing.jl src/SecondQuantizedAlgebra.jl test/printing_test.jl
git commit -m "feat: add REPL display for all types"
```

---

### Task 11: LaTeX Recipes

**Files:**
- Create: `src/latexify_recipes.jl`
- Modify: `src/SecondQuantizedAlgebra.jl`
- Test: `test/latexify_test.jl`
- Reference: `src/legacy/latexify_recipes.jl`

- [ ] **Step 1: Write the test file**

```julia
# test/latexify_test.jl
using SecondQuantizedAlgebra
using Latexify
using Test

@testset "latexify" begin
    h = FockSpace(:c)
    a = Destroy(h, :a)
    ad = a'

    @testset "Operators" begin
        @test latexify(a) == L"a"
        @test latexify(ad) == L"a^{\dagger}"
    end

    @testset "QMul" begin
        s = latexify(3 * ad * a)
        @test occursin("3", s)
        @test occursin("a^{\\dagger}", s)
        @test occursin("a", s)

        # Unit prefactor omitted
        s2 = latexify(1 * ad * a)
        @test !startswith(s2, "\$1")

        # Scalar
        @test latexify(QMul(5, QSym[])) == L"5"
    end

    @testset "QAdd" begin
        s = latexify(ad * a + 1)
        @test occursin("a^{\\dagger}", s)
        @test occursin("+", s)
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project -e 'using Pkg; Pkg.test(test_args=["latexify_test"])'`
Expected: FAIL — no latexrecipe defined

- [ ] **Step 3: Write the implementation**

```julia
# src/latexify_recipes.jl

@latexrecipe function f(x::Destroy)
    return Expr(:latexifymerge, x.name)
end

@latexrecipe function f(x::Create)
    return Expr(:latexifymerge, "$(x.name)^{\\dagger}")
end

@latexrecipe function f(x::QMul)
    if isempty(x.args_nc)
        return x.arg_c
    end
    parts = []
    if x.arg_c == -1
        push!(parts, :(-))
    elseif !isone(x.arg_c)
        push!(parts, x.arg_c)
    end
    for op in x.args_nc
        push!(parts, op)
    end
    return Expr(:latexifymerge, parts...)
end

@latexrecipe function f(x::QAdd)
    return Expr(:call, :+, x.arguments...)
end

# LaTeX MIME show
const QLaTeX = Union{<:QField}
Base.show(io::IO, ::MIME"text/latex", x::QLaTeX) = write(io, latexify(x))
```

Update `src/SecondQuantizedAlgebra.jl` — add `include("latexify_recipes.jl")`.

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project -e 'using Pkg; Pkg.test(test_args=["latexify_test"])'`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/latexify_recipes.jl src/SecondQuantizedAlgebra.jl test/latexify_test.jl
git commit -m "feat: add LaTeX rendering via Latexify.jl"
```

---

### Task 12: Numerical Conversion (to_numeric, numeric_average)

**Files:**
- Create: `src/numeric.jl`
- Modify: `src/SecondQuantizedAlgebra.jl`
- Test: `test/numeric_test.jl`
- Reference: `src/legacy/utils.jl` (lines 134–442)

- [ ] **Step 1: Write the test file**

```julia
# test/numeric_test.jl
using SecondQuantizedAlgebra
using QuantumOpticsBase
using Test

@testset "numeric conversion" begin
    @testset "Single space — basic" begin
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)

        @test to_numeric(a, b) == destroy(b)
        @test to_numeric(a', b) == create(b)
    end

    @testset "QMul" begin
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)

        @test to_numeric(a' * a, b) == create(b) * destroy(b)
        @test to_numeric(2 * a, b) == 2 * destroy(b)
    end

    @testset "QAdd" begin
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)

        result = to_numeric(a + a', b)
        expected = destroy(b) + create(b)
        @test result == expected
    end

    @testset "Scalar" begin
        h = FockSpace(:fock)
        b = FockBasis(7)

        @test to_numeric(3, b) == 3 * one(b)
    end

    @testset "Product space" begin
        h = FockSpace(:a) ⊗ FockSpace(:b)
        @qnumbers a1::Destroy(h, 1) a2::Destroy(h, 2)
        b = FockBasis(3) ⊗ FockBasis(3)

        a1_num = to_numeric(a1, b)
        @test a1_num isa LazyTensor

        prod_num = to_numeric(a1' * a2, b)
        @test prod_num isa LazyTensor
    end

    @testset "numeric_average" begin
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)
        α = 0.1 + 0.2im
        ψ = coherentstate(b, α)

        @test numeric_average(a, ψ) ≈ α
        @test numeric_average(a' * a, ψ) ≈ abs2(α)
    end

    @testset "to_numeric and numeric_average export" begin
        @test isdefined(SecondQuantizedAlgebra, :to_numeric)
        @test isdefined(SecondQuantizedAlgebra, :numeric_average)
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project -e 'using Pkg; Pkg.test(test_args=["numeric_test"])'`
Expected: FAIL — `to_numeric` not defined

- [ ] **Step 3: Write the implementation**

```julia
# src/numeric.jl

"""
    to_numeric(op, basis; kwargs...)

Convert a symbolic operator to its numerical matrix form on the given basis.
"""
function to_numeric(op::Destroy, b::QuantumOpticsBase.FockBasis; kwargs...)
    return QuantumOpticsBase.destroy(b)
end
function to_numeric(op::Create, b::QuantumOpticsBase.FockBasis; kwargs...)
    return QuantumOpticsBase.create(b)
end

# Composite basis — embed single operator
function to_numeric(op::QSym, b::QuantumOpticsBase.CompositeBasis; kwargs...)
    idx = op.space_index
    op_num = to_numeric(op, b.bases[idx]; kwargs...)
    return QuantumOpticsBase.LazyTensor(b, [idx], (op_num,))
end

# QMul
function to_numeric(m::QMul, b::QuantumOpticsBase.Basis; kwargs...)
    if isempty(m.args_nc)
        return m.arg_c * _lazy_one(b)
    end
    ops_num = [to_numeric(op, b; kwargs...) for op in m.args_nc]
    return m.arg_c * prod(ops_num)
end

# QAdd
function to_numeric(s::QAdd, b::QuantumOpticsBase.Basis; kwargs...)
    terms_num = [to_numeric(t, b; kwargs...) for t in s.arguments]
    return sum(terms_num)
end

# Number
function to_numeric(x::Number, b::QuantumOpticsBase.Basis; kwargs...)
    return x * _lazy_one(b)
end

# State-based dispatch
function to_numeric(op, state; kwargs...)
    return to_numeric(op, QuantumOpticsBase.basis(state); kwargs...)
end

# Lazy identity helper
_lazy_one(b::QuantumOpticsBase.Basis) = one(b)
function _lazy_one(b::QuantumOpticsBase.CompositeBasis)
    return QuantumOpticsBase.LazyTensor(
        b, collect(1:length(b.bases)), Tuple(one(bi) for bi in b.bases)
    )
end

"""
    numeric_average(op, state; kwargs...)

Compute the expectation value of a symbolic operator with a quantum state.
"""
function numeric_average(op::QField, state; kwargs...)
    op_num = to_numeric(op, state; kwargs...)
    return QuantumOpticsBase.expect(op_num, state)
end
function numeric_average(x::Number, state; kwargs...)
    return x
end
```

Update `src/SecondQuantizedAlgebra.jl` — add `include("numeric.jl")` and `export to_numeric, numeric_average`.

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project -e 'using Pkg; Pkg.test(test_args=["numeric_test"])'`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/numeric.jl src/SecondQuantizedAlgebra.jl test/numeric_test.jl
git commit -m "feat: add to_numeric and numeric_average via QuantumOpticsBase"
```

---

### Task 13: Code Quality Tests

**Files:**
- Create: `test/code_quality_test.jl`

- [ ] **Step 1: Write the test file**

```julia
# test/code_quality_test.jl
using SecondQuantizedAlgebra
using Test

@testset "best practices" begin
    using Aqua
    Aqua.test_all(SecondQuantizedAlgebra; ambiguities=false, piracies=false)
end

@testset "ExplicitImports" begin
    using ExplicitImports
    @test check_no_implicit_imports(SecondQuantizedAlgebra) === nothing
    @test check_all_explicit_imports_via_owners(SecondQuantizedAlgebra) === nothing
    @test check_no_stale_explicit_imports(SecondQuantizedAlgebra) === nothing
    @test check_all_qualified_accesses_via_owners(SecondQuantizedAlgebra) === nothing
    @test check_no_self_qualified_accesses(SecondQuantizedAlgebra) === nothing
end

if isempty(VERSION.prerelease)
    @testset "Code linting" begin
        using JET
        rep = report_package(SecondQuantizedAlgebra)
        @show rep
        @test length(JET.get_reports(rep)) == 0
    end
end

@testset "Concretely typed" begin
    using CheckConcreteStructs
    all_concrete(SecondQuantizedAlgebra.FockSpace)
    all_concrete(SecondQuantizedAlgebra.ProductSpace)
    all_concrete(SecondQuantizedAlgebra.Destroy)
    all_concrete(SecondQuantizedAlgebra.Create)
    all_concrete(SecondQuantizedAlgebra.QMul)
    all_concrete(SecondQuantizedAlgebra.QAdd)
end
```

- [ ] **Step 2: Run test**

Run: `julia --project -e 'using Pkg; Pkg.test(test_args=["code_quality_test"])'`
Expected: PASS — zero JET reports, all structs concrete, no import issues.

- [ ] **Step 3: Fix any issues found**

Address any Aqua, JET, ExplicitImports, or CheckConcreteStructs failures.

- [ ] **Step 4: Commit**

```bash
git add test/code_quality_test.jl
git commit -m "test: add code quality tests (Aqua, JET, ExplicitImports, CheckConcreteStructs)"
```

---

### Task 14: Integration Test

**Files:**
- Create: `test/integration_test.jl`

This test exercises the full workflow from the old `test/legacy/fock_test.jl`, adapted for the new API (lazy `*`, explicit `normal_order` + `simplify`).

- [ ] **Step 1: Write the test file**

```julia
# test/integration_test.jl
using SecondQuantizedAlgebra
using SymbolicUtils
using Symbolics
using QuantumOpticsBase
using Test

@testset "integration" begin
    @testset "Basic operator creation and algebra" begin
        h = FockSpace(:c)
        @qnumbers a::Destroy(h)
        ad = a'

        # Hash uniqueness
        @test hash(a) != hash(ad)
        @test isequal(a, ad')

        # Lazy multiplication
        @test a * ad isa QMul
        @test ad * a isa QMul
    end

    @testset "Commutation via normal_order" begin
        h = FockSpace(:c)
        @qnumbers a::Destroy(h)

        # [a, a†] = 1 via normal ordering
        # a * a† normal-ordered → a†a + 1
        m = a * a'
        result = simplify(normal_order(m))
        # Should have a†a term and scalar 1 term
        @test result isa QAdd
        @test length(result.arguments) == 2
    end

    @testset "Hamiltonian evolution" begin
        h = FockSpace(:c)
        @qnumbers a::Destroy(h)

        ωc = 0.1313
        H = ωc * a' * a

        # [H, a] = ωc * (a'a*a - a*a'a) normal-ordered = -ωc * a
        Ha = H * a
        aH = a * H
        comm = Ha + (-aH)  # QAdd
        result = simplify(normal_order(comm))
        # Should simplify to -ωc * a
        @test result isa QAdd
    end

    @testset "Numeric conversion round-trip" begin
        h = FockSpace(:c)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)

        @test to_numeric(a, b) == destroy(b)
        @test to_numeric(a', b) == create(b)
        @test to_numeric(a' * a, b) == create(b) * destroy(b)

        α = 0.1 + 0.2im
        ψ = coherentstate(b, α)
        @test numeric_average(a, ψ) ≈ α
        @test numeric_average(a' * a, ψ) ≈ abs2(α)
    end

    @testset "Multi-space" begin
        h = FockSpace(:a) ⊗ FockSpace(:b)
        @qnumbers a::Destroy(h, 1) b::Destroy(h, 2)

        # Operators on different spaces — normal_order doesn't commute them
        m = a * b'
        result = normal_order(m)
        @test length(result.arguments) == 1  # no commutation

        # Same-space commutation still works
        m2 = a * a'
        result2 = simplify(normal_order(m2))
        @test length(result2.arguments) == 2  # a'a + 1
    end

    @testset "Symbolic parameters" begin
        @variables g ω
        h = FockSpace(:c)
        @qnumbers a::Destroy(h)

        expr = g * a' * a + ω * a' * a
        result = simplify(expr)
        @test length(result.arguments) == 1
        # Prefactor should be g + ω
    end
end
```

- [ ] **Step 2: Run test**

Run: `julia --project -e 'using Pkg; Pkg.test(test_args=["integration_test"])'`
Expected: PASS

- [ ] **Step 3: Fix any issues found**

Address failures from integration testing.

- [ ] **Step 4: Commit**

```bash
git add test/integration_test.jl
git commit -m "test: add integration tests covering full workflow"
```

---

### Task 15: Update Project.toml and Final Cleanup

**Files:**
- Modify: `Project.toml`
- Modify: `src/SecondQuantizedAlgebra.jl`

- [ ] **Step 1: Update dependency compat bounds**

Update `Project.toml` compat to latest releases. At minimum:
- `SymbolicUtils = "3.26"` (match KeldyshContraction)
- `Symbolics = "6"` (keep current)
- `TermInterface = "2"`
- `Latexify = "0.16"` (match KeldyshContraction)
- Remove `LinearAlgebra` from `[deps]` if unused in new code

Verify with: `julia --project -e 'using Pkg; Pkg.update(); Pkg.resolve()'`

- [ ] **Step 2: Verify final module file**

Ensure `src/SecondQuantizedAlgebra.jl` has all includes in correct order and complete export list:

```julia
module SecondQuantizedAlgebra

using SymbolicUtils: SymbolicUtils, BasicSymbolic, arguments, iscall, operation, substitute
using Symbolics: Symbolics
using TermInterface: TermInterface

using Combinatorics: combinations

using SciMLBase: SciMLBase
using QuantumOpticsBase: QuantumOpticsBase
import QuantumOpticsBase: ⊗, tensor

using LaTeXStrings: LaTeXStrings, @L_str, latexstring
using Latexify: Latexify, latexify, @latexrecipe
using MacroTools: MacroTools

include("types.jl")
include("hilbertspace.jl")
include("fock.jl")
include("qmul.jl")
include("qadd.jl")
include("macros.jl")
include("interface.jl")
include("normal_order.jl")
include("simplify.jl")
include("printing.jl")
include("latexify_recipes.jl")
include("numeric.jl")

export HilbertSpace, FockSpace, ProductSpace,
    ⊗, tensor,
    QField, QSym, QTerm,
    Destroy, Create,
    QMul, QAdd,
    @qnumbers,
    normal_order, simplify,
    to_numeric, numeric_average

end
```

- [ ] **Step 3: Run all tests**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: ALL PASS

- [ ] **Step 4: Commit**

```bash
git add Project.toml src/SecondQuantizedAlgebra.jl
git commit -m "chore: update deps, finalize module structure and exports"
```
