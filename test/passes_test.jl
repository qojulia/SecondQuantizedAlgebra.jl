using Test
using SecondQuantizedAlgebra
using SecondQuantizedAlgebra: _canonicalize_to_dict!, QTermDict, _CNUM_ONE,
    _CNUM_ZERO, _EMPTY_NE, QSym, CNum,
    _reduce_ops, _commute_ops, _expand_gs_ops, _substitute_ops,
    _stream!, _canonicalize!

@testset "canonicalize_to_dict! basic insert" begin
    h = FockSpace(:f)
    a = Destroy(h, :a)
    out = QTermDict()
    _canonicalize_to_dict!(out, QSym[a], _CNUM_ONE, _EMPTY_NE)
    @test length(out) == 1
end

@testset "canonicalize_to_dict! like-term collection" begin
    h = FockSpace(:f)
    a = Destroy(h, :a)
    out = QTermDict()
    _canonicalize_to_dict!(out, QSym[a], _CNUM_ONE, _EMPTY_NE)
    _canonicalize_to_dict!(out, QSym[a], _CNUM_ONE, _EMPTY_NE)
    @test length(out) == 1
    @test first(values(out)) == 2 + 0im
end

@testset "canonicalize_to_dict! zero coeff dropped" begin
    h = FockSpace(:f)
    a = Destroy(h, :a)
    out = QTermDict()
    _canonicalize_to_dict!(out, QSym[a], _CNUM_ZERO, _EMPTY_NE)
    @test isempty(out)
end

@testset "_reduce_ops: Transition composition" begin
    h = NLevelSpace(:a, 3)
    σ12 = Transition(h, :σ, 1, 2)
    σ23 = Transition(h, :σ, 2, 3)
    emitted = Tuple{Vector{QSym}, CNum}[]
    _reduce_ops(QSym[σ12, σ23], _CNUM_ONE) do o, c
        push!(emitted, (copy(o), c))
    end
    @test length(emitted) == 1
    @test length(emitted[1][1]) == 1
    @test emitted[1][1][1].i == 1 && emitted[1][1][1].j == 3
end

@testset "_reduce_ops: zero from incompatible composition" begin
    h = NLevelSpace(:a, 3)
    σ12 = Transition(h, :σ, 1, 2)
    σ31 = Transition(h, :σ, 3, 1)
    emitted = Tuple{Vector{QSym}, CNum}[]
    _reduce_ops(QSym[σ12, σ31], _CNUM_ONE) do o, c
        push!(emitted, (copy(o), c))
    end
    @test isempty(emitted)
end

@testset "_reduce_ops: no-op input passes through" begin
    h = FockSpace(:f)
    a = Destroy(h, :a)
    ad = adjoint(a)
    emitted = Tuple{Vector{QSym}, CNum}[]
    _reduce_ops(QSym[ad, a], _CNUM_ONE) do o, c
        push!(emitted, (copy(o), c))
    end
    @test length(emitted) == 1
    @test emitted[1][1] == QSym[ad, a]
end

@testset "_commute_ops: Fock aa† → a†a + 1" begin
    h = FockSpace(:f)
    a = Destroy(h, :a)
    ad = adjoint(a)
    emitted = Tuple{Vector{QSym}, CNum}[]
    _commute_ops(QSym[a, ad], _CNUM_ONE) do o, c
        push!(emitted, (copy(o), c))
    end
    @test length(emitted) == 2
    sort!(emitted, by = e -> length(e[1]))
    @test isempty(emitted[1][1]) && emitted[1][2] == _CNUM_ONE
    @test emitted[2][1] == QSym[ad, a]
end

@testset "_commute_ops: no-op on already-ordered pair" begin
    h = FockSpace(:f)
    a = Destroy(h, :a)
    ad = adjoint(a)
    emitted = Tuple{Vector{QSym}, CNum}[]
    _commute_ops(QSym[ad, a], _CNUM_ONE) do o, c
        push!(emitted, (copy(o), c))
    end
    @test length(emitted) == 1
    @test emitted[1][1] == QSym[ad, a]
end

@testset "_expand_gs_ops: σ¹¹ → 1 - σ²²" begin
    h = NLevelSpace(:a, 2)
    σ11 = Transition(h, :σ, 1, 1)
    emitted = Tuple{Vector{QSym}, CNum}[]
    _expand_gs_ops(QSym[σ11], _CNUM_ONE) do o, c
        push!(emitted, (copy(o), c))
    end
    @test length(emitted) == 2
    sort!(emitted, by = e -> length(e[1]))
    @test isempty(emitted[1][1]) && emitted[1][2] == _CNUM_ONE
    @test length(emitted[2][1]) == 1 && emitted[2][2] == -_CNUM_ONE
    @test emitted[2][1][1].i == 2 && emitted[2][1][1].j == 2
end

@testset "_expand_gs_ops: passthrough when no σᵍᵍ" begin
    h = FockSpace(:f)
    a = Destroy(h, :a)
    emitted = Tuple{Vector{QSym}, CNum}[]
    _expand_gs_ops(QSym[a], _CNUM_ONE) do o, c
        push!(emitted, (copy(o), c))
    end
    @test length(emitted) == 1
    @test emitted[1][1] == QSym[a]
end

@testset "_substitute_ops: operator → scalar" begin
    h = FockSpace(:f)
    a = Destroy(h, :a)
    emitted = Tuple{Vector{QSym}, CNum}[]
    _substitute_ops(QSym[a, adjoint(a)], _CNUM_ONE, Dict(a => 2)) do o, c
        push!(emitted, (copy(o), c))
    end
    @test length(emitted) == 1
    @test emitted[1][2] == 2 + 0im
    @test emitted[1][1] == QSym[adjoint(a)]
end

@testset "_stream!: idempotent on canonical input" begin
    h = FockSpace(:f)
    a = Destroy(h, :a)
    out = QTermDict()
    _stream!(out, QSym[adjoint(a), a], _CNUM_ONE, _EMPTY_NE)
    @test length(out) == 1
end

@testset "_canonicalize!: aa† → a†a + 1" begin
    h = FockSpace(:f)
    a = Destroy(h, :a)
    out = QTermDict()
    _canonicalize!(out, QSym[a, adjoint(a)], _CNUM_ONE, _EMPTY_NE)
    @test length(out) == 2
end
