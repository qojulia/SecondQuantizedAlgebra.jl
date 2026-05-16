using SecondQuantizedAlgebra
using Test
using JET
using Symbolics: @variables

import SecondQuantizedAlgebra: CNum, QSym,
    _mul_cnum, _add_cnum, _neg_cnum, _to_cnum, _iszero_num, _iszero_cnum,
    _CNUM_ZERO, _CNUM_ONE, _CNUM_NEG1, _CNUM_IM, _CNUM_NEG_IM,
    _NUM_ZERO, _NUM_ONE,
    _site_compare, _can_commute, _commute_pair, _reduce_pair, _ground_state_expand,
    _EMPTY_NE,
    SiteCmp, ReduceKind

# Type-stability tests for the canonicalization hot path.
#
# `Vector{QSym}` stores abstract elements, so 6 root dispatches (`_site_compare`,
# `_can_commute`, `_commute_pair`, `_reduce_pair` × 2 passes) are intrinsic to
# the design. Everything ELSE on the hot path is type-stable and the tests below
# pin that down so regressions show up immediately. The cnum arithmetic in
# particular is fully concrete: every branch of `_mul_cnum`, `_add_cnum`,
# `_neg_cnum` returns `Complex{Num}` with no `Any` leakage.

@testset "Type Stability" begin

    @testset "CNum arithmetic returns CNum (no Any leakage)" begin
        @test @inferred(_mul_cnum(_CNUM_ONE, _CNUM_ONE)) isa CNum
        @test @inferred(_mul_cnum(_CNUM_ONE, _CNUM_ZERO)) isa CNum
        @test @inferred(_mul_cnum(_CNUM_NEG1, _CNUM_IM)) isa CNum
        @test @inferred(_add_cnum(_CNUM_ONE, _CNUM_ONE)) isa CNum
        @test @inferred(_add_cnum(_CNUM_ONE, _CNUM_NEG1)) isa CNum
        @test @inferred(_neg_cnum(_CNUM_ONE)) isa CNum
        @test @inferred(_neg_cnum(_CNUM_IM)) isa CNum

        # With symbolic prefactors the slow path must also stay typed.
        @variables x y
        c1 = Complex(SecondQuantizedAlgebra.Num(x), _NUM_ZERO)
        c2 = Complex(SecondQuantizedAlgebra.Num(y), _NUM_ZERO)
        @test @inferred(_mul_cnum(c1, c2)) isa CNum
        @test @inferred(_add_cnum(c1, c2)) isa CNum
        @test @inferred(_neg_cnum(c1)) isa CNum
    end

    @testset "_to_cnum is type-stable for numeric inputs" begin
        @test @inferred(_to_cnum(1)) isa CNum
        @test @inferred(_to_cnum(0)) isa CNum
        @test @inferred(_to_cnum(-1)) isa CNum
        @test @inferred(_to_cnum(1.5)) isa CNum
        @test @inferred(_to_cnum(im)) isa CNum
        @test @inferred(_to_cnum(-im)) isa CNum
        @test @inferred(_to_cnum(1 + 2im)) isa CNum
    end

    @testset "Operator hook return types are concrete" begin
        h = FockSpace(:f)
        a = Destroy(h, :a)
        ad = Create(h, :a)

        # Same-type dispatch: fully concrete, no Vector{QSym} indirection.
        @test @inferred(_site_compare(a, ad, _EMPTY_NE)) isa SiteCmp
        @test @inferred(_can_commute(a, ad)) isa Bool
        @test @inferred(_can_commute(ad, a)) isa Bool
        @test @inferred(_commute_pair(a, ad)) isa Tuple{QSym, QSym, CNum, Vector{QSym}}

        ha = NLevelSpace(:atom, 3)
        σ12 = Transition(ha, :σ, 1, 2)
        σ21 = Transition(ha, :σ, 2, 1)
        @test @inferred(_reduce_pair(σ12, σ21)) isa Tuple{ReduceKind, QSym, CNum}
        @test @inferred(_ground_state_expand(σ12)) isa Tuple{Bool, Int, Int, Int}

        hp = PauliSpace(:s)
        px = Pauli(hp, :σ, 1)
        py = Pauli(hp, :σ, 2)
        @test @inferred(_reduce_pair(px, py)) isa Tuple{ReduceKind, QSym, CNum}

        hs = SpinSpace(:s)
        sx = Spin(hs, :S, 1)
        sy = Spin(hs, :S, 2)
        @test @inferred(_commute_pair(sy, sx)) isa Tuple{QSym, QSym, CNum, Vector{QSym}}
    end

    @testset "Public arithmetic on Fock operators" begin
        h = FockSpace(:f)
        a = Destroy(h, :a)
        ad = Create(h, :a)

        # Each of these is a hot-path entry; the outer call returns QAdd.
        @test (a * ad) isa SecondQuantizedAlgebra.QAdd
        @test (ad * a) isa SecondQuantizedAlgebra.QAdd
        @test (a + ad) isa SecondQuantizedAlgebra.QAdd
        @test (a - ad) isa SecondQuantizedAlgebra.QAdd
        @test (a^2) isa SecondQuantizedAlgebra.QAdd
        @test (2 * a) isa SecondQuantizedAlgebra.QAdd
        @test commutator(a, ad) isa SecondQuantizedAlgebra.QAdd
    end

    @testset "Public arithmetic on NLevel operators" begin
        h = NLevelSpace(:atom, 3)
        σ12 = Transition(h, :σ, 1, 2)
        σ21 = Transition(h, :σ, 2, 1)
        @test (σ12 * σ21) isa SecondQuantizedAlgebra.QAdd
        @test (σ21 * σ12) isa SecondQuantizedAlgebra.QAdd
        @test (σ12 + σ21) isa SecondQuantizedAlgebra.QAdd
    end

    @testset "Public arithmetic on Pauli/Spin" begin
        hp = PauliSpace(:s)
        px = Pauli(hp, :σ, 1)
        py = Pauli(hp, :σ, 2)
        @test (px * py) isa SecondQuantizedAlgebra.QAdd
        @test (px + py) isa SecondQuantizedAlgebra.QAdd

        hs = SpinSpace(:s)
        sx = Spin(hs, :S, 1)
        sy = Spin(hs, :S, 2)
        @test (sx * sy) isa SecondQuantizedAlgebra.QAdd
        @test commutator(sx, sy) isa SecondQuantizedAlgebra.QAdd
    end

    @testset "Public arithmetic on PhaseSpace" begin
        hp = PhaseSpace(:osc)
        x = Position(hp, :x)
        p = Momentum(hp, :p)
        @test (x * p) isa SecondQuantizedAlgebra.QAdd
        @test (p * x) isa SecondQuantizedAlgebra.QAdd
        @test commutator(x, p) isa SecondQuantizedAlgebra.QAdd
    end

    @static if isempty(VERSION.prerelease)
        @testset "JET report_package (correctness)" begin
            # report_package finds actual type errors (method-not-found, etc.).
            # `ignore_missing_comparison = true` follows QuantumToolbox.jl's pattern —
            # silences `Missing` union-split warnings that bubble up from Symbolics.
            # report_opt is intentionally omitted: 4 root dispatches + 7 cascade from
            # `Vector{QSym}` abstract elements are intrinsic to the design.
            result = JET.report_package(
                SecondQuantizedAlgebra;
                target_modules = (SecondQuantizedAlgebra,),
                ignore_missing_comparison = true,
            )
            @test isempty(JET.get_reports(result))
        end

        @testset "JET report_call on entry points (no errors)" begin
            # @report_call catches actual MethodError/UndefVar issues at runtime.
            # We exercise each operator family's hot path. `a^n` is skipped because
            # JET descends through Symbolics' `Num^Integer` and surfaces hundreds of
            # macro-generated UndefVarError reports from SymbolicUtils' BasicSymbolic
            # internals — outside our control.
            hf = FockSpace(:f); a = Destroy(hf, :a); ad = Create(hf, :a)
            hn = NLevelSpace(:atom, 3); σ12 = Transition(hn, :σ, 1, 2); σ21 = Transition(hn, :σ, 2, 1)
            hp = PauliSpace(:s); px = Pauli(hp, :σ, 1); py = Pauli(hp, :σ, 2)
            hs = SpinSpace(:s); sx = Spin(hs, :S, 1); sy = Spin(hs, :S, 2)
            hph = PhaseSpace(:osc); xx = Position(hph, :x); pp = Momentum(hph, :p)

            for (name, expr) in [
                    ("a' * a", () -> ad * a),
                    ("σ12 * σ21", () -> σ12 * σ21),
                    ("px * py", () -> px * py),
                    ("sx * sy", () -> sx * sy),
                    ("x * p", () -> xx * pp),
                    ("commutator(a, a')", () -> commutator(a, ad)),
                    ("normal_order(a'*a*a')", () -> normal_order(ad * a * ad)),
                ]
                rep = JET.@report_call target_modules=(SecondQuantizedAlgebra,) ignore_missing_comparison=true expr()
                @testset "$name" begin
                    @test isempty(JET.get_reports(rep))
                end
            end
        end
    end
end
