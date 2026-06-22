using Test
using SecondQuantizedAlgebra
using Symbolics: @variables, Num
import SecondQuantizedAlgebra: Coeff, CNum, Monomial, Poly, _to_cnum, _to_complex, to_num,
    _is_native, _is_poly, _is_symbolic_cnum, _conj_cnum, _mul_cnum, _add_cnum, _neg_cnum,
    _iszero_cnum, _CNUM_ONE, _CNUM_ZERO, _CNUM_NEG1, _CNUM_IM, _NUM_ZERO

# Coefficients carry a native `ComplexF64` fast path and a `Complex{Num}` symbolic
# fallback. These tests pin the invariants the rest of the package relies on:
# concrete numbers stay native, genuine symbols stay symbolic, materialization is
# faithful, and equality / hashing are canonical.

@testset "Coeff" begin
    @testset "native classification" begin
        for x in (0, 1, -1, 2, 1.5, 0.7, im, -im, 2 + 3im, 1 // 2, -3 // 4)
            @test _is_native(_to_cnum(x))
        end
        @test CNum === Coeff
    end

    @testset "exactness gate keeps inexact rationals / bignums symbolic" begin
        # 1//3 has no exact ComplexF64, so it must stay on the symbolic path.
        @test !_is_native(_to_cnum(1 // 3))
        # 2^70 + 1 exceeds Float64 precision (does not round-trip), so it stays symbolic;
        # an exactly representable bignum (e.g. 2^70) would correctly go native.
        @test !_is_native(_to_cnum(big(2)^70 + 1))
        @test _is_native(_to_cnum(big(2)^70))
        # round-trip is still numerically faithful
        @test _to_complex(_to_cnum(1 // 3)) ≈ 1 / 3
    end

    @testset "symbolic fallback" begin
        @variables g
        c = _to_cnum(Complex(Num(g), _NUM_ZERO))
        @test !_is_native(c)
        @test isequal(to_num(c), Complex(Num(g), _NUM_ZERO))
        # a symbolic expression that folds to a number lands back on the native path
        @test _is_native(_to_cnum(Num(g) - Num(g) + 4))
    end

    @testset "to_num round-trip and integer display" begin
        @test isequal(to_num(_to_cnum(2)), Complex(Num(2), Num(0)))
        @test isequal(to_num(_to_cnum(2)), Complex(Num(2), Num(0)))   # integer, not 2.0
        @test string(real(to_num(_to_cnum(2)))) == "2"
        @test string(real(to_num(_to_cnum(0.7)))) == "0.7"
        @test isequal(to_num(_to_cnum(2 + 3im)), Complex(Num(2), Num(3)))
    end

    @testset "_to_complex matches the value" begin
        @test _to_complex(_to_cnum(2)) === ComplexF64(2)
        @test _to_complex(_to_cnum(2 + 3im)) === ComplexF64(2, 3)
        @test _to_complex(_CNUM_ZERO) === ComplexF64(0)
    end

    @testset "conjugation and signed-zero normalization" begin
        # conj(2) produces a -0.0 imaginary part; it must normalize so structurally
        # equal coefficients stay isequal AND hash identically (dict correctness).
        c = _to_cnum(2)
        @test isequal(_conj_cnum(c), c)
        @test hash(_conj_cnum(c)) == hash(c)
        @test isequal(_neg_cnum(_neg_cnum(c)), c)
        @test hash(_neg_cnum(_neg_cnum(c))) == hash(c)
        # complex conjugation flips the sign of the imaginary part
        @test isequal(_conj_cnum(_to_cnum(2 + 3im)), _to_cnum(2 - 3im))
    end

    @testset "native arithmetic stays native and exact" begin
        @test _is_native(_mul_cnum(_to_cnum(2), _to_cnum(3)))
        @test _mul_cnum(_to_cnum(2), _to_cnum(3)) == 6
        @test _add_cnum(_to_cnum(2), _to_cnum(3)) == 5
        @test _mul_cnum(_CNUM_IM, _CNUM_IM) == -1
        @test _iszero_cnum(_add_cnum(_to_cnum(2), _neg_cnum(_to_cnum(2))))
    end

    @testset "mixed comparison with Number / Complex{Num}" begin
        @test _to_cnum(2) == 2
        @test 2 == _to_cnum(2)
        @test _to_cnum(2) == 2 + 0im
        @test _to_cnum(2) == Complex(Num(2), Num(0))
        @test isequal(_to_cnum(2), 2)
        @test !(_to_cnum(2) == 3)
    end

    @testset "scalar arithmetic operators on Coeff" begin
        @test _to_cnum(2) + _to_cnum(3) == 5
        @test _to_cnum(2) + 3 == 5
        @test 3 + _to_cnum(2) == 5
        @test _to_cnum(5) - 2 == 3
        @test -_to_cnum(2) == -2
        @test _to_cnum(2) * 3 == 6
        @test _to_cnum(6) / 2 == 3
    end

    @testset "QAdd equality is canonical across construction paths" begin
        h = FockSpace(:f)
        a = Destroy(h, :a)
        m = 2 * a' * a
        @test isequal(m', m)               # adjoint reproduces the same coefficient
        @test hash(m') == hash(m)
        @test isequal(simplify(m), m)
    end

    @testset "Poly tier" begin
        @variables ω g κ
        @variables gc::Number   # complex-symtype parameter

        @testset "recognition" begin
            # single monomials, sums, powers, and quotients all land on the native Poly path
            for x in (g, ω * g, 2 * ω * g, 0.5g, g^3, g / κ, im * g, ω + g, (g + κ)^2)
                @test _is_poly(_to_cnum(x))
            end
            # an irreducible one-argument call on an atom (`exp`, `sin`,
            # `real`/`imag` of a complex parameter, ...) is kept native as an opaque
            # Poly atom rather than escalating the whole coefficient to the symbolic path
            for x in (exp(g), sin(g), real(gc), imag(gc), conj(gc))
                @test _is_poly(_to_cnum(x))
            end
            # a radical of a single atom is a Poly with a *rational* exponent
            # (`sqrt(g) = g^(1//2)`), so radicals canonicalize on the fast path too
            for x in (sqrt(g), cbrt(g), g^(1 // 2), g^(3 // 2), g^(-1 // 2), sqrt(sqrt(g)))
                @test _is_poly(_to_cnum(x))
            end
            # genuinely non-atomic expressions stay on the symbolic (Complex{Num}) path:
            # a non-atom argument, a radical of a product, a float exponent (only exact
            # `Rational` exponents fold), or a division by a sum
            for x in (exp(g + κ), sin(g * κ), sqrt(g * κ), (g + κ)^(1 // 2), g^0.5, 1 / (g + κ))
                @test _is_symbolic_cnum(_to_cnum(x))
                @test !_is_poly(_to_cnum(x))
            end
            # inexact scalars keep the coefficient symbolic (exactness gate), never a Poly
            @test !_is_native(_to_cnum((1 // 3) * g))
            @test !_is_poly(_to_cnum((1 // 3) * g))
        end

        @testset "materialization round-trip" begin
            @test isequal(to_num(_to_cnum(2 * ω * g)), Complex(Num(2 * ω * g), Num(0)))
            @test isequal(to_num(_to_cnum(g^3)), Complex(Num(g^3), Num(0)))
            @test isequal(to_num(_to_cnum(im * g)), Complex(Num(0), Num(g)))
            @test isequal(to_num(_to_cnum((1 // 3) * g)), Complex(Num((1 // 3) * g), Num(0)))
        end

        @testset "multiply merges factors, no CAS" begin
            p = _mul_cnum(_to_cnum(ω * g), _to_cnum(ω * κ))   # -> ω^2 g κ
            @test _is_poly(p)
            @test isequal(to_num(p), Complex(Num(ω^2 * g * κ), Num(0)))
            # scalars combine; (2g)(3κ) = 6 g κ
            @test isequal(to_num(_mul_cnum(_to_cnum(2g), _to_cnum(3κ))), Complex(Num(6 * g * κ), Num(0)))
            # factors that cancel collapse back to native
            @test _is_native(_mul_cnum(_to_cnum(g), _to_cnum(1 / g)))
        end

        @testset "add stays a native Poly (no CAS escalation)" begin
            # same factor set: scalars add, one term
            s = _add_cnum(_to_cnum(2g), _to_cnum(3g))
            @test _is_poly(s) && isequal(to_num(s), Complex(Num(5g), Num(0)))
            # different factor sets: two-term Poly, still native (this is the Design C win)
            d = _add_cnum(_to_cnum(g), _to_cnum(κ))
            @test _is_poly(d) && isequal(to_num(d), Complex(Num(g + κ), Num(0)))
        end

        @testset "(sum)^n is stored in canonical expanded form" begin
            # the package-wide 'always expand' invariant extends to polynomial coefficients
            c = _to_cnum((g + κ)^2)
            @test _is_poly(c)
            @test isequal(to_num(c), Complex(Num(g^2 + 2 * g * κ + κ^2), Num(0)))
        end

        @testset "conjugation" begin
            # real parameters: conj is identity
            @test isequal(to_num(_conj_cnum(_to_cnum(ω * g))), Complex(Num(ω * g), Num(0)))
            # complex-symtype parameter: conj reaches the factor
            @test isequal(to_num(_conj_cnum(_to_cnum(gc))), Complex(Num(conj(gc)), Num(0)))
            @test isequal(to_num(_conj_cnum(_to_cnum(im * g))), Complex(Num(0), Num(-g)))
        end

        @testset "complex parameters and irreducible couplings stay native" begin
            # `@variables _::Complex` is stored as `Complex(real(_), imag(_))`; its
            # real/imag parts are recognized atoms, so the parameter stays on the Poly
            # path instead of poisoning downstream arithmetic with a Complex{Num} tail.
            @variables gv::Complex γ::Real
            @test _is_poly(_to_cnum(gv))
            @test _is_poly(_to_cnum(√(γ) * gv))
            # conjugation is native and faithful: real γ is self-conjugate, the complex
            # parameter's imaginary part flips sign.
            @test isequal(to_num(_conj_cnum(_to_cnum(gv))), conj(gv))
            @test isequal(to_num(_conj_cnum(_to_cnum(√(γ) * gv))), conj(√(γ) * gv))
        end

        @testset "radicals canonicalize via rational exponents" begin
            @variables gv::Complex γ::Real
            # squaring a radical folds to its radicand, matching Symbolics' `sqrt(g)^2 -> g`
            @test isequal(_mul_cnum(_to_cnum(sqrt(g)), _to_cnum(sqrt(g))), _to_cnum(g))
            @test isequal(_to_cnum(sqrt(g)^2), _to_cnum(g))
            @test hash(_to_cnum(sqrt(g)^2)) == hash(_to_cnum(g))
            three = _mul_cnum(_mul_cnum(_to_cnum(sqrt(g)), _to_cnum(sqrt(g))), _to_cnum(sqrt(g)))
            @test isequal(three, _to_cnum(g * sqrt(g)))
            @test isequal(three, _to_cnum(g^(3 // 2)))
            @test isequal(_mul_cnum(_to_cnum(sqrt(γ)), _to_cnum(γ)), _to_cnum(γ^(3 // 2)))
            @test isequal(_to_cnum(sqrt(sqrt(g))), _to_cnum(g^(1 // 4)))
            # sqrt(g) and g^(1//2) are the same coefficient
            @test isequal(_to_cnum(sqrt(g)), _to_cnum(g^(1 // 2)))
            for c in (_to_cnum(sqrt(g)), _to_cnum(γ^(3 // 2)), _to_cnum(g^(1 // 4)))
                @test isequal(_to_cnum(to_num(c)), c)   # materialization round-trips
            end
            @test isequal(to_num(_conj_cnum(_to_cnum(√(γ) * gv))), conj(√(γ) * gv))
            # radical of a product / scaled atom is not distributed: symbolic leaf
            @test _is_symbolic_cnum(_to_cnum(sqrt(g * κ)))
            @test _is_symbolic_cnum(_to_cnum((2g)^(1 // 2)))
            for e in (-1 // 2, -3 // 2)   # negative fractional exps use the divide branch
                c = _to_cnum(g^e)
                @test _is_poly(c) && isequal(_to_cnum(to_num(c)), c)
            end
            @test isequal(_mul_cnum(_to_cnum(sqrt(g)), _to_cnum(1 / g)), _to_cnum(g^(-1 // 2)))
            @test _is_native(_to_cnum(Num(4)^(1 // 2)))   # numeric base stays native
            @test to_num(_to_cnum(Num(4)^(1 // 2))) == 2
            # only exact `Rational` exponents fold; a float exponent (even `1/3`) stays
            # a symbolic leaf, so the `^` recognizer needs no float-to-rational guesswork
            @test _is_symbolic_cnum(_to_cnum(g^(1 / 3)))
            @test _is_symbolic_cnum(_to_cnum(g^0.5))
        end

        @testset "radical coefficients dedup in a QAdd" begin
            @variables γ::Real
            h = FockSpace(:f)
            a = Destroy(h, :a)
            @test isequal((sqrt(γ) * sqrt(γ)) * (a' * a), γ * (a' * a))
            @test hash((sqrt(γ) * sqrt(γ)) * (a' * a)) == hash(γ * (a' * a))
            @test isequal(sqrt(γ) * sqrt(γ) * a + γ * a, 2γ * a)
        end

        @testset "negation and exact cancellation" begin
            @test isequal(to_num(_neg_cnum(_to_cnum(2 * ω * g))), Complex(Num(-2 * ω * g), Num(0)))
            @test _iszero_cnum(_add_cnum(_to_cnum(ω * g), _neg_cnum(_to_cnum(ω * g))))
        end

        @testset "array-indexed parameters are recognized" begin
            @variables ωs[1:3] gs[1:3]
            @test _is_poly(_to_cnum(ωs[1]))
            @test _is_poly(_to_cnum(ωs[1] * gs[2]))
            @test isequal(to_num(_to_cnum(ωs[1] * gs[2])), Complex(Num(ωs[1] * gs[2]), Num(0)))
            # products merge array-indexed factors like any other atom
            p = _mul_cnum(_to_cnum(ωs[1] * gs[2]), _to_cnum(ωs[1]))   # -> ωs[1]^2 gs[2]
            @test isequal(to_num(p), Complex(Num(ωs[1]^2 * gs[2]), Num(0)))
        end

        @testset "poly coefficients in a QAdd" begin
            h = FockSpace(:f)
            a = Destroy(h, :a)
            # symbolic-coefficient product builds and substitutes correctly
            expr = (ω * g) * (a' * a)
            @test isequal(prefactor(expr), Complex(Num(ω * g), Num(0)))
            sub = substitute(expr, Dict(ω => 2, g => 3))
            @test isequal(prefactor(sub), Complex(Num(6), Num(0)))
        end
    end
end
