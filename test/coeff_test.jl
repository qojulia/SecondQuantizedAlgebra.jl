using Test
using SecondQuantizedAlgebra
using Symbolics: Symbolics, @variables, Num
using SymbolicUtils: SymbolicUtils
import SecondQuantizedAlgebra: Coeff, CNum, Monomial, Poly, _to_cnum, _to_complex, to_num,
    _is_native, _is_poly, _is_symbolic_cnum, _conj_cnum, _mul_cnum, _add_cnum, _neg_cnum,
    _iszero_cnum, _CNUM_ONE, _CNUM_ZERO, _CNUM_NEG1, _CNUM_IM, _NUM_ZERO, expim,
    _phase_coeff, _is_phase, exponential_form, trigonometric_form

function _phase_allocations(f)
    f()
    return @allocations f()
end

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

        @testset "wide sums fold in one pass (batched _rec_sum)" begin
            @variables p[1:6]
            pv = [SecondQuantizedAlgebra.SymbolicUtils.unwrap(Num(p[k])) for k in 1:6]
            pairwise(terms) = foldl(
                (c, a) -> _add_cnum(c, _to_cnum(a)), terms; init = _CNUM_ZERO,
            )

            # distinct monomials: a 6-term Poly, matching the pairwise fold exactly
            distinct = [Num(k) * p[k] for k in 1:6]
            cd = _to_cnum(sum(distinct))
            @test _is_poly(cd) && length(cd.tail.terms) == 6
            @test isequal(cd, pairwise(distinct))

            # repeated factor sets coalesce: 1·g + 2·g + 3·g -> 6g (one term)
            cc = _to_cnum(g + 2g + 3g)
            @test _is_poly(cc) && length(cc.tail.terms) == 1
            @test isequal(to_num(cc), Complex(Num(6g), Num(0)))

            # native + poly mixed in one sum: constants collapse to a single term
            cm = _to_cnum(2 + g + 3 + κ)
            @test isequal(to_num(cm), Complex(Num(5 + g + κ), Num(0)))

            # full cancellation across many terms -> exact zero
            @test _iszero_cnum(_to_cnum(g + κ - g - κ))

            # a sym-leaf term (irreducible call) drops to the fallback add and the
            # polynomial part still coalesces around it
            cs = _to_cnum(g + sin(ω) + g)
            @test isequal(to_num(cs), Complex(Num(2g + sin(ω)), Num(0)))
        end

        @testset "wide products fold in one pass (batched _rec_prod)" begin
            # A `*`-headed coefficient of single-monomial factors collapses to one
            # monomial; the factor lists merge in a single pass, not pairwise.
            # Result must be byte-identical to the pairwise `_mul_cnum` fold.
            @variables r[1:6]
            pairwise(fs) = foldl(
                (c, a) -> _mul_cnum(c, _to_cnum(a)), fs; init = _CNUM_ONE,
            )

            # distinct atoms: one monomial with 6 factors, matching the pairwise fold
            facs = [Num(r[k]) for k in 1:6]
            cp = _to_cnum(prod(facs))
            @test _is_poly(cp) && length(cp.tail.terms) == 1
            @test length(cp.tail.terms[1].syms) == 6
            @test isequal(cp, pairwise(facs))

            # scalars multiply, exponents accumulate across repeated factors
            ce = _to_cnum(2g * 3g * κ)            # 6 g^2 κ
            @test isequal(to_num(ce), Complex(Num(6 * g^2 * κ), Num(0)))

            # a single-monomial part times a multi-term factor distributes (intrinsic)
            cdist = _to_cnum(g * κ * (g + κ))
            @test isequal(to_num(cdist), Complex(Num(g^2 * κ + g * κ^2), Num(0)))

            # factor cancellation drops the atom entirely
            @test isequal(to_num(_to_cnum(g * κ / g)), Complex(Num(κ), Num(0)))

            # a zero factor collapses the whole product to native zero
            @test _iszero_cnum(_to_cnum(g * κ * 0 * ω))

            # a sym-leaf factor drops to the fallback multiply; atoms still merge
            csl = _to_cnum(g * exp(κ) * g)
            @test isequal(to_num(csl), Complex(Num(g^2 * exp(κ)), Num(0)))
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

        @testset "conjugation involution folds (regression #7cc3ad7)" begin
            # `conj` is an involution: conjugating twice must return the original
            # factor, not nest a `conj(conj(x))` that never folds and survives
            # downstream.
            c = _to_cnum(gc)
            @test isequal(_conj_cnum(_conj_cnum(c)), c)
            @test isequal(to_num(_conj_cnum(_conj_cnum(c))), Complex(Num(gc), Num(0)))
            @test hash(_conj_cnum(_conj_cnum(c))) == hash(c)
            # double-conjugating a scaled complex factor also returns the original
            cs = _to_cnum(g * gc)
            @test isequal(_conj_cnum(_conj_cnum(cs)), cs)
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

        @testset "elementary functions of a literal zero fold to their exact value" begin
            @variables ω t
            # Euler-expanding exp(im*ω*t) leaves exp(0) factors; they must fold away.
            c = _to_cnum(Complex{Num}(exp(1.0im * ω * t)))
            @test !occursin("exp(0)", repr(to_num(c)))
            @test isequal(to_num(c), Complex(Num(cos(t * ω)), Num(sin(t * ω))))
            # exact identities at 0 fold to native
            @test _is_native(_to_cnum(Num(exp(Num(0)))))
            @test _to_complex(_to_cnum(Num(exp(Num(0))))) ≈ 1
            @test _to_complex(_to_cnum(Num(cos(Num(0))))) ≈ 1
            @test _to_complex(_to_cnum(Num(sin(Num(0))))) ≈ 0
            # symbol / non-zero-constant / irrational args stay symbolic (exactness gate)
            @test !_is_native(_to_cnum(Num(exp(ω))))
            @test !_is_native(_to_cnum(Num(exp(Num(2)))))
            @test !_is_native(_to_cnum(Num(exp(Num(0.1)))))
            @test !_is_native(_to_cnum(Num(cos(Num(1)))))
            @test !_is_native(_to_cnum(cos(Num(π))))
            @test !_is_native(_to_cnum(sin(Num(π))))
        end
    end

    @testset "unit phases" begin
        @variables ω J t θ
        p = _phase_coeff(ω * t)
        @test _is_poly(p)
        @test _is_phase(p.tail.terms[1].syms[1])
        # conjugation is exponent negation, so a phase cancels against its own conjugate
        @test _mul_cnum(p, _conj_cnum(p)) == _CNUM_ONE
        pp = _mul_cnum(p, p)
        @test _mul_cnum(pp, _conj_cnum(pp)) == _CNUM_ONE
        # `expim(0)` is 1, not an atom
        @test _phase_coeff(0) == _CNUM_ONE
        @test _to_cnum(expim(Num(0))) == _CNUM_ONE
        # one atom per frequency however the argument is spelled, and whatever its sign
        @test isequal(_phase_coeff((ω + 2J) * t), _phase_coeff(ω * t + 2 * J * t))
        @test _mul_cnum(_phase_coeff(-ω * t), _phase_coeff(ω * t)) == _CNUM_ONE
        # commensurate frequencies are exponent arithmetic, not an angle-addition identity
        @variables E1 E2 E3
        ph(i, j) = _mul_cnum(_phase_coeff(i * t), _conj_cnum(_phase_coeff(j * t)))
        @test _iszero_cnum(
            _add_cnum(_mul_cnum(ph(E1, E2), ph(E2, E3)), _neg_cnum(ph(E1, E3)))
        )
        # lowering round-trips back onto the same atom
        @test isequal(p, _to_cnum(to_num(p)))
        @test isequal(_conj_cnum(p), _to_cnum(to_num(_conj_cnum(p))))
        # a numeric argument still evaluates
        @test (@inferred expim(0.5)) ≈ exp(0.5im)
        @test (@inferred expim(ω * t)) isa Coeff
        @test_throws ArgumentError expim(1 + 0im)
        @test_throws ArgumentError expim(1 + 2im)
        @variables z::Complex
        @test_throws ArgumentError expim(z)
        @test_throws ArgumentError expim(SymbolicUtils.unwrap(z))
    end

    @testset "canonical phase group" begin
        @variables x y z g h ω t
        px = expim(x)
        py = expim(y)

        @test (@inferred px * py) == expim(x + y)
        @test px * py == py * px
        @test px * px == expim(2x)
        @test px * py * expim(z) == expim(x + y + z)
        @test (@inferred px / py) == expim(x - y)
        @test (@inferred px^3) == expim(3x)
        @test (@inferred px^(-3)) == expim(-3x)
        @test inv(px) == conj(px) == expim(-x)
        @test px * inv(px) == _CNUM_ONE
        @test expim(x) * expim(y - x) == expim(y)
        @test expim((ω + 2g) * t) * expim(-2g * t) == expim(ω * t)

        combined = g * px * py
        combined_term = only((combined.tail::Poly).terms)
        @test count(_is_phase, combined_term.syms) == 1
        @test only(combined_term.exps[findall(_is_phase, combined_term.syms)]) in (-1, 1)

        hs = FockSpace(:phase_collection)
        a = Destroy(hs, :a)
        @test iszero(px * py * a - expim(x + y) * a)

        @test_throws MethodError px^(1 // 2)
        @test_throws MethodError sqrt(px)
        phase_symbol = only(only((px.tail::Poly).terms).syms)
        phase_root = SymbolicUtils.term(sqrt, phase_symbol; type = Complex{Real})
        @test_throws ArgumentError _to_cnum(phase_root)

        cg, ch = _to_cnum(g), _to_cnum(h)
        @test _phase_allocations(() -> cg * ch) <= 10
        @test _phase_allocations(() -> px * conj(px)) <= 18
        @test _phase_allocations(() -> px * px) <= 165
        @test _phase_allocations(() -> px * py) <= 285
    end

    @testset "phase substitution, calculus, and projections" begin
        @variables ω t θ g z::Complex
        p = expim(ω * t)

        @test substitute(p, Dict()) === p
        @test (@inferred substitute(p, Dict(ω => 2ω))) == expim(2ω * t)
        @test substitute(p, Dict(t => 0)) == _CNUM_ONE
        @test substitute(p, Dict(g => 2g)) === p
        @test _to_complex(substitute(p, Dict(ω => 2.0, t => 0.25))) ≈ exp(0.5im)
        @test_throws ArgumentError substitute(p, Dict(ω => z))
        @test_throws ArgumentError substitute(p, Dict(ω => 1 + 2im))
        @test_throws ArgumentError substitute(p, Dict(ω => 1 + 0im))

        D = Symbolics.Differential(t)
        differentiated = @inferred Symbolics.derivative(p, t)
        @test differentiated == im * ω * p
        @test _to_cnum(Symbolics.expand_derivatives(D(p))) == differentiated

        phase = expim(θ)
        @test isequal(@inferred(real(phase)), cos(θ))
        @test isequal(@inferred(imag(phase)), sin(θ))
        @test isequal(@inferred(real(conj(phase))), cos(θ))
        @test isequal(@inferred(imag(conj(phase))), -sin(θ))
        @test isequal(real(im * phase), -sin(θ))
        @test isequal(imag(im * phase), cos(θ))
        @test (@inferred abs(phase)) == 1
        @test (@inferred abs2(phase)) == 1
        @test abs(conj(phase)) == 1
        @test abs(-im * phase) == 1
        @test abs2(-im * phase) == 1
        @test_throws MethodError abs(_to_cnum(g))
        @test_throws MethodError abs2(_to_cnum(g))
    end

    @testset "explicit phase representations" begin
        @variables θ ω t z::Complex
        cosine_phase = @inferred exponential_form(cos(θ))
        sine_phase = @inferred exponential_form(sin(θ))
        composite_phase = @inferred exponential_form(cos(ω * t) + sin(ω * t))
        intact_exp = SymbolicUtils.term(
            exp, im * SymbolicUtils.unwrap(θ); type = Complex{Real},
        )
        intact_cis = SymbolicUtils.term(
            cis, SymbolicUtils.unwrap(θ); type = Complex{Real},
        )
        arbitrary_exp = SymbolicUtils.term(
            exp, SymbolicUtils.unwrap(z); type = Complex{Real},
        )

        @test cosine_phase == (expim(θ) + expim(-θ)) / 2
        @test sine_phase == (expim(θ) - expim(-θ)) / (2im)
        @test trigonometric_form(cosine_phase) == _to_cnum(cos(θ))
        @test trigonometric_form(sine_phase) == _to_cnum(sin(θ))
        @test trigonometric_form(composite_phase) == _to_cnum(cos(ω * t) + sin(ω * t))
        @test trigonometric_form(expim(θ)^2) == _to_cnum(cos(2θ) + im * sin(2θ))
        @test (@inferred exponential_form(exp(im * θ))) == expim(θ)
        @test (@inferred exponential_form(cis(θ))) == expim(θ)
        @test (@inferred exponential_form(intact_exp)) == expim(θ)
        @test (@inferred exponential_form(intact_cis)) == expim(θ)
        @test isequal(exponential_form(arbitrary_exp), _to_cnum(arbitrary_exp))
        @test exponential_form(exp(θ)) == _to_cnum(exp(θ))

        h = FockSpace(:phase_forms)
        a = Destroy(h, :a)
        q = (cos(ω * t) + sin(ω * t)) * a + cos(θ) * a' * a
        phase_q = @inferred exponential_form(q)
        @test Set(keys(phase_q.arguments)) == Set(keys(q.arguments))
        @test iszero(simplify(trigonometric_form(phase_q) - q))

        nested_phase = expim(cos(θ))
        @test exponential_form(nested_phase) == nested_phase
        @test exponential_form(3) == 3
        @test trigonometric_form(3) == 3
    end

    # The public `expim` used to return a `Num`. `Base.conj(::Num)` is the identity and
    # `Complex * Num` splits into real/imag halves, so every law below silently failed.
    @testset "public expim is a coefficient" begin
        @variables ω t
        h = FockSpace(:f)
        a = Destroy(h, :a)
        p = expim(ω * t)
        @test p isa Coeff

        @test conj(p) * p == _CNUM_ONE
        @test (@inferred inv(p)) == conj(p)
        @test (@inferred p^(-1)) == conj(p)
        @test !isequal(conj(p), p)
        @test (im * p) * conj(p) == _CNUM_IM
        @test _is_poly(im * p)
        @test isequal(expim(ω * t) * expim(-ω * t) * a, 1 * a)
        @test isequal(p * a, expim(ω * t) * a)

        # a phase over a literal is its value, so numerically cancelling terms fold
        @test _is_native(expim(Num(1.0)))
        @test iszero(expim(Num(1.0)) * a - exp(im * 1.0) * a)

        # `type`/`shape` take part in hash-consing, so an unregistered head would come back
        # from any rebuild as a different atom whose `conj` is the identity
        b = SecondQuantizedAlgebra._expim(SymbolicUtils.unwrap(ω * t))
        rebuilt = SymbolicUtils.maketerm(
            typeof(b), SymbolicUtils.operation(b), SymbolicUtils.arguments(b),
            SymbolicUtils.metadata(b),
        )
        @test rebuilt === b
        @test SymbolicUtils.symtype(b) === Complex{Real}
        @test Symbolics.substitute(b, Dict(SymbolicUtils.unwrap(ω) => 2.0)) ===
            SecondQuantizedAlgebra._expim(SymbolicUtils.unwrap(2.0 * t))

        # conjugation negates the argument rather than conjugating it
        @test isequal(
            SecondQuantizedAlgebra.qadjoint(Num(b)),
            SecondQuantizedAlgebra._expim(SymbolicUtils.unwrap(-ω * t)),
        )

        # a registered derivative, not an inert `Differential` node
        d = Symbolics.expand_derivatives(Symbolics.Differential(t)(Num(b)))
        @test isequal(_to_cnum(d), _mul_cnum(_mul_cnum(_CNUM_IM, _to_cnum(ω)), p))
    end
end
