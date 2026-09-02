using Test
using SecondQuantizedAlgebra
using Symbolics: Symbolics, @variables, Num
using SymbolicUtils: SymbolicUtils
import SecondQuantizedAlgebra: Coeff, CNum, Monomial, Poly, to_cnum, to_complex, to_num,
    is_native, is_poly, is_symbolic_cnum, conj_cnum, mul_cnum, add_cnum, neg_cnum,
    iszero_cnum, CNUM_ONE, CNUM_ZERO, CNUM_NEG1, CNUM_IM, NUM_ZERO, expim,
    phase_coeff, is_phase, exponential_form, trigonometric_form, phase_terms, from_raw,
    CNUM_HALF, cnum_sym, raw_parts

function phase_allocations(f)
    f()
    return @allocations f()
end

# Coefficients carry a native `ComplexF64` fast path and a `Complex{Num}` symbolic
# fallback. These tests pin the invariants the rest of the package relies on:
# concrete numbers stay native, genuine symbols stay symbolic, materialization is
# faithful, and equality / hashing are canonical.

@testset "Coeff" begin
    @testset "native classification" begin
        for x in (0, 1, -1, 2, 1.5, 0.7, im, -im, 2 + 3im)
            @test is_native(to_cnum(x))
        end
        @test !is_native(to_cnum(1 // 2))
        @test !is_native(to_cnum(-3 // 4))
        @test CNum === Coeff
    end

    @testset "exact rational coefficients avoid Float64" begin
        # Every rational literal stays exact, including ones that happen to round-trip
        # through Float64 without loss.
        @test is_poly(to_cnum(1 // 4))
        @test isequal(to_num(to_cnum(1 // 4)), Complex(Num(1 // 4), Num(0)))
        exact_complex = to_cnum(Complex(1 // 4, 1 // 2))
        @test is_poly(exact_complex)
        @test isequal(to_num(exact_complex), Complex(Num(1 // 4), Num(1 // 2)))
        @test is_native(to_cnum(0.25))
        @test to_complex(to_cnum(1 // 4)) === ComplexF64(0.25)
        # A rational coefficient remains exact when multiplied by a symbolic parameter.
        @variables r
        @test is_poly(to_cnum((1 // 4) * r))
        @test occursin("1//4", string(to_num(to_cnum((1 // 4) * r))))
        @test occursin("1//2", string(to_num(to_cnum(r / 2))))
    end

    @testset "exactness gate keeps large values symbolic" begin
        # 1//3 remains exact as well.
        @test !is_native(to_cnum(1 // 3))
        # An exactly integral Gaussian rational still needs a round-trip check: Int-sized
        # values above 2^53 are not all exactly representable by ComplexF64.
        large_exact = Complex((2^53 + 1) // 1, 0 // 1)
        large_coeff = to_cnum(large_exact)
        @test !is_native(large_coeff)
        @test isequal(to_num(large_coeff), Complex(Num(2^53 + 1), Num(0)))
        # 2^70 + 1 exceeds Float64 precision (does not round-trip), so it stays symbolic;
        # an exactly representable bignum (e.g. 2^70) would correctly go native.
        @test !is_native(to_cnum(big(2)^70 + 1))
        @test is_native(to_cnum(big(2)^70))
        # round-trip is still numerically faithful
        @test to_complex(to_cnum(1 // 3)) ≈ 1 / 3
    end

    @testset "raw numeric constants fold back to native" begin
        raw_const = SymbolicUtils.term(
            complex, 3, 0; type = Complex{Real}, shape = UnitRange{Int}[],
        )
        folded = from_raw(raw_const)
        @test is_native(folded)
        @test to_num(folded) == Complex(Num(3), Num(0))
    end

    @testset "raw Number atoms retain complex projections" begin
        @variables z::Number
        raw_atom = from_raw(z; normalize = false)
        @test SymbolicUtils._iszero(real(raw_atom) - real(z))
        @test SymbolicUtils._iszero(imag(raw_atom) - imag(z))
        materialized = to_num(raw_atom)
        @test SymbolicUtils._iszero(real(materialized) - real(z))
        @test SymbolicUtils._iszero(imag(materialized) - imag(z))
    end

    @testset "internal half stays on the native fast path" begin
        @test is_native(CNUM_HALF)
    end

    @testset "symbolic fallback" begin
        @variables g
        c = to_cnum(Complex(Num(g), NUM_ZERO))
        @test !is_native(c)
        @test isequal(to_num(c), Complex(Num(g), NUM_ZERO))
        # a symbolic expression that folds to a number lands back on the native path
        @test is_native(to_cnum(Num(g) - Num(g) + 4))
    end

    @testset "to_num round-trip and integer display" begin
        @test isequal(to_num(to_cnum(2)), Complex(Num(2), Num(0)))
        @test isequal(to_num(to_cnum(2)), Complex(Num(2), Num(0)))   # integer, not 2.0
        @test string(real(to_num(to_cnum(2)))) == "2"
        @test string(real(to_num(to_cnum(0.7)))) == "0.7"
        @test isequal(to_num(to_cnum(2 + 3im)), Complex(Num(2), Num(3)))
    end

    @testset "to_complex matches the value" begin
        @test to_complex(to_cnum(2)) === ComplexF64(2)
        @test to_complex(to_cnum(2 + 3im)) === ComplexF64(2, 3)
        @test to_complex(CNUM_ZERO) === ComplexF64(0)
    end

    @testset "conjugation and signed-zero normalization" begin
        # conj(2) produces a -0.0 imaginary part; it must normalize so structurally
        # equal coefficients stay isequal AND hash identically (dict correctness).
        c = to_cnum(2)
        @test isequal(conj_cnum(c), c)
        @test hash(conj_cnum(c)) == hash(c)
        @test isequal(neg_cnum(neg_cnum(c)), c)
        @test hash(neg_cnum(neg_cnum(c))) == hash(c)
        # complex conjugation flips the sign of the imaginary part
        @test isequal(conj_cnum(to_cnum(2 + 3im)), to_cnum(2 - 3im))
    end

    @testset "native arithmetic stays native and exact" begin
        @test is_native(mul_cnum(to_cnum(2), to_cnum(3)))
        @test mul_cnum(to_cnum(2), to_cnum(3)) == 6
        @test add_cnum(to_cnum(2), to_cnum(3)) == 5
        @test mul_cnum(CNUM_IM, CNUM_IM) == -1
        @test iszero_cnum(add_cnum(to_cnum(2), neg_cnum(to_cnum(2))))
    end

    @testset "mixed comparison with Number / Complex{Num}" begin
        @test to_cnum(2) == 2
        @test 2 == to_cnum(2)
        @test to_cnum(2) == 2 + 0im
        @test to_cnum(2) == Complex(Num(2), Num(0))
        @test isequal(to_cnum(2), 2)
        @test !(to_cnum(2) == 3)
    end

    @testset "scalar arithmetic operators on Coeff" begin
        @test to_cnum(2) + to_cnum(3) == 5
        @test to_cnum(2) + 3 == 5
        @test 3 + to_cnum(2) == 5
        @test to_cnum(5) - 2 == 3
        @test -to_cnum(2) == -2
        @test to_cnum(2) * 3 == 6
        @test to_cnum(6) / 2 == 3
    end

    @testset "QAdd equality is canonical across construction paths" begin
        h = FockSpace(:f)
        a = Destroy(h, :a)
        m = 2 * a' * a
        @test isequal(m', m)               # adjoint reproduces the same coefficient
        @test hash(m') == hash(m)
        @test isequal(simplify(m), m)

        @variables g k
        simplified_fraction = simplify(((g + g * k) / (1 + k)) * a)
        @test isequal(simplified_fraction, g * a)
        @test hash(simplified_fraction) == hash(g * a)
    end

    @testset "Poly tier" begin
        @variables ω g κ
        @variables gc::Number   # complex-symtype parameter

        @testset "recognition" begin
            # single monomials, sums, powers, and quotients all land on the native Poly path
            for x in (g, ω * g, 2 * ω * g, 0.5g, g^3, g / κ, im * g, ω + g, (g + κ)^2)
                @test is_poly(to_cnum(x))
            end
            # an irreducible one-argument call on an atom (`exp`, `sin`,
            # `real`/`imag` of a complex parameter, ...) is kept native as an opaque
            # Poly atom rather than escalating the whole coefficient to the symbolic path
            for x in (exp(g), sin(g), real(gc), imag(gc), conj(gc))
                @test is_poly(to_cnum(x))
            end
            # a radical of a single atom is a Poly with a *rational* exponent
            # (`sqrt(g) = g^(1//2)`), so radicals canonicalize on the fast path too
            for x in (sqrt(g), cbrt(g), g^(1 // 2), g^(3 // 2), g^(-1 // 2), sqrt(sqrt(g)))
                @test is_poly(to_cnum(x))
            end
            # genuinely non-atomic expressions stay on the symbolic (Complex{Num}) path:
            # a non-atom argument, a radical of a product, a float exponent (only exact
            # `Rational` exponents fold), or a division by a sum
            for x in (exp(g + κ), sin(g * κ), sqrt(g * κ), (g + κ)^(1 // 2), g^0.5, 1 / (g + κ))
                @test is_symbolic_cnum(to_cnum(x))
                @test !is_poly(to_cnum(x))
            end
            # Exact rational scalars remain on the polynomial path.
            @test !is_native(to_cnum((1 // 3) * g))
            @test is_poly(to_cnum((1 // 3) * g))
        end

        @testset "materialization round-trip" begin
            @test isequal(to_num(to_cnum(2 * ω * g)), Complex(Num(2 * ω * g), Num(0)))
            @test isequal(to_num(to_cnum(g^3)), Complex(Num(g^3), Num(0)))
            @test isequal(to_num(to_cnum(im * g)), Complex(Num(0), Num(g)))
            @test isequal(to_num(to_cnum((1 // 3) * g)), Complex(Num((1 // 3) * g), Num(0)))

            # Number-symtype atoms occupy the real coefficient slot even when mixed
            # arithmetic sends the coefficient through the raw symbolic tier.
            @test isequal(real(to_cnum(gc)), gc)
            @test iszero(imag(to_cnum(gc)))
            mixed = to_cnum(gc) / to_cnum(sqrt(g * κ))
            @test isequal(to_num(mixed), Complex(Num(gc / sqrt(g * κ)), Num(0)))

            # `Symbolics.IM` is the actual imaginary unit, not a Number-symtype
            # coefficient slot. It must take the native path before generic atoms.
            symbolic_im = to_cnum(Symbolics.IM)
            @test isequal(symbolic_im, CNUM_IM)
            @test iszero(real(symbolic_im))
            @test isequal(imag(symbolic_im), Num(1))
        end

        @testset "multiply merges factors, no CAS" begin
            p = mul_cnum(to_cnum(ω * g), to_cnum(ω * κ))   # -> ω^2 g κ
            @test is_poly(p)
            @test isequal(to_num(p), Complex(Num(ω^2 * g * κ), Num(0)))
            # scalars combine; (2g)(3κ) = 6 g κ
            @test isequal(to_num(mul_cnum(to_cnum(2g), to_cnum(3κ))), Complex(Num(6 * g * κ), Num(0)))
            # factors that cancel collapse back to native
            @test is_native(mul_cnum(to_cnum(g), to_cnum(1 / g)))
        end

        @testset "add stays a native Poly (no CAS escalation)" begin
            # same factor set: scalars add, one term
            s = add_cnum(to_cnum(2g), to_cnum(3g))
            @test is_poly(s) && isequal(to_num(s), Complex(Num(5g), Num(0)))
            # different factor sets: two-term Poly, still native (this is the Design C win)
            d = add_cnum(to_cnum(g), to_cnum(κ))
            @test is_poly(d) && isequal(to_num(d), Complex(Num(g + κ), Num(0)))
        end

        @testset "wide sums fold in one pass (batched rec_sum)" begin
            @variables p[1:6]
            pv = [SecondQuantizedAlgebra.SymbolicUtils.unwrap(Num(p[k])) for k in 1:6]
            pairwise(terms) = foldl(
                (c, a) -> add_cnum(c, to_cnum(a)), terms; init = CNUM_ZERO,
            )

            # distinct monomials: a 6-term Poly, matching the pairwise fold exactly
            distinct = [Num(k) * p[k] for k in 1:6]
            cd = to_cnum(sum(distinct))
            @test is_poly(cd) && length(cd.tail.terms) == 6
            @test isequal(cd, pairwise(distinct))

            # repeated factor sets coalesce: 1·g + 2·g + 3·g -> 6g (one term)
            cc = to_cnum(g + 2g + 3g)
            @test is_poly(cc) && length(cc.tail.terms) == 1
            @test isequal(to_num(cc), Complex(Num(6g), Num(0)))

            # native + poly mixed in one sum: constants collapse to a single term
            cm = to_cnum(2 + g + 3 + κ)
            @test isequal(to_num(cm), Complex(Num(5 + g + κ), Num(0)))

            # full cancellation across many terms -> exact zero
            @test iszero_cnum(to_cnum(g + κ - g - κ))

            # a sym-leaf term (irreducible call) drops to the fallback add and the
            # polynomial part still coalesces around it
            cs = to_cnum(g + sin(ω) + g)
            @test isequal(to_num(cs), Complex(Num(2g + sin(ω)), Num(0)))
        end

        @testset "wide products fold in one pass (batched rec_prod)" begin
            # A `*`-headed coefficient of single-monomial factors collapses to one
            # monomial; the factor lists merge in a single pass, not pairwise.
            # Result must be byte-identical to the pairwise `mul_cnum` fold.
            @variables r[1:6]
            pairwise(fs) = foldl(
                (c, a) -> mul_cnum(c, to_cnum(a)), fs; init = CNUM_ONE,
            )

            # distinct atoms: one monomial with 6 factors, matching the pairwise fold
            facs = [Num(r[k]) for k in 1:6]
            cp = to_cnum(prod(facs))
            @test is_poly(cp) && length(cp.tail.terms) == 1
            @test length(cp.tail.terms[1].syms) == 6
            @test isequal(cp, pairwise(facs))

            # scalars multiply, exponents accumulate across repeated factors
            ce = to_cnum(2g * 3g * κ)            # 6 g^2 κ
            @test isequal(to_num(ce), Complex(Num(6 * g^2 * κ), Num(0)))

            # a single-monomial part times a multi-term factor distributes (intrinsic)
            cdist = to_cnum(g * κ * (g + κ))
            @test isequal(to_num(cdist), Complex(Num(g^2 * κ + g * κ^2), Num(0)))

            # factor cancellation drops the atom entirely
            @test isequal(to_num(to_cnum(g * κ / g)), Complex(Num(κ), Num(0)))

            # a zero factor collapses the whole product to native zero
            @test iszero_cnum(to_cnum(g * κ * 0 * ω))

            # a sym-leaf factor drops to the fallback multiply; atoms still merge
            csl = to_cnum(g * exp(κ) * g)
            @test isequal(to_num(csl), Complex(Num(g^2 * exp(κ)), Num(0)))
        end

        @testset "(sum)^n is stored in canonical expanded form" begin
            # the package-wide 'always expand' invariant extends to polynomial coefficients
            c = to_cnum((g + κ)^2)
            @test is_poly(c)
            @test isequal(to_num(c), Complex(Num(g^2 + 2 * g * κ + κ^2), Num(0)))
        end

        @testset "conjugation" begin
            # real parameters: conj is identity
            @test isequal(to_num(conj_cnum(to_cnum(ω * g))), Complex(Num(ω * g), Num(0)))
            # complex-symtype parameter: conj reaches the factor
            @test isequal(to_num(conj_cnum(to_cnum(gc))), Complex(Num(conj(gc)), Num(0)))
            @test isequal(to_num(conj_cnum(to_cnum(im * g))), Complex(Num(0), Num(-g)))
        end

        @testset "conjugation involution folds (regression #7cc3ad7)" begin
            # `conj` is an involution: conjugating twice must return the original
            # factor, not nest a `conj(conj(x))` that never folds and survives
            # downstream.
            c = to_cnum(gc)
            @test isequal(conj_cnum(conj_cnum(c)), c)
            @test isequal(to_num(conj_cnum(conj_cnum(c))), Complex(Num(gc), Num(0)))
            @test hash(conj_cnum(conj_cnum(c))) == hash(c)
            # double-conjugating a scaled complex factor also returns the original
            cs = to_cnum(g * gc)
            @test isequal(conj_cnum(conj_cnum(cs)), cs)
        end

        @testset "complex parameters and irreducible couplings stay native" begin
            # `@variables _::Complex` is stored as `Complex(real(_), imag(_))`; its
            # real/imag parts are recognized atoms, so the parameter stays on the Poly
            # path instead of poisoning downstream arithmetic with a Complex{Num} tail.
            @variables gv::Complex γ::Real
            @test is_poly(to_cnum(gv))
            @test is_poly(to_cnum(√(γ) * gv))
            # conjugation is native and faithful: real γ is self-conjugate, the complex
            # parameter's imaginary part flips sign.
            @test isequal(to_num(conj_cnum(to_cnum(gv))), conj(gv))
            @test isequal(to_num(conj_cnum(to_cnum(√(γ) * gv))), conj(√(γ) * gv))
        end

        @testset "radicals canonicalize via rational exponents" begin
            @variables gv::Complex γ::Real
            # squaring a radical folds to its radicand, matching Symbolics' `sqrt(g)^2 -> g`
            @test isequal(mul_cnum(to_cnum(sqrt(g)), to_cnum(sqrt(g))), to_cnum(g))
            @test isequal(to_cnum(sqrt(g)^2), to_cnum(g))
            @test hash(to_cnum(sqrt(g)^2)) == hash(to_cnum(g))
            three = mul_cnum(mul_cnum(to_cnum(sqrt(g)), to_cnum(sqrt(g))), to_cnum(sqrt(g)))
            @test isequal(three, to_cnum(g * sqrt(g)))
            @test isequal(three, to_cnum(g^(3 // 2)))
            @test isequal(mul_cnum(to_cnum(sqrt(γ)), to_cnum(γ)), to_cnum(γ^(3 // 2)))
            @test isequal(to_cnum(sqrt(sqrt(g))), to_cnum(g^(1 // 4)))
            # sqrt(g) and g^(1//2) are the same coefficient
            @test isequal(to_cnum(sqrt(g)), to_cnum(g^(1 // 2)))
            for c in (to_cnum(sqrt(g)), to_cnum(γ^(3 // 2)), to_cnum(g^(1 // 4)))
                @test isequal(to_cnum(to_num(c)), c)   # materialization round-trips
            end
            @test isequal(to_num(conj_cnum(to_cnum(√(γ) * gv))), conj(√(γ) * gv))
            # radical of a product / scaled atom is not distributed: symbolic leaf
            @test is_symbolic_cnum(to_cnum(sqrt(g * κ)))
            @test is_symbolic_cnum(to_cnum((2g)^(1 // 2)))
            for e in (-1 // 2, -3 // 2)   # negative fractional exps use the divide branch
                c = to_cnum(g^e)
                @test is_poly(c) && isequal(to_cnum(to_num(c)), c)
            end
            @test isequal(mul_cnum(to_cnum(sqrt(g)), to_cnum(1 / g)), to_cnum(g^(-1 // 2)))
            @test is_native(to_cnum(Num(4)^(1 // 2)))   # numeric base stays native
            @test to_num(to_cnum(Num(4)^(1 // 2))) == 2
            # only exact `Rational` exponents fold; a float exponent (even `1/3`) stays
            # a symbolic leaf, so the `^` recognizer needs no float-to-rational guesswork
            @test is_symbolic_cnum(to_cnum(g^(1 / 3)))
            @test is_symbolic_cnum(to_cnum(g^0.5))
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
            @test isequal(to_num(neg_cnum(to_cnum(2 * ω * g))), Complex(Num(-2 * ω * g), Num(0)))
            @test iszero_cnum(add_cnum(to_cnum(ω * g), neg_cnum(to_cnum(ω * g))))
            # Raw symbolic terms do not go through a CAS during addition, but exact opposites
            # must still cancel for scalar identities (e.g. symbolic orthogonal matrices).
            @variables u v
            raw = to_cnum(cos(u * v))
            @test iszero_cnum(add_cnum(raw, neg_cnum(raw)))
        end

        @testset "array-indexed parameters are recognized" begin
            @variables ωs[1:3] gs[1:3]
            @test is_poly(to_cnum(ωs[1]))
            @test is_poly(to_cnum(ωs[1] * gs[2]))
            @test isequal(to_num(to_cnum(ωs[1] * gs[2])), Complex(Num(ωs[1] * gs[2]), Num(0)))
            # products merge array-indexed factors like any other atom
            p = mul_cnum(to_cnum(ωs[1] * gs[2]), to_cnum(ωs[1]))   # -> ωs[1]^2 gs[2]
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
            c = to_cnum(Complex{Num}(exp(1.0im * ω * t)))
            @test !occursin("exp(0)", repr(to_num(c)))
            @test isequal(to_num(c), Complex(Num(cos(t * ω)), Num(sin(t * ω))))
            # exact identities at 0 fold to native
            @test is_native(to_cnum(Num(exp(Num(0)))))
            @test to_complex(to_cnum(Num(exp(Num(0))))) ≈ 1
            @test to_complex(to_cnum(Num(cos(Num(0))))) ≈ 1
            @test to_complex(to_cnum(Num(sin(Num(0))))) ≈ 0
            # symbol / non-zero-constant / irrational args stay symbolic (exactness gate)
            @test !is_native(to_cnum(Num(exp(ω))))
            @test !is_native(to_cnum(Num(exp(Num(2)))))
            @test !is_native(to_cnum(Num(exp(Num(0.1)))))
            @test !is_native(to_cnum(Num(cos(Num(1)))))
            @test !is_native(to_cnum(cos(Num(π))))
            @test !is_native(to_cnum(sin(Num(π))))
        end

        @testset "explicit complex slots differentiate component-wise" begin
            @variables z::Number
            c = to_cnum(Complex(Num(z), Num(z)))
            expected = to_cnum(1 + im)
            @test isequal(Symbolics.derivative(c, z), expected)
            D = Symbolics.Differential(z)
            @test isequal(to_cnum(Symbolics.expand_derivatives(D(c))), expected)
        end
    end

    @testset "unit phases" begin
        @variables ω J t θ
        p = phase_coeff(ω * t)
        @test is_poly(p)
        @test is_phase(p.tail.terms[1].syms[1])
        # conjugation is exponent negation, so a phase cancels against its own conjugate
        @test mul_cnum(p, conj_cnum(p)) == CNUM_ONE
        pp = mul_cnum(p, p)
        @test mul_cnum(pp, conj_cnum(pp)) == CNUM_ONE
        # `expim(0)` is 1, not an atom
        @test phase_coeff(0) == CNUM_ONE
        @test to_cnum(expim(Num(0))) == CNUM_ONE
        # one atom per frequency however the argument is spelled, and whatever its sign
        @test isequal(phase_coeff((ω + 2J) * t), phase_coeff(ω * t + 2 * J * t))
        @test mul_cnum(phase_coeff(-ω * t), phase_coeff(ω * t)) == CNUM_ONE
        # commensurate frequencies are exponent arithmetic, not an angle-addition identity
        @variables E1 E2 E3
        ph(i, j) = mul_cnum(phase_coeff(i * t), conj_cnum(phase_coeff(j * t)))
        @test iszero_cnum(
            add_cnum(mul_cnum(ph(E1, E2), ph(E2, E3)), neg_cnum(ph(E1, E3)))
        )
        # lowering round-trips back onto the same atom
        @test isequal(p, to_cnum(to_num(p)))
        @test isequal(conj_cnum(p), to_cnum(to_num(conj_cnum(p))))
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
        @test px * inv(px) == CNUM_ONE
        @test expim(x) * expim(y - x) == expim(y)
        @test expim((ω + 2g) * t) * expim(-2g * t) == expim(ω * t)

        combined = g * px * py
        combined_term = only((combined.tail::Poly).terms)
        @test count(is_phase, combined_term.syms) == 1
        @test only(combined_term.exps[findall(is_phase, combined_term.syms)]) in (-1, 1)

        hs = FockSpace(:phase_collection)
        a = Destroy(hs, :a)
        @test iszero(px * py * a - expim(x + y) * a)

        @test_throws MethodError px^(1 // 2)
        @test_throws MethodError sqrt(px)
        phase_symbol = only(only((px.tail::Poly).terms).syms)
        phase_root = SymbolicUtils.term(sqrt, phase_symbol; type = Complex{Real})
        @test_throws ArgumentError to_cnum(phase_root)

        cg, ch = to_cnum(g), to_cnum(h)
        @test phase_allocations(() -> cg * ch) <= 10
        @test phase_allocations(() -> px * conj(px)) <= 18
        @test phase_allocations(() -> px * px) <= 165
        @test phase_allocations(() -> px * py) <= 285
    end

    @testset "phase substitution, calculus, and projections" begin
        @variables ω t θ g z::Complex
        p = expim(ω * t)

        @test substitute(p, Dict()) === p
        @test (@inferred substitute(p, Dict(ω => 2ω))) == expim(2ω * t)
        @test substitute(p, Dict(t => 0)) == CNUM_ONE
        @test substitute(p, Dict(g => 2g)) === p
        @test substitute(p, Dict(g => z)) === p
        @test to_complex(substitute(p, Dict(ω => 2.0, t => 0.25))) ≈ exp(0.5im)
        @test_throws ArgumentError substitute(p, Dict(ω => z))
        @test_throws ArgumentError substitute(p, Dict(ω * t => z))
        @test_throws ArgumentError substitute(p, Dict(ω => 1 + 2im))
        @test_throws ArgumentError substitute(p, Dict(ω => 1 + 0im))

        D = Symbolics.Differential(t)
        differentiated = @inferred Symbolics.derivative(p, t)
        @test differentiated == im * ω * p
        @test to_cnum(Symbolics.expand_derivatives(D(p))) == differentiated

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
        @test_throws MethodError abs(to_cnum(g))
        @test_throws MethodError abs2(to_cnum(g))
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
        @test all(
            term -> isequal(to_num(term.amplitude), Complex(Num(1 // 2), Num(0))),
            phase_terms(cosine_phase),
        )
        @test trigonometric_form(cosine_phase) == to_cnum(cos(θ))
        @test trigonometric_form(sine_phase) == to_cnum(sin(θ))
        @test trigonometric_form(composite_phase) == to_cnum(cos(ω * t) + sin(ω * t))
        @test trigonometric_form(expim(θ)^2) == to_cnum(cos(2θ) + im * sin(2θ))
        phase_over_root = cos(ω * t) * expim(-ω * t) / sqrt(2 * ω)
        exponential_phase_over_root = @inferred exponential_form(phase_over_root)
        @test !occursin("cos(", string(exponential_phase_over_root))
        expected_phase_over_root = (1 + expim(-2 * ω * t)) / (2 * sqrt(2 * ω))
        @test iszero(simplify(exponential_phase_over_root - expected_phase_over_root))
        @test (@inferred exponential_form(exp(im * θ))) == expim(θ)
        @test (@inferred exponential_form(cis(θ))) == expim(θ)
        @test (@inferred exponential_form(intact_exp)) == expim(θ)
        @test (@inferred exponential_form(intact_cis)) == expim(θ)
        @test isequal(exponential_form(arbitrary_exp), to_cnum(arbitrary_exp))
        @test exponential_form(exp(θ)) == to_cnum(exp(θ))

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

    @testset "finite phase decomposition" begin
        @variables θ ω t g κ

        reconstruct(c) = sum(
            term.amplitude * expim(term.phase) for term in phase_terms(c);
            init = CNUM_ZERO,
        )

        coefficients = (
            to_cnum(3 // 4),
            2 * expim(-ω * t) + 5 * expim(2 * ω * t),
            exponential_form(cos(ω * t + θ)),
            (g + g * expim(-2 * ω * t)) / sqrt(2 * κ),
            (1 + expim(-ω * t)) * (1 + expim(2 * ω * t)),
        )
        for coefficient in coefficients
            terms = @inferred phase_terms(coefficient)
            @test iszero(simplify(reconstruct(coefficient) - exponential_form(coefficient)))
            @test all(term -> term.amplitude isa Coeff && term.phase isa Num, terms)
        end

        @test isempty(@inferred phase_terms(CNUM_ZERO))
        @test_throws ArgumentError phase_terms(CNUM_ONE / (1 + expim(ω * t)))
        phase_symbol = only(only((expim(θ).tail::Poly).terms).syms)
        nested_phase = SymbolicUtils.term(exp, phase_symbol; type = Complex{Real})
        @test_throws ArgumentError phase_terms(to_cnum(nested_phase))
    end

    # The public `expim` used to return a `Num`. `Base.conj(::Num)` is the identity and
    # `Complex * Num` splits into real/imag halves, so every law below silently failed.
    @testset "public expim is a coefficient" begin
        @variables ω t
        h = FockSpace(:f)
        a = Destroy(h, :a)
        p = expim(ω * t)
        @test p isa Coeff

        @test conj(p) * p == CNUM_ONE
        @test (@inferred inv(p)) == conj(p)
        @test (@inferred p^(-1)) == conj(p)
        @test !isequal(conj(p), p)
        @test (im * p) * conj(p) == CNUM_IM
        @test is_poly(im * p)
        @test isequal(expim(ω * t) * expim(-ω * t) * a, 1 * a)
        @test isequal(p * a, expim(ω * t) * a)

        # a phase over a literal is its value, so numerically cancelling terms fold
        @test is_native(expim(Num(1.0)))
        @test iszero(expim(Num(1.0)) * a - exp(im * 1.0) * a)

        # `type`/`shape` take part in hash-consing, so an unregistered head would come back
        # from any rebuild as a different atom whose `conj` is the identity
        b = SecondQuantizedAlgebra.expim_symbolic(SymbolicUtils.unwrap(ω * t))
        rebuilt = SymbolicUtils.maketerm(
            typeof(b), SymbolicUtils.operation(b), SymbolicUtils.arguments(b),
            SymbolicUtils.metadata(b),
        )
        @test rebuilt === b
        @test SymbolicUtils.symtype(b) === Complex{Real}
        @test Symbolics.substitute(b, Dict(SymbolicUtils.unwrap(ω) => 2.0)) ===
            SecondQuantizedAlgebra.expim_symbolic(SymbolicUtils.unwrap(2.0 * t))

        # conjugation negates the argument rather than conjugating it
        @test isequal(
            SecondQuantizedAlgebra.qadjoint(Num(b)),
            SecondQuantizedAlgebra.expim_symbolic(SymbolicUtils.unwrap(-ω * t)),
        )

        # a registered derivative, not an inert `Differential` node
        d = Symbolics.expand_derivatives(Symbolics.Differential(t)(Num(b)))
        @test isequal(to_cnum(d), mul_cnum(mul_cnum(CNUM_IM, to_cnum(ω)), p))
    end

    @testset "raw symbolic coefficient fallbacks" begin
        @variables θ φ ψ ω t g κ r α::Number
        phase = expim(θ)
        raw_power = phase^2 * exp(ω * t)
        raw_product = phase * expim(φ) * exp(ω * t)
        raw_denominator = phase / (expim(φ) * exp(ω * t))

        @test isequal(phase_terms(raw_power)[1].phase, 2θ)
        @test isequal(phase_terms(raw_product)[1].phase, θ + φ)
        @test isequal(phase_terms(raw_denominator)[1].phase, θ - φ)
        @test occursin("cos", string(trigonometric_form(raw_product)))
        @test occursin("sin", string(trigonometric_form(raw_product)))
        @test !occursin("expim", string(to_num(exp(ω * t) / (phase * expim(φ)))))
        @test iszero(raw_product + (-raw_product))
        @test !isequal(conj(raw_product), raw_product)

        # The raw phase atom is constructed only to exercise the boundary normalizer; the
        # behavior is checked through the public `to_num` representation.
        raw_phase = SymbolicUtils.term(
            expim, SymbolicUtils.unwrap(θ); type = Complex{Real}, shape = UnitRange{Int}[],
        )
        raw_phase_atom(x) = SymbolicUtils.term(
            expim, SymbolicUtils.unwrap(x); type = Complex{Real}, shape = UnitRange{Int}[],
        )
        raw_product_tree = SymbolicUtils.term(
            *, raw_phase, raw_phase_atom(φ), SymbolicUtils.unwrap(exp(ω * t));
            type = Complex{Real}, shape = UnitRange{Int}[],
        )
        raw_denominator_tree = SymbolicUtils.term(
            /,
            raw_phase,
            SymbolicUtils.term(
                *, raw_phase_atom(φ), raw_phase_atom(ψ);
                type = Complex{Real}, shape = UnitRange{Int}[],
            );
            type = Complex{Real}, shape = UnitRange{Int}[],
        )
        raw_nonphase_power_tree = SymbolicUtils.term(
            ^, SymbolicUtils.unwrap(r), -1; type = Complex{Real}, shape = UnitRange{Int}[],
        )
        raw_conjugate_tree = SymbolicUtils.term(
            conj, raw_phase; type = Complex{Real}, shape = UnitRange{Int}[],
        )
        @test isequal(
            phase_terms(exponential_form(Num(raw_phase^2)))[1].phase, 2θ,
        )
        @test isequal(
            phase_terms(exponential_form(Num(raw_product_tree)))[1].phase, θ + φ,
        )
        @test isequal(
            phase_terms(exponential_form(Num(raw_denominator_tree)))[1].phase,
            θ - φ - ψ,
        )
        @test length(phase_terms(exponential_form(Num(raw_nonphase_power_tree)))) == 1
        @test isequal(
            phase_terms(exponential_form(Num(raw_conjugate_tree)))[1].phase, -θ,
        )
        @test isequal(
            to_num(from_raw(raw_phase^2; normalize = false)),
            Complex(cos(2θ), sin(2θ)),
        )
        @test isequal(
            to_num(from_raw(conj(raw_phase); normalize = false)),
            Complex(cos(-θ), sin(-θ)),
        )
        @test isequal(
            to_num(from_raw(real(raw_phase); normalize = false)),
            Complex(cos(θ), Num(0)),
        )
        @test isequal(
            to_num(from_raw(imag(raw_phase); normalize = false)),
            Complex(sin(θ), Num(0)),
        )
        @test to_num(from_raw(abs(raw_phase); normalize = false)) == 1
        @test to_num(from_raw(abs2(raw_phase); normalize = false)) == 1

        large_q = complex(big(2)^70, big(3)^70) * Destroy(FockSpace(:raw), :a)
        large = prefactor(large_q)
        @test isequal(real(large), Num(big(2)^70))
        @test isequal(imag(large), Num(big(3)^70))

        exact_slots = cnum_sym(Num(1 // 2), Num(1 // 3))
        @test isequal(to_num(exact_slots), Complex(Num(1 // 2), Num(1 // 3)))
        raw_slots = cnum_sym(g, κ)
        @test isequal(real(to_num(raw_slots)), g)
        @test isequal(imag(to_num(raw_slots)), κ)
        numeric_slots = cnum_sym(Num(0.5), Num(0.25))
        @test isequal(to_num(numeric_slots), Complex(Num(0.5), Num(0.25)))

        raw_complex_expr = SymbolicUtils.term(
            complex,
            SymbolicUtils.unwrap(g),
            SymbolicUtils.unwrap(κ);
            type = Complex{Real},
            shape = UnitRange{Int}[],
        )
        raw_complex = from_raw(raw_complex_expr; normalize = false)
        raw_division = from_raw(
            raw_complex_expr / SymbolicUtils.unwrap(r); normalize = false,
        )
        @test isequal(real(conj(raw_complex)), g)
        @test isequal(imag(conj(raw_complex)), -κ)
        @test isequal(raw_parts(raw_complex), (g, κ))
        @test !isequal(conj(raw_division), raw_division)
        renamed = substitute(raw_complex, Dict(SymbolicUtils.unwrap(g) => 2))
        @test !isequal(renamed, raw_complex)

        raw_number_product = exponential_form(
            Num(
                SymbolicUtils.term(
                    *, raw_phase, SymbolicUtils.unwrap(α);
                    type = Complex{Real}, shape = UnitRange{Int}[],
                ),
            ),
        )
        @test !isequal(conj(raw_number_product), raw_number_product)

        unit_phase = Complex(3 // 5, 4 // 5) * phase
        @test isequal(real(unit_phase), (3 // 5) * cos(θ) - 0.8sin(θ))
        @test inv(exponential_form((1 // 2) * g)) isa Coeff
        @test inv(exponential_form(cos(θ) + sin(θ))) isa Coeff
        mixed_scalars = exponential_form(0.7cos(θ)) + exponential_form((1 // 2) * cos(θ))
        @test isequal(mixed_scalars, exponential_form(1.2cos(θ)))

        @test length(phase_terms(phase^(-1) * exp(ω * t))) == 1
        @test length(phase_terms(exponential_form(exp(ω * t)))) == 1
        @test length(phase_terms(exponential_form(exp(ω * t)^(1 // 2)))) == 1
        @test_throws ArgumentError substitute(raw_product, Dict(SymbolicUtils.unwrap(θ) => 1im))
    end
end
