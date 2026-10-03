using SecondQuantizedAlgebra
using Symbolics: Symbolics, @variables, Num
using SymbolicUtils: SymbolicUtils
using Test
using LinearAlgebra: I
import SecondQuantizedAlgebra: Coeff, to_cnum, Monomial, ExactComplex, BigExactComplex

@testset "Symbolic coefficients" begin
    h = FockSpace(:coefficients)
    a = Destroy(h, :a)

    coefficient(x) = get_prefactor(x * a)
    stored_coefficient(x) = only(x).second

    @testset "numeric and exact coefficients" begin
        @test coefficient(2) == 2
        @test coefficient(2 + 3im) == 2 + 3im
        @test isequal(coefficient(1 // 4), Complex(Num(1 // 4), Num(0)))
        @test isequal(
            coefficient(Complex(1 // 4, 1 // 2)),
            Complex(Num(1 // 4), Num(1 // 2))
        )
        @test string(real(coefficient(2))) == "2"
        @test occursin("1//4", string(real(coefficient(1 // 4))))

        large = Complex((2^53 + 1) // 1, 0 // 1)
        @test isequal(real(coefficient(large)), Num(2^53 + 1))
        @test isequal(imag(coefficient(large)), Num(0))
        @test occursin("1//4", string(real(coefficient((1 // 4) * a))))
        @test string(real(coefficient(0.25 * a))) == "0.25"
        @test iszero(simplify((1 // 4) * a + (1 // 4) * a - (1 // 2) * a))
        @test iszero(simplify((1 // 3) * a - (1 // 3) * a))
    end

    @testset "large and rational exactness boundaries" begin
        @variables r
        @test occursin("1//2", string(real(coefficient(r / 2 * a))))
        @test occursin("1//4", string(real(coefficient((1 // 4) * r * a))))
        large2 = big(2)^70
        @test coefficient(large2 * a) == large2
        @test isequal(real(coefficient((large2 + 1) * a)), Num(large2 + 1))
        @test string(real(coefficient(2))) == "2"
        @test string(real(coefficient(0.7))) == "0.7"
        @test isequal(coefficient(Complex(1 // 3, 0)), Complex(Num(1 // 3), Num(0)))
        @test iszero(simplify(((1 // 3) * a) * 3 - a))
    end

    @testset "conjugation and hash stability are observable" begin
        c = coefficient(2)
        @test isequal(conj(c), c)
        @test hash(conj(c)) == hash(c)
        doubled = conj(conj(coefficient(2 + 3im)))
        @test isequal(doubled, coefficient(2 + 3im))
        @test hash(doubled) == hash(coefficient(2 + 3im))
        @test isequal(conj(coefficient(2 + 3im)), coefficient(2 - 3im))
        q = 2 * a' * a
        @test isequal(q', q)
        @test hash(q') == hash(q)
        @variables g k
        simplified = simplify(((g + g * k) / (1 + k)) * a)
        @test isequal(simplified, g * a)
        @test hash(simplified) == hash(g * a)
    end

    @testset "wide sums fold and dedup through public algebra" begin
        @variables p[1:6] ω g κ
        distinct = [Num(k) * p[k] for k in 1:6]
        s = sum(distinct[k] * a for k in 1:6)
        @test length(s) == 1
        @test occursin("p[1]", string(get_prefactor(s)))
        @test occursin("p[6]", string(get_prefactor(s)))
        coalesced = g * a + 2g * a + 3g * a
        @test iszero(simplify(coalesced - 6g * a))
        @test length(coalesced) == 1
        different = g * a + κ * a
        @test length(different) == 1
        @test occursin("g", string(get_prefactor(different)))
        @test occursin("κ", string(get_prefactor(different)))
        @test iszero(simplify(different - (g + κ) * a))
        product = (2g * a) * (3κ * a')
        @test iszero(simplify(product - 6 * g * κ * a * a'))
        cancelling = g * a + κ * a - g * a - κ * a
        @test iszero(cancelling)
        @variables q
        two_ops = p[1] * a + p[2] * a'
        @test length(two_ops) == 2
    end

    @testset "radicals and exact powers stay canonical" begin
        @variables ω g κ
        @test iszero(simplify(sqrt(g) * sqrt(g) - g))
        @test iszero(simplify(g^(1 // 2) * g^(1 // 2) - g))
        @test iszero(simplify(sqrt(4) * a - 2 * a))
        @test iszero(simplify(g^2 * a - g * g * a))
    end

    @testset "radicals of exact numbers stay exact" begin
        function contains_float(x)
            x = SymbolicUtils.unwrap(x)
            x isa AbstractFloat && return true
            x isa Complex && return contains_float(real(x)) || contains_float(imag(x))
            x isa Number && return false
            SymbolicUtils.isconst(x) && return contains_float(x.val)
            SymbolicUtils.iscall(x) || return false
            return any(contains_float, SymbolicUtils.arguments(x))
        end
        exact(x) = !contains_float(coefficient(x))

        root_two = sqrt(Num(2))
        @test exact(root_two)
        @test isequal(real(coefficient(root_two)), root_two)
        @test exact(cbrt(Num(2)))

        @test exact(get_prefactor(simplify(sqrt(Num(1 // 2)) * a)))
        @test iszero(simplify(root_two * root_two * a - 2 * a))

        @test coefficient(sqrt(Num(4))) == 2
        @test isequal(coefficient(sqrt(Num(1 // 4))), Complex(Num(1 // 2), Num(0)))
        @test coefficient(sqrt(Num(2.0))) ≈ sqrt(2.0)
    end

    @testset "exact coefficients beyond 2^53 and Int64" begin
        b = Destroy(FockSpace(:wide), :b)
        exact_parts(c) = (Symbolics.value(real(c)), Symbolics.value(imag(c)))

        re, _ = exact_parts(get_prefactor((3^17 * a) * 3^17))
        @test re isa Integer && re == 3^34
        re, ip = exact_parts(
            get_prefactor(((2^27 + 1 + 2^27 * im) * a) * (2^27 + 1 + (2^27 + 2) * im)),
        )
        @test re == 1 && ip isa Integer && ip == 2 * (2^27 + 1)^2

        re, _ = exact_parts(get_prefactor(2 * ((1 // big(3)^25) * b * a)))
        @test re isa Rational && re == 2 // big(3)^25
        re, _ = exact_parts(get_prefactor((1 // big(2)^40) * b * a))
        @test re isa Rational && re == 1 // 2^40
        re, _ = exact_parts(
            get_prefactor((1 // Int128(3)^25) * b * a + (1 // Int128(2)^40) * b * a),
        )
        @test re isa Rational && re == 1 // big(3)^25 + 1 // big(2)^40
        @test isequal(((1 // big(3)^25) * b * a) * big(3)^25, b * a)

        wide = (1 // 3^39) * a + (1 // 2^40) * a
        back = wide - (1 // 2^40) * a
        @test isequal(back, (1 // 3^39) * a)
        @test hash(back) == hash((1 // 3^39) * a)
        # A constant that fits the small tier returns to the native tier.
        @test SecondQuantizedAlgebra.is_native(stored_coefficient(back))
        # Conjugating a `typemin(Int)` imaginary part overflows `Rational{Int}`.
        edge = to_cnum(Complex(1 // 3, typemin(Int) // 1))
        @test isequal(conj(edge), to_cnum(Complex(big(1) // 3, -big(typemin(Int)) // 1)))

        # Native Gaussian integers overflow into the exact tiers instead of wrapping.
        top = typemax(Int)
        @test isequal(to_cnum(top) * 2, to_cnum(2 * big(top)))
        @test isequal(to_cnum(top) + 1, to_cnum(big(top) + 1))
        @test isequal(-to_cnum(typemin(Int)), to_cnum(-big(typemin(Int))))
        @test isequal(conj(to_cnum(Complex(0, typemin(Int)))), to_cnum(Complex(0, -big(typemin(Int)))))
        @test isequal(to_cnum(top) / to_cnum(3 + im), to_cnum(Complex(3 * big(top) // 10, -big(top) // 10)))
        # A big-tier polynomial combined with a native integer.
        @variables g
        big_poly = to_cnum(g) * to_cnum(1 // big(3)^40)
        @test isequal((big_poly + 5) - big_poly, to_cnum(5))
        @test isequal(big_poly / 5, to_cnum(g) * to_cnum(1 // (5 * big(3)^40)))

        @test isequal(to_cnum(1) / to_cnum(3), to_cnum(1 // 3))
        @test isequal(inv(to_cnum(3 + 4im)), to_cnum(3 // 25 - 4 // 25 * im))
        @test isequal(to_cnum(1.0) / to_cnum(0.3), to_cnum(1.0 / 0.3))
        @test isequal(to_cnum(5) / to_cnum(1 // 3), to_cnum(15))

        @variables θ
        phase = SecondQuantizedAlgebra.expim(-θ) * to_cnum(Complex(3 // 5, 4 // 5))
        @test isequal(real(phase), (3 // 5) * cos(θ) + (4 // 5) * sin(θ))
        @test isequal(imag(phase), (4 // 5) * cos(θ) - (3 // 5) * sin(θ))
    end

    @testset "numeric radicals multiply exactly" begin
        inverse_root_six = 1 / sqrt(Num(6))
        weighted = inverse_root_six * a
        @test isequal(stored_coefficient(commutator(weighted, weighted')), to_cnum(1 // 6))

        @test isequal(to_cnum(sqrt(Num(2))) * to_cnum(sqrt(Num(3))), to_cnum(sqrt(Num(6))))
        @test hash(to_cnum(sqrt(Num(2))) * to_cnum(sqrt(Num(3)))) == hash(to_cnum(sqrt(Num(6))))
        @test isequal(to_cnum(sqrt(Num(8))), 2 * to_cnum(sqrt(Num(2))))
        @test isequal(to_cnum(1 / sqrt(Num(2))), to_cnum(sqrt(Num(2))) / 2)
        @test isequal(to_cnum(cbrt(Num(2)))^3, to_cnum(2))

        @variables g
        raw = to_cnum(sin(g + 1))
        half = to_cnum(1 / sqrt(Num(2)))
        @test isequal((half * raw) * half, to_cnum(1 // 2) * raw)

        # A prime above the trial bound stays an unevaluated radical.
        large_prime = 4611686018427387847
        @test SecondQuantizedAlgebra.to_complex(to_cnum(sqrt(Num(large_prime)))) ≈ sqrt(large_prime)
    end

    @testset "radicals reduce to prime atoms and compose with the big tier" begin
        via_product = to_cnum(sqrt(Num(2))) * to_cnum(sqrt(Num(3)))
        via_radicand = to_cnum(sqrt(Num(6)))
        via_folded = to_cnum(sqrt(Num(24))) / to_cnum(2)
        @test isequal(via_folded, via_radicand)
        @test hash(via_folded) == hash(via_radicand)
        via_division = to_cnum(sqrt(Num(12))) / to_cnum(2)
        @test isequal(via_division, to_cnum(sqrt(Num(3))))
        @test hash(via_division) == hash(to_cnum(sqrt(Num(3))))

        small = to_cnum(sqrt(Num(2))) * to_cnum(2)^40
        @test isequal(small^2, to_cnum(2)^81)
        big_tier = to_cnum(sqrt(Num(2))) * to_cnum(2)^63
        @test isequal(big_tier^2, to_cnum(2)^127)

        @test isequal(to_cnum(sqrt(Num(big(2)^71))), to_cnum(big(2))^35 * to_cnum(sqrt(Num(2))))

        @test isequal(inv(to_cnum(sqrt(Num(2)))), to_cnum(sqrt(Num(2))) / to_cnum(2))

        @test isequal(to_cnum(sqrt(Num(2)))^200, to_cnum(2)^100)
    end

    @testset "cross-tier stress: prime-radical sums forced through the big tier" begin
        primes6 = [2, 3, 5, 7, 11, 13]
        s = sum(to_cnum(sqrt(Num(p))) for p in primes6)
        s2 = s^2
        big_scale = to_cnum(2)^70
        scaled = s2 * big_scale

        # One constant term plus the 15 cross terms `√p·√q`, all in the big tier.
        @test scaled.tail.terms isa Vector{Monomial{BigExactComplex}}
        @test length(scaled.tail.terms) == 16

        built = sum(
            to_cnum(2)^70 * (to_cnum(sqrt(Num(p))) * to_cnum(sqrt(Num(q))))
                for p in primes6, q in primes6
        )
        @test isequal(scaled, built)
        @test hash(scaled) == hash(built)

        back = scaled * inv(big_scale)
        @test back.tail.terms isa Vector{Monomial{ExactComplex}}
    end

    @testset "symbolic arithmetic stays faithful" begin
        @variables g κ r
        @test isequal(coefficient(g), Complex(Num(g), Num(0)))
        @test isequal(coefficient(im * g), Complex(Num(0), Num(g)))
        @test isequal(coefficient(g + im * κ), Complex(Num(g), Num(κ)))
        @test isequal(coefficient(g * κ), Complex(Num(g * κ), Num(0)))
        @test isequal(coefficient(g + g), Complex(Num(2g), Num(0)))
        @test isequal(coefficient(g * κ / g), Complex(Num(κ), Num(0)))
        @test iszero(simplify((g - g) * a))
        @test isequal(conj(coefficient(2 + 3im)), coefficient(2 - 3im))
        @test isequal(conj(conj(coefficient(g + im * κ))), coefficient(g + im * κ))

        expression = (g + g * r) / (1 + r) * a
        @test iszero(simplify(expression - g * a))
    end

    @testset "public arithmetic collects and cancels coefficients" begin
        @variables γ δ β
        opposite = (γ / δ) * a + (-γ / δ) * a
        collected = (γ / δ) * a + (β / δ) * a

        @test iszero(simplify(opposite))
        @test iszero(simplify(collected - ((γ + β) / δ) * a))
        @test !iszero(simplify(collected))
    end

    @testset "additive and multiplicative identities" begin
        @variables x y

        @test one(Coeff) isa Coeff
        @test zero(Coeff) isa Coeff
        @test isone(one(Coeff))
        @test iszero(zero(Coeff))
        @test isequal(one(to_cnum(x)), one(Coeff))
        @test isequal(zero(to_cnum(x)), zero(Coeff))
        @test isequal(oneunit(Coeff), one(Coeff))

        c = to_cnum(x)
        @test isequal(c * one(Coeff), c)
        @test isequal(c + zero(Coeff), c)

        # `Coeff` is not a `Number`, so generic reductions need the identities above.
        @test isequal(sum(Coeff[to_cnum(x), to_cnum(y)]), to_cnum(x + y))
        @test isequal(prod(Coeff[to_cnum(x), to_cnum(y)]), to_cnum(x * y))
        @test isequal(sum(Coeff[]), zero(Coeff))
        @test isequal(prod(Coeff[]), one(Coeff))

        # Array constructors reach for the identities; `I` goes through the type itself.
        @test all(iszero, zeros(Coeff, 2, 2))
        @test all(isone, ones(Coeff, 2, 2))
        @test isequal(Coeff(2), to_cnum(2))
        @test isequal(
            Matrix{Coeff}(I, 2, 2),
            Coeff[one(Coeff) zero(Coeff); zero(Coeff) one(Coeff)],
        )
    end
end
