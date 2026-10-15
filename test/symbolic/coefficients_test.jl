using SecondQuantizedAlgebra
using Symbolics: Symbolics, @variables, Num
using SymbolicUtils: SymbolicUtils
using Test
using LinearAlgebra: I
import SecondQuantizedAlgebra: Coeff, CNum, to_cnum, to_num, Monomial, ExactComplex, BigExactComplex, expim

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

    @testset "positive numeric factors leave symbolic root branches intact" begin
        @variables ω ν Z::Number
        @test isequal(sqrt(ω / 2) * sqrt(2ω) * a, ω * a)
        @test isequal(sqrt(2ω) * a, sqrt(Num(2)) * sqrt(ω) * a)
        @test isequal(sqrt(4ω) * a, 2 * sqrt(ω) * a)
        @test iszero(a / sqrt(2ω) - sqrt(ω / 2) * a / ω)
        @test isequal(sqrt(2.0ω) * a, sqrt(2.0) * sqrt(ω) * a)
        @test isequal(sqrt(2Z) * a, sqrt(Num(2)) * sqrt(Z) * a)
        @test isequal(sqrt(Z) * sqrt(Z) * a, Z * a)
        @test !isequal(sqrt(ω * ν) * a, sqrt(ω) * sqrt(ν) * a)
        @test !isequal(sqrt(ω * Z) * a, sqrt(ω) * sqrt(Z) * a)
        @test !isequal(sqrt(-ω) * a, im * sqrt(ω) * a)
        @test !isequal(sqrt(Z^2) * a, Z * a)
        @test !isequal(sqrt(-Z) * a, im * sqrt(Z) * a)
        @test !isequal(sqrt(1 / Z) * a, a / sqrt(Z))
        @test !isequal((-ω)^(1 // 3) * a, -cbrt(ω) * a)
        @test iszero(sqrt(-Z) * a - sqrt(-Z) * a)

        value(expr, rules) = SecondQuantizedAlgebra.to_complex(stored_coefficient(substitute(expr, rules)))
        @test value(sqrt(ω * ν) * a, Dict(ω => -1, ν => -1)) ≈ 1
        square_root = Num(SymbolicUtils.Term{SymbolicUtils.SymReal}(sqrt, [Symbolics.unwrap(ω^2)]; type = Real))
        @test !isequal(square_root * a, ω * a)
        @test value(square_root * a, Dict(ω => -2)) ≈ 2
        @test value(cbrt(ω) * a, Dict(ω => -8)) ≈ -2
        @test value(sqrt(1 / Z) * a, Dict(Z => -2)) ≈ im / sqrt(2)
        for z in (-1 + 0.25im, -1 - 0.25im, 1 + 0.25im)
            rules = Dict(Z => z)
            @test value(sqrt(2Z) * a, rules) ≈ sqrt(2z)
            @test value(sqrt(Z^2) * a, rules) ≈ sqrt(z^2)
            @test value(sqrt(-Z) * a, rules) ≈ sqrt(-z)
            @test value(sqrt(1 / Z) * a, rules) ≈ sqrt(1 / z)
        end
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
        @test isequal(to_cnum(cbrt(Num(-8))), to_cnum(-2))
        @test isequal(to_cnum(sqrt(Num(-4))), to_cnum(2im))
        @test isequal(to_cnum(sqrt(Num(-2)))^2, to_cnum(-2))
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
        edge = to_cnum(Complex(1 // 3, typemin(Int) // 1))
        @test isequal(conj(edge), to_cnum(Complex(big(1) // 3, -big(typemin(Int)) // 1)))

        top = typemax(Int)
        @test isequal(to_cnum(top) * 2, to_cnum(2 * big(top)))
        @test isequal(to_cnum(top) + 1, to_cnum(big(top) + 1))
        @test isequal(-to_cnum(typemin(Int)), to_cnum(-big(typemin(Int))))
        @test isequal(conj(to_cnum(Complex(0, typemin(Int)))), to_cnum(Complex(0, -big(typemin(Int)))))
        @test isequal(to_cnum(top) / to_cnum(3 + im), to_cnum(Complex(3 * big(top) // 10, -big(top) // 10)))
        @variables g
        big_poly = to_cnum(g) * to_cnum(1 // big(3)^40)
        @test isequal((big_poly + 5) - big_poly, to_cnum(5))
        @test isequal(big_poly / 5, to_cnum(g) * to_cnum(1 // (5 * big(3)^40)))

        @test isequal(to_cnum(1) / to_cnum(3), to_cnum(1 // 3))
        @test isequal(inv(to_cnum(3 + 4im)), to_cnum(3 // 25 - 4 // 25 * im))
        re, _ = exact_parts(to_cnum(5) / to_cnum(1 // 3))
        @test re isa Integer && re == 15

        @variables θ
        phase = SecondQuantizedAlgebra.expim(-θ) * to_cnum(Complex(3 // 5, 4 // 5))
        @test isequal(real(phase), (3 // 5) * cos(θ) + (4 // 5) * sin(θ))
        @test isequal(imag(phase), (4 // 5) * cos(θ) - (3 // 5) * sin(θ))
    end

    @testset "floats stay floats and exact values stay exact" begin
        exact_parts(c) = (Symbolics.value(real(c)), Symbolics.value(imag(c)))
        @test exact_parts(to_cnum(2.0) * to_cnum(1 // 3))[1] isa AbstractFloat
        @test exact_parts(to_cnum(2) * to_cnum(1 // 3))[1] isa Rational
        for x in (Inf, -Inf, NaN, 1.0e-300, -0.0)
            @test isequal(SecondQuantizedAlgebra.to_complex(to_cnum(x)), ComplexF64(x) + 0.0)
        end
        @test isequal(to_cnum(1 // 0), to_cnum(Inf))
        @test isequal(to_cnum(-1 // 0), to_cnum(-Inf))
        @test isequal(ExactComplex(2, 1)^-2, ExactComplex(3 // 25, -4 // 25))
        @test isequal(ExactComplex(2, 0)^62, ExactComplex(2^62, 0))
        @test_throws OverflowError ExactComplex(2, 0)^63
    end

    @testset "float scalars of parameter polynomials cancel in both tiers" begin
        @variables g d
        for (x, y) in ((0.5g + d, -0.5g), (big(2)^70 * g + 0.5d, -0.5d), (-0.25g * d + 0.0g, 0.25g * d))
            sum_ = to_cnum(x) + to_cnum(y)
            expected = to_cnum(Symbolics.simplify(x + y))
            @test isequal(sum_, expected)
            @test hash(sum_) == hash(expected)
        end
        c = to_cnum(big(2)^70 * g + 0.5d)
        @test iszero(c - c)
        @test isequal(c + c, 2 * c)
        @test isequal(to_cnum(2) * c, c * to_cnum(2))
    end

    @testset "real symbolic quotients keep exact constants" begin
        @variables x::Real y::Real
        constants(e) =
            SymbolicUtils.iscall(e) ? mapreduce(constants, vcat, SymbolicUtils.arguments(e); init = Any[]) :
            SymbolicUtils.isconst(e) ? Any[SymbolicUtils.unwrap_const(e)] : Any[]
        for q in (1 / (1 + x), (y - 3x) / (1 + x^2), (1 // 4) * x / (1 + y), (y - 3x / y^4) / (1 + x^4 / y^8))
            for c in (convert(CNum, q), SecondQuantizedAlgebra.simplify(convert(CNum, q)))
                expr = c.tail.expr
                @test SymbolicUtils.symtype(expr) <: Real
                @test all(v -> v isa Union{Integer, Rational}, constants(expr))
                @test isequal(imag(to_num(c)), 0)
            end
        end
    end

    @testset "fractional powers preserve phase restrictions" begin
        @variables θ::Real
        phase = SymbolicUtils.Term{SymbolicUtils.SymReal}(expim, [Symbolics.unwrap(θ)]; type = Number)
        unit_root = Num(SymbolicUtils.Term{SymbolicUtils.SymReal}(sqrt, [phase]; type = Number))
        @test_throws ArgumentError unit_root * a
        scaled_root = Num(SymbolicUtils.Term{SymbolicUtils.SymReal}(sqrt, [2phase]; type = Number))
        retained = scaled_root * a
        @test isequal(get_variables(retained), [Symbolics.unwrap(θ)])
        @test occursin("sqrt", string(get_prefactor(retained)))
    end

    @testset "numeric radicals multiply exactly" begin
        inverse_root_six = 1 / sqrt(Num(6))
        weighted = inverse_root_six * a
        @test isequal(stored_coefficient(commutator(weighted, weighted')), to_cnum(1 // 6))

        @test isequal(to_cnum(sqrt(Num(2))) * to_cnum(sqrt(Num(3))), to_cnum(sqrt(Num(6))))
        @test hash(to_cnum(sqrt(Num(2))) * to_cnum(sqrt(Num(3)))) == hash(to_cnum(sqrt(Num(6))))
        @test isequal(to_cnum(sqrt(Num(8)))^2, to_cnum(8))
        @test isequal(to_cnum(1 / sqrt(Num(2))), to_cnum(sqrt(Num(2))) / 2)
        @test isequal(to_cnum(cbrt(Num(2)))^3, to_cnum(2))

        @variables g
        raw = to_cnum(sin(g + 1))
        half = to_cnum(1 / sqrt(Num(2)))
        @test isequal((half * raw) * half, to_cnum(1 // 2) * raw)

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

    @testset "radicals with prime factors beyond the trial bound stay canonical" begin
        @variables g k
        P = big(2)^127 - 1             # prime wider than Int
        q = 4294967311                 # smallest prime above 2^32
        q2 = 4294967357
        c = cos(g + k)
        # Issue #291: each spelling of the same value cancels.
        e = sqrt(Num(P)) * c * a
        @test iszero(e * e - (P * c^2) * a * a)
        f = sqrt(Num(q)) * c * a
        @test isequal(f * f, (q * c^2) * a * a)
        both = sqrt(Num(P)) * sqrt(Num(q)) * c * a
        @test iszero(both * both - (P * q * c^2) * a * a)

        root(n) = to_cnum(sqrt(Num(n)))
        # Composite radicands split by Pollard rho, below and above typemax(Int).
        # The first small semiprime makes the batched gcd overshoot to n and replay; the second
        # needs a second rho constant.
        for (n, parts) in (
                (big(65537) * 65551, (65537, 65551)), (big(65563) * 66413, (65563, 66413)),
                (big(q) * 65537 * 65539, (q, 65537, 65539)), (big(q) * q2, (q, q2)), (P * q, (P, q)),
            )
            @test isequal(root(n), prod(root, parts))
            @test hash(root(n)) == hash(prod(root, parts))
        end
        # Lowering keeps primes above the trial bound as separate roots, so it never builds a
        # radicand that only rho could split again.
        two_large = root(3040000009) * root(3040000039)
        @test isequal(to_cnum(real(to_num(two_large))), two_large)
        @test isequal(to_num(root(2) * root(3)), to_num(root(6)))
        # Perfect powers of a large prime fold into the scalar, in the big tier when needed.
        @test isequal(root(big(q)^3), q * root(q))
        @test isequal(root(P^3 * 2), P * root(P) * root(2))
        @test isequal(root(P)^2, to_cnum(P))
        @test isequal(to_cnum(sqrt(Num(1 // P))) * root(P), to_cnum(1))
        @test isequal(to_cnum(cbrt(Num(P)))^3, to_cnum(P))
        @test isequal(to_cnum(sqrt(Num(-P))), im * root(P))

        # Lowering keeps one root per wide prime and rebuilds the same coefficient.
        wide = sqrt(Num(P)) * sqrt(Num(q)) * g * a
        @test sort(split(string(real(get_prefactor(wide))), '*')) == sort(["g", "sqrt($P)", "sqrt($q)"])
        @test isequal(get_prefactor(wide) * a, wide)
    end

    @testset "cross-tier stress: prime-radical sums forced through the big tier" begin
        primes6 = [2, 3, 5, 7, 11, 13]
        s = sum(to_cnum(sqrt(Num(p))) for p in primes6)
        s2 = s^2
        big_scale = to_cnum(2)^70
        scaled = s2 * big_scale

        @test scaled.tail.terms isa Vector{Monomial{BigExactComplex}}

        built = sum(
            to_cnum(2)^70 * (to_cnum(sqrt(Num(p))) * to_cnum(sqrt(Num(q))))
                for p in primes6, q in primes6
        )
        @test isequal(scaled, built)
        @test hash(scaled) == hash(built)

        @test isequal(scaled * inv(big_scale), s2)
    end

    @testset "exact radicals and tiers survive lowering and public arithmetic" begin
        @variables g k
        roundtrip(e) = isequal(get_prefactor(e) * a, e)

        # Radicals lower to one grouped root per exponent and rebuild the same coefficient.
        six = sqrt(Num(2)) * sqrt(Num(3)) * a
        @test string(real(get_prefactor(six))) == "sqrt(6)"
        @test roundtrip(six)
        mixed = sqrt(Num(2)) * cbrt(Num(9)) * a
        @test string(real(get_prefactor(mixed))) == "(3^(2//3))*sqrt(2)"
        @test roundtrip(mixed)

        # Two radical primes whose product exceeds Int stay separate roots.
        p, q = 4294967291, 4294967279
        wide = sqrt(Num(p)) * sqrt(Num(q)) * a
        @test string(real(get_prefactor(wide))) == "sqrt($p)*sqrt($q)"
        @test roundtrip(wide)
        @test get_prefactor(wide * wide) == big(p) * q

        # Integer parts folded out of a radical overflow into the big tier.
        root = sqrt(Num(big(2)^127)) * a
        @test isequal(root, big(2)^63 * sqrt(Num(2)) * a)
        @test get_prefactor(root * root) == big(2)^127
        @test roundtrip(root)
        product = sqrt(Num(2^61)) * sqrt(Num(big(3)^61)) * g * a
        @test isequal(product, big(6)^30 * sqrt(Num(6)) * g * a)
        @test roundtrip(product)

        # A big Gaussian constant times a non-polynomial coefficient.
        huge = big(3)^50
        left = (huge * im * a) * (cos(g + k) * a)
        @test isequal(left, (cos(g + k) * a) * (huge * im * a))
        @test isequal(get_prefactor(left), Complex(Num(0), huge * cos(g + k)))
        @test isequal(get_prefactor(left) * a * a, left)
        both = ((huge + 2im) * a) * (cos(g + k) * a)
        @test isequal(get_prefactor(both) * a * a, both)

        # Roots of a negative radicand that cannot be factored: a product of two primes near
        # 2^64 is beyond the Pollard-rho budget.
        P = big(18446744073709551557) * big(18446744073709551533)
        real_root = cbrt(Num(-P)) * a
        @test isequal(get_prefactor(real_root), Complex(cbrt(Num(-P)), Num(0)))
        @test roundtrip(real_root)
        @test get_prefactor(sqrt(Num(-P)) * a) ≈ im * sqrt(Float64(P))

        # An unfactorable radical raised to a power stays a symbolic factor of the product.
        leaf = sqrt(Num(P)) * cos(g + k) * a
        partner = sin(g + k) * a
        @test isequal(leaf^3 * partner, partner * leaf^3)
        @test isequal(get_prefactor(leaf^3 * partner), get_prefactor(leaf^3) * sin(g + k))

        # Cube roots inside non-polynomial coefficients still combine exactly.
        root2 = cbrt(Num(2)) * cos(g + k) * a
        root4 = cbrt(Num(4)) * cos(g + k) * a
        @test isequal(root2 * root4, 2 * cos(g + k)^2 * a * a)
        @test isequal(root4 * root4, 2 * cbrt(Num(2)) * cos(g + k)^2 * a * a)

        # Only a unit-modulus phase has a magnitude, in either tier.
        @test abs(stored_coefficient((0.6 + 0.8im) * expim(g) * a)) == 1
        @test abs2(stored_coefficient((0.6 + 0.8im) * expim(g) * a)) == 1
        @test_throws MethodError abs(stored_coefficient(huge * expim(g) * a))

        # Floats combine with exact values in either order.
        half = stored_coefficient(0.5 * a)
        @test half / 2.0 == 0.25
        @test inv(half) == 2
        @test isequal(g * a / 3 + 0.5 * g * a, 0.5 * g * a + g * a / 3)
        mixed_sum = get_prefactor(substitute(g * a / 3 + 0.5 * g * a, Dict(g => 1)))
        @test Symbolics.value(real(mixed_sum)) ≈ 5 / 6
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

    @testset "quotients with different denominators cancel without simplify" begin
        @variables x::Real y::Real g w::Complex
        c = to_cnum

        @test iszero(c(g / (x + im * y)) * c(g / (x - im * y)) - c(g^2 / (x^2 + y^2)))
        @test iszero(c(1 / (x + 1) + 1 / (x + 2)) - c((2x + 3) / ((x + 1) * (x + 2))))
        @test iszero(c(sqrt(x) / (x + 1)) - c(x / (sqrt(x) * (x + 1))))
        @test iszero(c(1 / (w + 1)) * c(1 / (conj(w) + 1)) - c(1 / (w * conj(w) + w + conj(w) + 1)))
        @test iszero((g / (x + im * y)) * a * (g / (x - im * y)) - (g^2 / (x^2 + y^2)) * a)

        # A nonzero difference is never reported as zero.
        @test !iszero(c(1 / (x + 1)) - c(1 / (x + 2)))
        @test !iszero(c(1 / (x + 1)) - c(1 / (x + 1)) * (1 + 1 // 10^12))
        @test !iszero(c(g / (w + 1)) - c(g / (conj(w) + 1)))
        @test !iszero((g / (x + 1)) * a - (g / (x + 2)) * a)
    end

    @testset "steady-state amplitude of a damped oscillator cancels exactly" begin
        @variables q0::Real qd::Real γ::Real F::Real n::Real
        ω0 = 2q0^2
        Δa = 2qd^2 - ω0
        L = hamiltonian_action(-Δa * a' * a - F / (4q0) * (a + a')) +
            γ * (n + 1) * dissipator(a) + γ * n * dissipator(a')
        α = conjugate(a, DisplacementFrame(a, L)) - a
        d = α * α' / q0^2 - F^2 / (4ω0^2 * (Δa^2 + γ^2 / 4))
        @test iszero(d)
        @test !iszero(d + F^2 / (4ω0^2 * (Δa^2 + γ^2)))
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

    @testset "lowering and identity checks infer concrete types" begin
        @variables g::Real z
        for c in (to_cnum(2 // 3), to_cnum(g), to_cnum(z), to_cnum(Complex(g, g)), to_cnum(sqrt(g) * z))
            @test (@inferred to_num(c)) isa Complex{Num}
        end
        one_term = one(a)
        @test @inferred isone(one_term)
        @test isone(1.0 * one_term)
        @test !(@inferred isone(2 * one_term))
        @test !isone(g * one_term)
        @test !isone(a + 0)
    end
end
