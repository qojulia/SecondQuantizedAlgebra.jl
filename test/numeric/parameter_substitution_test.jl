using SecondQuantizedAlgebra
import SecondQuantizedAlgebra: QSym, get_prefactor
using QuantumOpticsBase
using Symbolics: @variables, substitute
using Test

dat(x) = dense(x).data

@testset "Numeric parameters and substitutions" begin
    @testset "to_numeric: Dict substitution" begin
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)

        custom_op = 2 * destroy(b)
        d = Dict(a => custom_op)

        @test to_numeric(a, b, d) == custom_op
        @test to_numeric(a', b, d) == create(b)

        result_mul = to_numeric(a' * a, b, d)
        @test dat(result_mul) == dat(create(b) * custom_op)
    end

    @testset "Number-symtype coefficient round-trip" begin
        @variables g::Number
        b = FockBasis(3)
        a = Destroy(FockSpace(:c), :a)

        for v in (2 + 3im, 1 + 1im, 0.5 - 0.25im)
            op = substitute((g * a)' * (g * a), Dict(g => v))
            @test dat(to_numeric(op, b)) ≈ abs2(v) * dat(to_numeric(a' * a, b))
        end

        @test dat(to_numeric(substitute(g * a, Dict(g => 2 + 3im)), b)) ≈ (2 + 3im) * dat(to_numeric(a, b))
    end

    @testset "Public coefficient lowering" begin
        h = FockSpace(:c)
        @qnumbers a::Destroy(h)
        @variables x::Real

        lowered = get_prefactor(x * a)
        @test isequal(real(lowered), x)
        @test iszero(imag(lowered))
    end

    @testset "to_numeric keyword translation" begin
        h = FockSpace(:c)
        @qnumbers a::Destroy(h)
        b = FockBasis(4)
        A = destroy(b)
        Ad = create(b)
        ψ = coherentstate(b, 0.3 - 0.1im)
        @variables x::Real ϕ::Real E::Number

        @test dat(to_numeric(substitute(sqrt(x) * a, Dict(x => 2.0)), b)) ≈ sqrt(2.0) * dat(A)
        @test dat(to_numeric(substitute(exp(im * ϕ) * a, Dict(ϕ => 0.5)), b)) ≈ exp(0.5im) * dat(A)

        @test dat(to_numeric(sqrt(x) * a, b; parameter = Dict(x => 2.0))) ≈ sqrt(2.0) * dat(A)
        @test dat(to_numeric(exp(im * ϕ) * a, b; parameter = Dict(ϕ => 0.5))) ≈ exp(0.5im) * dat(A)

        custom = 2 * A
        @test dat(to_numeric(a', b; operators = Dict(a => custom))) == dat(custom')
        @test dat(to_numeric(a', b; operators = Dict(a => custom), adjoint_ops = false)) == dat(Ad)
        @test [dat(x) for x in to_numeric([a, a'], b; operators = Dict(a => custom))] == [dat(custom), dat(custom')]

        f = to_numeric(E * a + conj(E) * a', b; time_parameter = Dict(E => t -> 1 + im * t))
        @test f isa TimeDependentSum
        @test expect(f(0.5), ψ) ≈ expect((1 + 0.5im) * A + (1 - 0.5im) * Ad, ψ)
    end

    @testset "keyword to_numeric on composite basis" begin
        hf = SecondQuantizedAlgebra.var"⊗"(FockSpace(:c), FockSpace(:d))
        a = Destroy(hf, :a, 1)
        @variables Δ::Real
        b = QuantumOpticsBase.var"⊗"(FockBasis(2), FockBasis(3))
        Ia = identityoperator(b)
        an = to_numeric(a, b)

        H = to_numeric(Δ * a' * a + 5, b; parameter = Dict(Δ => 2.0))
        @test size(dense(H).data) == (length(b), length(b))
        @test dense(H).data ≈ dense(2.0 * an' * an + 5 * Ia).data

        custom = 2 * an
        H2 = to_numeric(a + 3, b; operators = Dict(a => custom))
        @test dense(H2).data ≈ dense(custom + 3 * Ia).data

        z = to_numeric(a - a, b)
        @test size(dense(z).data) == (length(b), length(b))
        @test dense(z).data ≈ zeros(length(b), length(b))

        r = to_numeric(a', b; operators = Dict(a => custom))
        @test dense(r).data ≈ dense(custom').data
        r2 = to_numeric(a', b; operators = Dict(a => custom), adjoint_ops = false)
        @test dense(r2).data ≈ dense(to_numeric(a', b)).data
    end

    @testset "time_parameter variants" begin
        h = FockSpace(:c)
        @qnumbers a::Destroy(h)
        b = FockBasis(4)
        A = destroy(b)
        ψ = coherentstate(b, 0.3 - 0.1im)
        @variables E::Number

        f0 = to_numeric(E * a, b; time_parameter = Dict(E => 3.0))
        @test expect(f0(0.0), ψ) ≈ expect(3.0 * A, ψ)
        @test expect(f0(10.0), ψ) ≈ expect(3.0 * A, ψ)

        f1 = to_numeric(2.0 * a, b; time_parameter = Dict(E => t -> 1.0 + 0im))
        @test expect(f1(7.0), ψ) ≈ expect(2.0 * A, ψ)

        f2 = to_numeric(conj(E) * a, b; time_parameter = Dict(conj(E) => t -> 2 + im * t))
        @test expect(f2(1.0), ψ) ≈ expect((2 + 1im) * A, ψ)

        @test to_numeric(E * a, b; time_parameter = Dict(E => t -> 1), op_type = identity) isa TimeDependentSum
        @test_throws ArgumentError to_numeric(E * a, b; time_parameter = Dict(E => t -> 1), op_type = dense)
        @test_throws ArgumentError to_numeric(E * a, b; time_parameter = Dict(E => (p, t) -> p[1]))
    end

    @testset "keyword to_numeric, scalar argument" begin
        b = FockBasis(4)
        Ib = one(b)
        ψ = coherentstate(b, 0.3 - 0.1im)
        @variables x::Real E::Number

        @test to_numeric(3, b; parameter = Dict(x => 1.0)) == 3 * Ib
        @test to_numeric(2.0 + 0im, b; parameter = Dict{Any, Any}()) == (2.0 + 0im) * Ib
        @test dense(to_numeric(2 * x, b; parameter = Dict(x => 1.5))) ≈ 3.0 * Ib
        @test_throws ArgumentError to_numeric(x, b)

        fconst = to_numeric(2.0, b; time_parameter = Dict(E => t -> 1.0 + 0im))
        @test expect(fconst(0.0), ψ) ≈ expect(2.0 * Ib, ψ)
        @test expect(fconst(9.0), ψ) ≈ expect(2.0 * Ib, ψ)

        ft = to_numeric(E, b; time_parameter = Dict(E => t -> 1 + im * t))
        @test expect(ft(0.0), ψ) ≈ expect((1 + 0im) * Ib, ψ)
        @test expect(ft(2.0), ψ) ≈ expect((1 + 2im) * Ib, ψ)
    end

    @testset "keyword to_numeric, state argument" begin
        h = FockSpace(:c)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)
        α = 0.1 + 0.2im
        ψ = coherentstate(b, α)
        @variables Δ::Real

        op_state = to_numeric(Δ * a, ψ; parameter = Dict(Δ => 2.0))
        @test dat(op_state) == dat(2.0 * destroy(b))
    end

    @testset "constant symbolic coefficient reduction" begin
        h = FockSpace(:c)
        @qnumbers a::Destroy(h)
        b = FockBasis(3)
        A = destroy(b)
        @variables z::Number

        c0 = 1.0 + 2.0im
        cases = (
            ("real", real(conj(z)), real(conj(c0))),
            ("imag", imag(conj(z)), imag(conj(c0))),
            ("conj", conj(real(z) + im * imag(z)), conj(c0)),
            ("plus", real(z) + real(conj(z)), real(c0) + real(conj(c0))),
            ("times", real(z) * real(conj(z)), real(c0) * real(conj(c0))),
            ("div", real(z) / real(conj(z)), real(c0) / real(conj(c0))),
            ("pow", real(z)^3, real(c0)^3),
        )
        for (label, coeff, expected) in cases
            op = substitute(coeff * a, Dict(z => c0))
            @test dat(to_numeric(op, b)) ≈ expected * dat(A)
        end

        op_sin = substitute(sin(z) * a, Dict(z => 0.5))
        @test dat(to_numeric(op_sin, b)) ≈ sin(0.5) * dat(A)
    end

    @testset "keyword to_numeric: vector, complex params, errors" begin
        h = FockSpace(:c)
        @qnumbers a::Destroy(h)
        b = FockBasis(4)
        A = destroy(b)
        Ad = create(b)
        @variables x::Real z::Number E::Number

        Hd = to_numeric(2.0 * a' * a, b; op_type = dense)
        @test Hd isa QuantumOpticsBase.Operator
        @test Hd ≈ dense(2.0 * Ad * A)
        @test dense(to_numeric(2.0 * a' * a, b)) ≈ Hd

        @test [dat(x) for x in to_numeric([a, a'], b; parameter = Dict(x => 1.0))] == [dat(A), dat(Ad)]
        @test dat(to_numeric(z * a, b; parameter = Dict(z => 2 + 3im))) ≈ (2 + 3im) * dat(A)

        @test_throws ArgumentError to_numeric(x * a, b)
        @test_throws ArgumentError to_numeric(a, b; operators = Dict(x => A))
        @test_throws ArgumentError to_numeric(
            x * E * a, b; time_parameter = Dict(E => t -> 1.0 + 0im),
        )
    end

    @testset "complex parameter key" begin
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(6)
        @variables κ::Complex
        H = κ * a + conj(κ) * a'
        M = to_numeric(H, b; parameter = Dict(κ => 1.0 + 2.0im))
        Mref = to_numeric((1.0 + 2.0im) * a + (1.0 - 2.0im) * a', b)
        @test dat(M) ≈ dat(Mref)
    end

    @testset "time_parameter — unsupported key" begin
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(6)
        @variables E::Number
        @test_throws ArgumentError to_numeric(E * a, b; time_parameter = Dict(2 * E => t -> 1.0 + 0im))
    end
end
