using Test
using SecondQuantizedAlgebra
using Symbolics: Symbolics, @variables, Num
using LinearAlgebra: I, diag, diagm, norm, Diagonal, Hermitian, eigvals
using QuantumOpticsBase: FockBasis, NLevelBasis, SpinBasis
using Random: MersenneTwister
import SymbolicUtils
import SecondQuantizedAlgebra: UnitaryTransform, canonicality_residuals, ParamRelation,
    Op, QAdd, _to_cnum, _to_complex, _substitute_cnum, _zero_qadd,
    _single_qadd, _CNUM_ONE, _site_generators, _site_key, simplify, to_num, IndexedOperator,
    expim

# Unitary transforms are stored as their action on the generators of a site. These tests pin
# the closed forms, the `is_canonical` self-test that certifies them, the composition and
# inversion laws, the coverage checks that stop a half-transformed answer from escaping, and
# the gauge terms of the time-dependent constructors (the displacement one against a matrix
# oracle, since its disputed piece is a c-number that observable evolution cannot see).

@testset "Unitary transforms" begin
    h = FockSpace(:f)
    a = Destroy(h, :a)
    hab = FockSpace(:a) ⊗ FockSpace(:b)
    ma = Destroy(hab, :a, 1)
    mb = Destroy(hab, :b, 2)
    hph = PhaseSpace(:q)
    x = Position(hph, :x)
    p = Momentum(hph, :p)
    hsp = SpinSpace(:S)
    Sx, Sy, Sz = Spin(hsp, :S, 1), Spin(hsp, :S, 2), Spin(hsp, :S, 3)
    hpa = PauliSpace(:qb)
    σx, σy, σz = Pauli(hpa, :σ, 1), Pauli(hpa, :σ, 2), Pauli(hpa, :σ, 3)
    hnl = NLevelSpace(:atom, 2)
    σ12 = Transition(hnl, :σ, 1, 2)
    σ11 = Transition(hnl, :σ, 1, 1)
    σ22 = Transition(hnl, :σ, 2, 2)

    @variables θ φ r ϕ ω Ω Δ η g t u v dx dp
    @variables αc::Number

    @testset "closed forms" begin
        @test isequal(conjugate(a, Displace(a, αc)), a + αc)
        @test isequal(conjugate(a', Displace(a, αc)), a' + conj(αc))
        # A unit-modulus factor is one phase atom, not an Euler pair, so conjugating it is
        # exact and `is_canonical` needs no `atol` and no declared relation.
        @test isequal(conjugate(a, Rotation(a, θ)), expim(-θ) * a)
        @test isequal(conjugate(a', Rotation(a, θ)), expim(θ) * a')
        @test isempty(Rotation(a, θ).reductions)
        @test is_canonical(Rotation(a, θ); atol = 0)
        @test isequal(conjugate(a, Squeeze(a, r)), cosh(r) * a + sinh(r) * a')
        @test isequal(
            conjugate(a, Squeeze(a, r, ϕ)),
            cosh(r) * a + expim(ϕ) * sinh(r) * a',
        )
        @test isequal(conjugate(ma, Rotation(ma, mb, θ)), cos(θ) * ma + sin(θ) * mb)
        @test isequal(conjugate(mb, Rotation(ma, mb, θ)), cos(θ) * mb - sin(θ) * ma)
        @test isequal(conjugate(ma, Squeeze(ma, mb, r)), cosh(r) * ma + sinh(r) * mb')
        @test isequal(conjugate(ma, Bogoliubov(ma, mb, u, v)), u * ma + v * mb')
        @test isequal(conjugate(x, Displace(x, p, dx, dp)), x + dx)
        @test isequal(conjugate(p, Displace(x, p, dx, dp)), p + dp)
        @test isequal(conjugate(x, Rotation(x, p, θ)), cos(θ) * x + sin(θ) * p)
        @test isequal(conjugate(p, Rotation(x, p, θ)), cos(θ) * p - sin(θ) * x)
        @test isequal(conjugate(x, Squeeze(x, p, r)), exp(r) * x)
        @test isequal(conjugate(p, Squeeze(x, p, r)), (1 / exp(r)) * p)
        @test isequal(conjugate(Sx, Rotation(Sx, 3, θ)), cos(θ) * Sx - sin(θ) * Sy)
        @test isequal(conjugate(Sy, Rotation(Sx, 3, θ)), cos(θ) * Sy + sin(θ) * Sx)
        @test isequal(conjugate(Sz, Rotation(Sx, 3, θ)), 1 * Sz)
        @test isequal(conjugate(σx, Rotation(σx, 3, θ)), cos(θ) * σx - sin(θ) * σy)
        @test isequal(conjugate(σz, Rotation(σx, 3, θ)), 1 * σz)
        # Every rotation axis rotates the *other* two and fixes its own.
        for ax in 1:3
            @test isequal(conjugate(Spin(hsp, :S, ax), Rotation(Sx, ax, θ)), 1 * Spin(hsp, :S, ax))
        end
    end

    @testset "phase-space and Fock conventions agree" begin
        # `a = (x + im*p)/sqrt(2)`, so `x -> exp(r)*x` is the `ϕ = 0` single-mode squeeze and
        # `x -> cos(θ)*x + sin(θ)*p` is the Fock phase rotation.
        cr = conjugate(a, Squeeze(a, r))
        @test isequal(real(cr[Op[a]]), cosh(r))
        @test isequal(real(cr[Op[a']]), sinh(r))
        rot = conjugate(a, Rotation(a, θ))
        @test isequal(real(rot[Op[a]]), expim(-θ))
    end

    @testset "adjoint rules follow the forward ones" begin
        for U in (
                Displace(a, αc), Rotation(a, θ), Squeeze(a, r, ϕ),
                Rotation(ma, mb, θ), Squeeze(ma, mb, r), Bogoliubov(ma, mb, u, v),
            )
            for gen in generators(U)
                @test isequal(conjugate(adjoint(gen), U), adjoint(conjugate(gen, U)))
            end
        end
    end

    @testset "is_canonical" begin
        for (name, U) in (
                "Displace" => Displace(a, αc),
                "Rotation" => Rotation(a, θ),
                "Squeeze" => Squeeze(a, r),
                "Squeeze(r,ϕ)" => Squeeze(a, r, ϕ),
                "Rotation composite" => Rotation(a, ω * t),
                "beamsplitter" => Rotation(ma, mb, θ),
                "two-mode squeeze" => Squeeze(ma, mb, r),
                "Bogoliubov" => Bogoliubov(ma, mb, u, v),
                "phase Displace" => Displace(x, p, dx, dp),
                "phase Rotation" => Rotation(x, p, θ),
                "phase Squeeze" => Squeeze(x, p, r),
                "spin Rotation" => Rotation(Sx, 3, θ),
                "spin Rotation y" => Rotation(Sx, 2, θ),
                "Pauli Rotation" => Rotation(σx, 1, θ),
                "NLevel Rotation" => Rotation(σ12, [cos(θ) -sin(θ); sin(θ) cos(θ)]),
                "NLevel swap" => Rotation(σ12, [0 1; 1 0]),
                "composed" => Rotation(a, θ) * Squeeze(a, r, ϕ),
            )
            @testset "$name" begin
                @test is_canonical(U)
            end
        end
    end

    @testset "is_canonical rejects a non-unitary map" begin
        # `a -> 2a` preserves linearity but not the commutator.
        rules = Dict{Op, QAdd}(a => 2 * a, a' => 2 * a')
        U = UnitaryTransform(rules, rules, _zero_qadd(), Num(0), ParamRelation[], Num[])
        @test !is_canonical(U)
        @test any(!iszero, canonicality_residuals(U))
        # A rotation that forgets to conjugate the phase breaks Hermiticity.
        bad = Dict{Op, QAdd}(a => cos(θ) * a, a' => sin(θ) * a')
        Ub = UnitaryTransform(bad, bad, _zero_qadd(), Num(0), ParamRelation[], Num[])
        @test !is_canonical(Ub)
        # A correct forward map with a wrong stored inverse: preserving the algebra says
        # nothing about `inv_rules`, so the round-trip residuals have to catch it.
        fwd = Dict{Op, QAdd}(
            a => (cos(θ) - im * sin(θ)) * a, a' => (cos(θ) + im * sin(θ)) * a',
        )
        Uf = UnitaryTransform(fwd, fwd, _zero_qadd(), Num(0), ParamRelation[], Num[])
        @test !is_canonical(Uf)
        @test !isequal(conjugate(conjugate(a, Uf), inv(Uf)), 1 * a)
        # A complex 3-level matrix that is not unitary fails as well.
        h3 = NLevelSpace(:three, 3)
        s3 = Transition(h3, :τ, 1, 2)
        @test !is_canonical(Rotation(s3, [1 1 0; 0 1 0; 0 0 1]))
        # `x -> x + im` preserves the CCR but not Hermiticity, so the phase-space residuals
        # must carry the Hermiticity test and not the commutator alone.
        @test !is_canonical(Displace(x, p, im, 0))
        @test any(!iszero, canonicality_residuals(Displace(x, p, im, 0)))
        @test is_canonical(Displace(x, p, dx, dp))
        # The N-level products are tested against the first row and column only, which is
        # equivalent to the full `n^4` table. Every one of these is rejected either way.
        for W in (
                [1 1 0; 0 1 0; 0 0 1], 2.0 * Matrix(I, 3, 3), Float64[1 0 0; 0 1 0; 0 0 2],
                Float64[1 0 0; 1 1 0; 0 0 1], Float64[1 0 0; 0 0 1; 0 1.0001 0],
            )
            @test !is_canonical(Rotation(s3, W))
        end
        for W in (
                Float64[0.6 0.8 0; 0.8 -0.6 0; 0 0 1], Float64[0 1 0; 0 0 1; 1 0 0],
            )
            @test is_canonical(Rotation(s3, W))
        end
        # `atol` waves through a residual it can *bound* below the tolerance. A `cos`/`sin`
        # factor is bounded by one, so a tiny scalar in front of it really is rounding; an
        # unbounded factor is not admitted whatever its scalar.
        @test SecondQuantizedAlgebra._bounded_by(_to_cnum(1.0e-17 * cos(ω * t)), 1.0e-12)
        @test !SecondQuantizedAlgebra._bounded_by(_to_cnum(1.0e-17 * sinh(r)), 1.0e-12)
        @test !SecondQuantizedAlgebra._bounded_by(_to_cnum(1.0e-17 * ω), 1.0e-12)
        @test !SecondQuantizedAlgebra._bounded_by(_to_cnum(0.5 * cos(ω * t)), 1.0e-12)
    end

    @testset "finite-dimensional oracle, complex unitary" begin
        # `U = Σ W[k,l]*σ^{kl}` built as a `QAdd` must reproduce `conjugate(A, U)` for every
        # matrix unit. The contraction runs over W's row index, and the transpose does not
        # agree, so the direction is pinned here rather than assumed.
        h3 = NLevelSpace(:three, 3)
        s3 = Transition(h3, :τ, 1, 2)
        c, s = 0.6, 0.8
        W = ComplexF64[
            c (s * im) 0
            (s * im) c 0
            0 0 (0.6 + 0.8im)
        ]
        @test norm(W' * W - I) < 1.0e-14
        Uq = _zero_qadd()
        for k in 1:3, l in 1:3
            iszero(W[k, l]) && continue
            Uq = Uq + W[k, l] * Transition(h3, :τ, k, l)
        end
        Ut = Rotation(s3, W)
        for i in 1:3, j in 1:3
            op = Transition(h3, :τ, i, j)
            @test isequal(simplify(conjugate(op, Ut)), simplify(Uq' * op * Uq))
        end
        @test is_canonical(Ut)
        # ... but only to rounding: the exact test rejects a Float64 matrix.
        @test !is_canonical(Ut; atol = 0)
        @test all(q -> all(p -> abs(p.second.z) < 1.0e-14, q.arguments), canonicality_residuals(Ut))
    end

    @testset "round trip" begin
        for U in (
                Displace(a, αc), Rotation(a, θ), Squeeze(a, r, ϕ),
                Rotation(x, p, θ), Squeeze(x, p, r), Rotation(Sx, 3, θ),
                Rotation(σ12, [cos(θ) -sin(θ); sin(θ) cos(θ)]),
                Squeeze(ma, mb, r), Bogoliubov(ma, mb, u, v),
            )
            for gen in generators(U)
                @test isequal(simplify(conjugate(conjugate(gen, U), inv(U))), 1 * gen)
            end
        end
        # Non-commuting factors: `inv_rules` composes in the reverse order, which a single-`U`
        # round trip cannot catch.
        Uc = Rotation(a, θ) * Squeeze(a, r, ϕ)
        @test isequal(simplify(conjugate(conjugate(a, Uc), inv(Uc))), 1 * a)
        @test isequal(simplify(conjugate(conjugate(a' * a, Uc), inv(Uc))), a' * a)
        @test isequal(inv(inv(Uc)).rules, Uc.rules)
        @test isequal(adjoint(Uc).rules, inv(Uc).rules)
    end

    @testset "composition" begin
        U1 = Rotation(a, θ)
        U2 = Squeeze(a, r, ϕ)
        H = a' * a + a + a'
        @test isequal(simplify(transform(H, U1 * U2)), simplify(transform(transform(H, U1), U2)))
        # ... and with a time-dependent factor, where the gauge has to compose too.
        Ut = Rotation(a, ω * t, t)
        Us = Squeeze(a, r, ϕ)
        @test isequal(
            simplify(transform(H, Ut * Us)), simplify(transform(transform(H, Ut), Us))
        )
        # A static factor adopts the other's time variable; genuinely different ones throw.
        @test isequal((Ut * Us).time, t)
        @test isequal((Us * Ut).time, t)
        @test_throws ArgumentError Rotation(a, ω * t, t) * Rotation(a, ω * φ, φ)
        # Two sites compose into one transform covering both.
        Uab = Rotation(ma, θ) * Rotation(mb, φ)
        @test length(generators(Uab)) == 4
        @test isequal(simplify(conjugate(ma' * ma * mb' * mb, Uab)), ma' * ma * mb' * mb)
    end

    @testset "coverage" begin
        # A second named Fock mode on the same space is a distinct site and needs no rule.
        b = Destroy(h, :b)
        @test isequal(conjugate(b' * b, Rotation(a, θ)), b' * b)
        # A partial rule set cannot even be constructed.
        part = Dict{Op, QAdd}(σ12 => 1 * σ12, σ12' => 1 * σ12')
        @test_throws ArgumentError UnitaryTransform(
            part, part, _zero_qadd(), Num(0), ParamRelation[], Num[]
        )
        # A phase-space rule set without the conjugate variable likewise.
        half = Dict{Op, QAdd}(x => 1 * x)
        @test_throws ArgumentError UnitaryTransform(
            half, half, _zero_qadd(), Num(0), ParamRelation[], Num[]
        )
        # A second operator family sharing a covered site's key throws instead of coming
        # back untransformed. A `Pauli` and a `Spin` with the same name on subspace 1 have
        # the same site key but generate different algebras.
        hpq = PauliSpace(:S)
        @test isequal(_site_key(Pauli(hpq, :S, 1)), _site_key(Sx))
        @test_throws ArgumentError conjugate(Pauli(hpq, :S, 1), Rotation(Sx, 3, θ))
        # A differently named `Position` on a covered site throws: `_site_compare` treats it
        # as the conjugate variable of `p`, so letting it through returns a half-transformed
        # expression with `p` displaced and it left alone.
        y = Position(hph, :y)
        @test isequal(_site_key(y), _site_key(x))
        @test_throws ArgumentError conjugate(y * p, Displace(x, p, dx, dp))
        # ... and says so, rather than asking for a rule that cannot exist.
        @test occursin(
            "one canonical pair", (
                try
                    conjugate(y * p, Displace(x, p, dx, dp))
                catch e
                    sprint(showerror, e)
                end
            )
        )
        # An indexed operator on a transformed subspace throws rather than passing through.
        i = Index(h, :i, 4, h)
        ai = IndexedOperator(a, i)
        @test_throws ArgumentError conjugate(ai' * ai, Rotation(a, θ))
        @test_throws ArgumentError conjugate(Σ(ai' * ai, i), Rotation(a, θ))
        # A transform keyed on a *free* index covers the whole family instead, so a
        # differently indexed member is instantiated rather than refused.
        j = Index(h, :j, 4, h)
        aj = IndexedOperator(a, j)
        Ui = Rotation(ai, θ)
        @test isequal(conjugate(ai, Ui), expim(-θ) * ai)
        @test isequal(conjugate(aj' * aj, Ui), aj' * aj)
        # The unindexed operator is still a different operator, not the representative.
        @test occursin(
            "carries no index", (
                try
                    conjugate(a, Ui)
                catch e
                    sprint(showerror, e)
                end
            )
        )
        # A per-slot index is a resolved site, so it keeps the exact semantics.
        U3 = Rotation(IndexedOperator(a, i(3)), θ)
        @test isequal(conjugate(IndexedOperator(a, i(3)), U3), expim(-θ) * IndexedOperator(a, i(3)))
        @test_throws ArgumentError conjugate(IndexedOperator(a, i(4)), U3)
        # An indexed family of an *unrelated* mode is not on a covered site, so it passes.
        bi = IndexedOperator(b, i)
        @test isequal(conjugate(bi' * bi, Rotation(a, θ)), bi' * bi)
        # Eager canonicalization turns `σ12*σ21` into the leaf `σ11`, so the N-level rule set
        # has to carry the diagonal: it does, and the result is not the input.
        Un = Rotation(σ12, [cos(θ) -sin(θ); sin(θ) cos(θ)])
        @test !isequal(conjugate(σ12 * σ12', Un), 1 * σ11)
        @test isequal(
            conjugate(σ11, Un),
            cos(θ)^2 * σ11 - cos(θ) * sin(θ) * σ12 - cos(θ) * sin(θ) * σ12' + sin(θ)^2 * σ22,
        )
    end

    @testset "accessors" begin
        U = Bogoliubov(ma, mb, u, v)
        @test isequal(constraints(U), Num[u^2 - v^2 - 1])
        @test isempty(constraints(Squeeze(ma, mb, r)))
        @test iszero(gauge_term(Rotation(a, θ)))
        @test generators(Rotation(a, θ)) == sort!(Op[a, a'])
        @test length(generators(Rotation(σ12, [0 1; 1 0]))) == 4
        @test occursin("UnitaryTransform", string(Rotation(a, θ)))
    end

    @testset "site generators" begin
        @test isequal(_site_generators(a), Op[a, a'])
        @test isequal(_site_generators(a'), Op[a, a'])
        @test isequal(_site_generators(Sx), Op[Sx, Sy, Sz])
        @test isequal(_site_generators(σy), Op[σx, σy, σz])
        @test length(_site_generators(σ12)) == 4
        @test isequal(_site_generators(σ12), Op[σ11, σ12, σ12', σ22])
        # Not recoverable for phase space: `x` and `p` carry different names, so the site
        # cannot describe itself from one member and reports an empty set. Those sites are
        # validated against the generators the caller actually supplied instead.
        @test isempty(_site_generators(x))
        @test isempty(_site_generators(p))
        # The name is part of the site key, so two Fock modes on one space differ.
        @test !isequal(_site_key(a), _site_key(Destroy(h, :b)))
    end

    @testset "constructor validation" begin
        @test_throws ArgumentError Displace(Sx, θ)
        @test_throws ArgumentError Rotation(Sx, 0, θ)
        @test_throws ArgumentError Rotation(Sx, 4, θ)
        # A Fock mode has no axis, so `(Op, scalar, scalar)` is always `(a, θ, t)` there.
        # This used to depend on whether the angle was spelled `2` or `2.0`.
        @test isequal(
            conjugate(a, Rotation(a, 2, t)), conjugate(a, Rotation(a, 2.0, t))
        )
        @test_throws ArgumentError Rotation(σ12, [1 0 0; 0 1 0; 0 0 1])
        @test_throws ArgumentError Rotation(a, [0 1; 1 0])
        @test_throws ArgumentError Rotation(ma, ma, θ)        # needs two distinct modes
        @test_throws ArgumentError Displace(p, x, dx, dp)     # wrong pair order
        # Validation errors name the constructor the caller actually used.
        @test occursin(
            "`Bogoliubov`", (
                try
                    Bogoliubov(ma, ma, u, v)
                catch e
                    sprint(showerror, e)
                end
            )
        )
        @test occursin(
            "`Squeeze`", (
                try
                    Squeeze(ma, ma, r)
                catch e
                    sprint(showerror, e)
                end
            )
        )
        # A numeric parameter folds a relation's members to literals, which must never reach
        # the reduction: `Squeeze(a, 0.5)` has `cosh(r)`/`sinh(r)` as `1.13`/`0.52`.
        @test all(SecondQuantizedAlgebra._is_usable_rel, Squeeze(a, r).reductions)
        @test all(SecondQuantizedAlgebra._is_usable_rel, Squeeze(a, 0.5).reductions)
        @test is_canonical(Squeeze(a, r))
        @test isequal(conjugate(a, Squeeze(a, r)), cosh(r) * a + sinh(r) * a')
    end

    @testset "gauge terms" begin
        # Static transforms have no gauge, so `transform` and `conjugate` agree.
        for U in (Displace(a, αc), Rotation(a, θ), Squeeze(ma, mb, r))
            gen = first(generators(U))
            @test isequal(transform(gen, U), conjugate(gen, U))
        end
        # `-∂ₜθ*G` for the rotation family.
        @test isequal(gauge_term(Rotation(a, ω * t, t)), -ω * a' * a)
        @test isequal(gauge_term(Rotation(Sx, 3, ω * t, t)), -ω * Sz)
        @test isequal(gauge_term(Rotation(σx, 3, ω * t, t)), (-1 // 2) * ω * σz)
        @test isequal(gauge_term(Rotation(x, p, ω * t, t)), (-1 // 2) * ω * (x * x + p * p))
        @test isequal(gauge_term(Rotation(ma, mb, ω * t, t)), -im * ω * (ma' * mb - mb' * ma))
        @test isequal(gauge_term(Squeeze(ma, mb, ω * t, t)), -im * ω * (ma' * mb' - mb * ma))
        # A time-independent parameter leaves the gauge zero even on the `t` path.
        @test iszero(gauge_term(Rotation(a, θ, t)))
        # `inv` negates the gauge for the commuting family.
        @test isequal(gauge_term(inv(Rotation(a, ω * t, t))), ω * a' * a)
        # An N-level frame's gauge is `-Σ Eₖ σ^{kk}`.
        Un = RotatingFrame(Δ * σ22, t)
        @test isequal(gauge_term(Un), -Δ * σ22)
        # `Rotation(σ, W, t)` derives that same gauge from `im*Ẇ'W`, for any `W`.
        Wd = Complex{Num}[cos(Δ * t) - im * sin(Δ * t) 0; 0 1]
        @test isequal(gauge_term(Rotation(σ12, Wd, t)), -Δ * σ11)
        @test is_canonical(Rotation(σ12, Wd, t))
        # A fourth argument still overrides it.
        @test isequal(gauge_term(Rotation(σ12, Wd, t, Ω * σ22)), Ω * σ22)
        # A static `Squeeze` phase keeps the `-∂ₜr*G` form of the rotation family.
        Gsq = im * (expim(ϕ) * a' * a' - expim(-ϕ) * a * a) / 2
        @test iszero(simplify(gauge_term(Squeeze(a, ω * t, ϕ, t)) + ω * Gsq))
        @test iszero(gauge_term(Squeeze(a, r, ϕ, t)))
        # A moving one adds the `a'*a + 1/2` piece that `-∂ₜθ*G` cannot carry.
        want = ω * sinh(r)^2 * (a' * a + 1 // 2) +
            (ω / 2) * sinh(r) * cosh(r) * (
            expim(ω * t) * a' * a' + expim(-ω * t) * a * a
        )
        @test iszero(simplify(gauge_term(Squeeze(a, r, ω * t, t)) - want))
    end

    @testset "Displace gauge against a matrix oracle" begin
        # The disputed piece of this gauge is a c-number, which contributes only a global
        # phase, so no observable-evolution test can see it. Compare the closed form with a
        # central difference of `im*(∂ₜD†)D` in a truncated Fock basis instead.
        N = 60
        A = diagm(1 => sqrt.(1.0:(N - 1)))
        Ad = collect(A')
        Ione = Matrix{ComplexF64}(I, N, N)
        keep = 5:45
        τ0, dτ = 0.41, 1.0e-5

        # Parameterized on the identity so the deeper squeeze truncation reuses it.
        function materialize(q, τ, ops, Id = Ione)
            M = zeros(ComplexF64, size(Id)...)
            for (term, c) in q
                z = _to_complex(_substitute_cnum(c, Dict(t => τ)))
                blk = copy(Id)
                for o in term.ops
                    blk = blk * ops[o.kind]
                end
                M .+= z .* blk
            end
            return M
        end
        cnumber(q, τ) = sum(
            (
                isempty(term.ops) ? _to_complex(_substitute_cnum(c, Dict(t => τ))) :
                    0.0 + 0.0im for (term, c) in q
            );
            init = 0.0 + 0.0im,
        )
        oracle(f, τ) = 1im * ((f(τ + dτ)' - f(τ - dτ)') / (2dτ)) * f(τ)

        @testset "Fock" begin
            αnum(τ) = 0.35 * exp(0.7im * τ)
            D(τ) = exp(Matrix(αnum(τ) * Ad - conj(αnum(τ)) * A))
            gq = gauge_term(Displace(a, 0.35 * (cos(0.7 * t) + im * sin(0.7 * t)), t))
            G = materialize(gq, τ0, Dict(SecondQuantizedAlgebra.OP_DESTROY => A, SecondQuantizedAlgebra.OP_CREATE => Ad))
            ref = oracle(D, τ0)
            z = cnumber(gq, τ0)
            @test norm(G[keep, keep] - ref[keep, keep]) < 1.0e-7
            # The c-number is load-bearing: flipping or dropping it fails.
            @test norm((G - 2z * Ione)[keep, keep] - ref[keep, keep]) > 0.1
            @test norm((G - z * Ione)[keep, keep] - ref[keep, keep]) > 0.1
        end

        @testset "phase space" begin
            X = (A + Ad) / sqrt(2)
            P = (A - Ad) / (1im * sqrt(2))
            dxn(τ) = 0.4 * cos(1.3τ)
            dpn(τ) = 0.25 * sin(0.9τ)
            T(τ) = exp(Matrix(1im * (dpn(τ) * X - dxn(τ) * P)))
            gq = gauge_term(Displace(x, p, 0.4 * cos(1.3 * t), 0.25 * sin(0.9 * t), t))
            G = materialize(gq, τ0, Dict(SecondQuantizedAlgebra.OP_POSITION => X, SecondQuantizedAlgebra.OP_MOMENTUM => P))
            ref = oracle(T, τ0)
            z = cnumber(gq, τ0)
            @test norm(G[keep, keep] - ref[keep, keep]) < 1.0e-7
            @test norm((G - 2z * Ione)[keep, keep] - ref[keep, keep]) > 0.1
            @test norm((G - z * Ione)[keep, keep] - ref[keep, keep]) > 0.1
        end

        @testset "squeeze with a moving phase" begin
            # `[A, ∂ₜA]` is an operator here, not a c-number, so nothing about this gauge is
            # a global phase: the oracle pins the `a'*a` rate and the `1/2` offset alike.
            # A squeeze spreads far up the ladder, so this needs its own deeper truncation.
            Nsq = 160
            Asq = diagm(1 => sqrt.(1.0:(Nsq - 1)))
            Adsq = collect(Asq')
            Isq = Matrix{ComplexF64}(I, Nsq, Nsq)
            ksq = 1:40
            rn(τ) = 0.3 + 0.17τ
            qn(τ) = 0.9 - 0.23τ + 0.11τ^2
            zn(τ) = rn(τ) * cis(qn(τ))
            Sq(τ) = exp(Matrix((zn(τ) * Adsq^2 - conj(zn(τ)) * Asq^2) / 2))
            gq = gauge_term(
                Squeeze(a, 0.3 + 0.17 * t, 0.9 - 0.23 * t + 0.11 * t^2, t)
            )
            G = materialize(
                gq, τ0,
                Dict(
                    SecondQuantizedAlgebra.OP_DESTROY => Asq,
                    SecondQuantizedAlgebra.OP_CREATE => Adsq,
                ),
                Isq,
            )
            ref = oracle(Sq, τ0)
            z = cnumber(gq, τ0)
            @test abs(z) > 1.0e-3
            @test norm(G[ksq, ksq] - ref[ksq, ksq]) < 1.0e-5
            @test norm((G - z * Isq)[ksq, ksq] - ref[ksq, ksq]) > 1.0e-3
        end
    end

    @testset "indexed families" begin
        hf = FockSpace(:ff)
        af = Destroy(hf, :af)
        ii = Index(hf, :ii, 4, hf)
        jj = Index(hf, :jj, 4, hf)
        ai, aj = IndexedOperator(af, ii), IndexedOperator(af, jj)
        # A family transform is `⊗ᵢUᵢ`, so instantiating the rules at another index has to
        # agree with building the same constructor there directly.
        for (name, U, V, gens) in (
                (
                    "Displace", Displace(ai, IndexedVariable(:αi, ii)),
                    Displace(aj, IndexedVariable(:αi, jj)), (aj, aj'),
                ),
                ("Rotation", Rotation(ai, θ), Rotation(aj, θ), (aj, aj')),
                (
                    "Squeeze", Squeeze(ai, IndexedVariable(:ri, ii)),
                    Squeeze(aj, IndexedVariable(:ri, jj)), (aj, aj'),
                ),
                ("Squeeze(r, ϕ)", Squeeze(ai, r, ϕ), Squeeze(aj, r, ϕ), (aj, aj')),
            )
            @testset "$name" begin
                @test is_canonical(U)
                for gen in gens
                    @test isequal(conjugate(gen, U), conjugate(gen, V))
                end
            end
        end
        # The relations move with the rules. Without that a family squeeze declares
        # `cosh(r(i))^2 - sinh(r(i))^2 = 1` and then fails to close at every other index.
        Us = Squeeze(ai, IndexedVariable(:rs, ii))
        @test iszero(simplify(commutator(conjugate(aj, Us), conjugate(aj', Us)) - 1))
        # Instances live on disjoint sites, so they commute after the transform with no
        # residual needed for it. Two free indices are `Undetermined` until declared unequal.
        @test iszero(
            assume_distinct_index(
                commutator(conjugate(ai, Us), conjugate(aj', Us)), [(ii, jj)]
            )
        )
        # Inside a sum the bound index is just another instantiation target.
        @test isequal(conjugate(Σ(aj' * aj, jj), Rotation(ai, θ)), Σ(aj' * aj, jj))
        # A family rule may not leave its site: `a_i -> cosh*a_i + sinh*b'` would pass every
        # same-site residual while breaking `[ã_i, ã_j']`.
        hf2 = FockSpace(:fa) ⊗ FockSpace(:fb)
        m1, m2 = Destroy(hf2, :m1, 1), Destroy(hf2, :m2, 2)
        i2 = Index(hf2, :i2, 4, FockSpace(:fa))
        @test occursin(
            "another site", (
                try
                    Rotation(IndexedOperator(m1, i2), IndexedOperator(m2, i2), θ)
                catch e
                    sprint(showerror, e)
                end
            )
        )
        # A time-dependent family gauge is the sum over the family, not one site's worth.
        Ut = Rotation(ai, ω * t, t)
        @test isequal(gauge_term(Ut), Σ(-ω * ai' * ai, ii))
        @test iszero(transform(Σ(ω * ai' * ai, ii), Ut))
        # The frame of an indexed `H0` is a family with no extra machinery: `ad` already
        # closes on the family generators.
        Uf = RotatingFrame(Σ(ω * ai' * ai, ii), t)
        @test is_canonical(Uf; atol = 0)
        @test isequal(conjugate(aj, Uf), expim(-ω * t) * aj)
        @test iszero(transform(Σ(ω * ai' * ai, ii), Uf))
        # ...including with a per-site frequency.
        Uω = RotatingFrame(Σ(IndexedVariable(:ωi, ii) * ai' * ai, ii), t)
        @test is_canonical(Uω; atol = 0)
        @test isequal(conjugate(aj, Uω), expim(-IndexedVariable(:ωi, jj) * t) * aj)
        # A frame whose blocks mix two sites of one index cannot be a `⊗ᵢUᵢ`, and the family
        # check is what catches it before it can answer wrongly.
        @test_throws ArgumentError RotatingFrame(
            Σ(
                g * (
                    IndexedOperator(m1, i2)' * IndexedOperator(m2, i2) +
                        IndexedOperator(m2, i2)' * IndexedOperator(m1, i2)
                ), i2,
            ), t,
        )
    end

    @testset "collective N-level" begin
        hcol = CollectiveNLevelSpace(:atc, 3)
        S(i, j) = CollectiveTransition(hcol, :S, i, j)
        Symbolics.@variables E1 E2 E3
        # The `U'S^{ij}U = Σ conj(W[i,k])*W[j,l]*S^{kl}` contraction holds verbatim for a
        # collective operator: it is linear in the atom sum and never touches the matrix-unit
        # product law that collective transitions do not obey.
        Wswap = [0.0 1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 1.0]
        Uswap = Rotation(S(1, 2), Wswap)
        @test is_canonical(Uswap)
        @test isequal(conjugate(S(1, 2), Uswap), 1 * S(2, 1))
        @test isequal(conjugate(S(3, 3), Uswap), 1 * S(3, 3))
        # A symbolic rotation in the 1-2 plane, checked against the contraction by hand.
        Wrot = [cos(θ) -sin(θ) 0; sin(θ) cos(θ) 0; 0 0 1]
        Urot = Rotation(S(1, 2), Wrot)
        @test is_canonical(Urot)
        want = _zero_qadd()
        for k in 1:3, l in 1:3
            want = want + (Wrot[1, k] * Wrot[2, l]) * S(k, l)
        end
        @test isequal(simplify(conjugate(S(1, 2), Urot)), simplify(want))
        # A non-unitary `W` fails, and the outer automorphism `S^{ij} -> -S^{ji}` is the
        # reason the number-invariance residual exists: it satisfies the bracket law, is
        # Hermitian-consistent and is involutive, so every other residual passes it.
        @test !is_canonical(Rotation(S(1, 2), diagm([1.0, 2.0, 1.0])))
        contra = Dict{Op, QAdd}(S(i, j) => -1 * S(j, i) for i in 1:3, j in 1:3)
        Ucon = UnitaryTransform(
            contra, contra, _zero_qadd(), Num(0), ParamRelation[], Num[],
        )
        @test !is_canonical(Ucon)
        @test count(!iszero, canonicality_residuals(Ucon)) == 1
        # A partial rule set is refused: the level count is on the space, so the keys are the
        # only witness of what the transform claims to cover.
        partial = Dict{Op, QAdd}(S(1, 1) => 1 * S(1, 1), S(1, 2) => 1 * S(1, 2))
        @test occursin(
            "full square", (
                try
                    UnitaryTransform(
                        partial, partial, _zero_qadd(), Num(0), ParamRelation[], Num[],
                    )
                catch e
                    sprint(showerror, e)
                end
            )
        )
        # A level-diagonal collective frame carries one phase per level, so the commensurate
        # frequencies of three levels close by exponent arithmetic and it certifies exactly.
        Hd = E1 * S(1, 1) + E2 * S(2, 2) + E3 * S(3, 3)
        Ud = RotatingFrame(Hd, t)
        @test is_canonical(Ud; atol = 0)
        @test isequal(conjugate(S(1, 3), Ud), expim(E1 * t) * expim(-E3 * t) * S(1, 3))
        @test isequal(gauge_term(Ud), -Hd)
        @test iszero(transform(Hd, Ud))
        # Off the level diagonal there is no dressed factorization to fall back on.
        @test occursin(
            "Rotation(S, W)", (
                try
                    RotatingFrame(E1 * S(1, 1) + Ω * (S(1, 2) + S(2, 1)), t)
                catch e
                    sprint(showerror, e)
                end
            )
        )
        # Alongside another site the generic path takes over, and the levels it covers are
        # inferred from `H0` since the operator does not carry them.
        hmix = CollectiveNLevelSpace(:atm, 2) ⊗ FockSpace(:cavm)
        Sm(i, j) = CollectiveTransition(hmix, :Sm, i, j)
        am = Destroy(hmix, :am, 2)
        Hmix = E1 * Sm(2, 2) + ω * am' * am
        Umix = RotatingFrame(Hmix, t)
        @test is_canonical(Umix)
        @test iszero(transform(Hmix, Umix))
        @test isequal(conjugate(Sm(1, 2), Umix), expim(-E1 * t) * Sm(1, 2))
        # Tavis-Cummings leaves the linear span, and is named as such.
        @test_throws ArgumentError RotatingFrame(
            E1 * Sm(2, 2) + ω * am' * am + Ω * (am' * Sm(1, 2) + am * Sm(2, 1)), t
        )
        # A collective operator cannot be indexed, so no family of them exists either.
        @test_throws ArgumentError IndexedOperator(S(1, 2), Index(hcol, :ic, 3, hcol))
    end

    @testset "operator-valued displacement" begin
        hpol = FockSpace(:cav) ⊗ FockSpace(:mech)
        ac, bm = Destroy(hpol, :ac, 1), Destroy(hpol, :bm, 2)
        Symbolics.@variables ω₀ ωₘ gp
        A = (-gp / ωₘ) * (ac' * ac)
        U = Displace(bm, A)
        @test is_canonical(U)
        @test isequal(conjugate(bm, U), bm - (gp / ωₘ) * ac' * ac)
        # The polaron transform in full: the coupling is eliminated and a Kerr term appears.
        H = ω₀ * ac' * ac + ωₘ * bm' * bm + gp * ac' * ac * (bm + bm')
        @test isequal(
            conjugate(H, U),
            (ω₀ - gp^2 / ωₘ) * ac' * ac + ωₘ * bm' * bm -
                (gp^2 / ωₘ) * ac' * ac' * ac * ac,
        )
        # The amplitude's site is gated, not free: what commutes with it passes, and what
        # does not would have a displacement operator for an image, so it throws.
        @test isequal(conjugate(ac' * ac, U), ac' * ac)
        @test occursin(
            "displacement operator", (
                try
                    conjugate(ac, U)
                catch e
                    sprint(showerror, e)
                end
            )
        )
        @test isequal(conjugate(conjugate(bm, U), inv(U)), 1 * bm)
        # A spin-conditional displacement is the same class, and the case where passing the
        # amplitude's site through untouched is provably wrong: `U'σ¹²U = σ¹²*D(2λ)`.
        hcd = FockSpace(:fc) ⊗ NLevelSpace(:atc, 2)
        acd = Destroy(hcd, :acd, 1)
        τ11 = Transition(hcd, :τc, 1, 1, 2)
        τ22 = Transition(hcd, :τc, 2, 2, 2)
        τ12 = Transition(hcd, :τc, 1, 2, 2)
        Uz = Displace(acd, r * (τ11 - τ22))
        @test is_canonical(Uz)
        @test isequal(conjugate(acd, Uz), acd + r * τ11 - r * τ22)
        @test isequal(conjugate(τ11, Uz), 1 * τ11)
        @test_throws ArgumentError conjugate(τ12, Uz)
        # Constructor invariants. Both failures are silent otherwise: a non-normal amplitude
        # leaves `[A, A']` in the commutator of the images, and an amplitude on the displaced
        # site does not truncate the Hadamard series at all.
        @test occursin(
            "own adjoint", (
                try
                    Displace(acd, r * τ12)
                catch e
                    sprint(showerror, e)
                end
            )
        )
        @test_throws ArgumentError Displace(acd, r * (acd' * acd))
        # The phase-space form, with the Hermiticity its quadratures need.
        hxq = PhaseSpace(:xq) ⊗ FockSpace(:fq)
        xq, pq = Position(hxq, :xq, 1), Momentum(hxq, :pq, 1)
        aq = Destroy(hxq, :aq, 2)
        Uq = Displace(xq, pq, r * (aq' * aq), θ * (aq' * aq))
        @test is_canonical(Uq)
        @test isequal(conjugate(xq, Uq), xq + r * aq' * aq)
        @test isequal(conjugate(pq, Uq), pq + θ * aq' * aq)
        @test_throws ArgumentError Displace(xq, pq, r * (aq' * aq), 1 * aq)
    end

    @testset "generic exponentiation" begin
        # Symbolic equality is the wrong test between the two spellings: the generic
        # constructor writes a rotation on phase atoms and the hand-written one on
        # `cos`/`sin`, and closing `(p + p⁻¹)/2 = cos(θ)` needs Euler's formula, which no
        # `ParamRelation` states. They agree as functions, which is what a numeric
        # substitution tests.
        agrees(q1, q2, pars) = all(
            e -> abs(_to_complex(_substitute_cnum(e.second, pars))) < 1.0e-12,
            (q1 - q2).arguments,
        )
        pars = Dict(θ => 0.37, r => 0.23)
        # Every hand-written constructor that has a `(G, θ)` pair is reproduced by the
        # generic one. `Squeeze(a, r, ϕ)` is absent because a `ϕ`-dependent generator is not
        # readable off the rules, and `Rotation(σ, W)`/`Bogoliubov` have no generator at all.
        for (name, U, G, par) in (
                ("Fock rotation", Rotation(a, θ), a' * a, θ),
                ("beamsplitter", Rotation(ma, mb, θ), im * (ma' * mb - mb' * ma), θ),
                ("quadrature rotation", Rotation(x, p, θ), (x * x + p * p) * (1 // 2), θ),
                ("quadrature squeeze", Squeeze(x, p, r), (x * p + p * x) * (1 // 2), r),
                ("two-mode squeeze", Squeeze(ma, mb, r), im * (ma' * mb' - mb * ma), r),
                ("single-mode squeeze", Squeeze(a, r), im * (a' * a' - a * a) * (1 // 2), r),
                ("spin rotation", Rotation(Sx, 3, θ), 1 * Sz, θ),
                ("Pauli rotation", Rotation(σx, 3, θ), (1 // 2) * σz, θ),
            )
            @testset "$name" begin
                E = UnitaryTransform(G, par)
                @test is_canonical(E)
                for gen in generators(U)
                    @test agrees(conjugate(gen, E), conjugate(gen, U), pars)
                end
            end
        end
        # A quadrature squeeze is a Hermitian generator on an imaginary rate, so it is a real
        # scaling and not a phase. Both directions must land on one atom at opposite
        # exponents, or the residual never closes: `exp(-r)*exp(r)` does not fold.
        Eq = UnitaryTransform((x * p + p * x) * (1 // 2), r)
        @test isequal(conjugate(x, Eq), conjugate(x, Squeeze(x, p, r)))
        @test isequal(conjugate(p, Eq), conjugate(p, Squeeze(x, p, r)))
        @test is_canonical(Eq; atol = 0)
        # The two rates a diagonal generator may carry are exclusive, and each is wrong for
        # the other kind: a phase on a Hermitian generator breaks Hermiticity, and a real
        # scaling of a ladder operator scales its commutator by the square.
        @test occursin(
            "Hermitian", (
                try
                    UnitaryTransform(im * (x * p + p * x) * (1 // 2), θ)
                catch e
                    sprint(showerror, e)
                end
            )
        )
        @test_throws ArgumentError UnitaryTransform(im * ω * a' * a, θ)
        # The static form carries no gauge; the three-argument form carries `-∂ₜθ*G`.
        @test iszero(gauge_term(UnitaryTransform(a' * a, θ)))
        Uf = UnitaryTransform(a' * a, ω * sin(t), t)
        @test isequal(conjugate(a, Uf), expim(-ω * sin(t)) * a)
        @test isequal(gauge_term(Uf), -ω * cos(t) * a' * a)
        @test is_canonical(Uf; atol = 0)
        # `RotatingFrame` is the `θ = t` case, so `-∂ₜt*H0` is exactly the gauge it had.
        @test isequal(gauge_term(RotatingFrame(ω * a' * a, t)), -ω * a' * a)
    end

    @testset "a moving H0 is refused" begin
        # `exp(-im*H0*t)` is not the propagator of a moving `H0`. This used to return
        # `exp(-im*cos(t)*t*ω)*a` where the phase is `exp(-im*ω*sin(t))`, build the wrong
        # gauge, and certify: every part of that was silent.
        msg = try
            RotatingFrame(ω * cos(t) * a' * a, t)
            ""
        catch e
            sprint(showerror, e)
        end
        @test occursin("time-independent", msg)
        @test occursin("UnitaryTransform", msg)
        # The escape hatch it names gives the right answer.
        U = UnitaryTransform(a' * a, ω * sin(t), t)
        @test isequal(conjugate(a, U), expim(-ω * sin(t)) * a)
        @test iszero(transform(ω * cos(t) * a' * a, U))
        # A `t`-free `H0` is untouched by the guard.
        @test is_canonical(RotatingFrame(ω * a' * a, t))
    end

    @testset "driven frames" begin
        # `[H0, a] = -η` is a c-number, so the adjoint action is affine rather than linear.
        # A drive alone rotates nothing, and the shift is linear in the frame parameter.
        U0 = RotatingFrame(η * (a + a'), t)
        @test is_canonical(U0)
        @test isequal(conjugate(a, U0), conjugate(a, Displace(a, -im * η * t)))
        @test isequal(conjugate(a', U0), conjugate(a', Displace(a, -im * η * t)))
        @test isequal(gauge_term(U0), -η * (a + a'))
        @test iszero(transform(η * (a + a'), U0))
        # With a detuning the shift is `(exp(-im*ω*t) - 1)*η/ω`: the homogeneous phase acting
        # on the completed square, minus the square itself.
        Ud = RotatingFrame(ω * a' * a + η * (a + a'), t)
        @test is_canonical(Ud)
        @test isequal(
            conjugate(a, Ud),
            expim(-ω * t) * a + (expim(-ω * t) * η) / ω - η / ω,
        )
        @test iszero(simplify(transform(ω * a' * a + η * (a + a'), Ud)))
        # A driven quadrature pair stays Hermitian, which the CCR alone does not pin.
        Uq = RotatingFrame(ω * (x * x + p * p) * (1 // 2) + η * x, t)
        @test is_canonical(Uq)
        @test iszero(simplify(transform(ω * (x * x + p * p) * (1 // 2) + η * x, Uq)))
        # A resonantly driven dark mode has no closed form here: the drive is outside the
        # range of the homogeneous part, so the frame grows without bound.
        h3c = FockSpace(:c1) ⊗ FockSpace(:c2) ⊗ FockSpace(:c3)
        c1, c2, c3 = Destroy(h3c, :c1, 1), Destroy(h3c, :c2, 2), Destroy(h3c, :c3, 3)
        chain = g * (c1' * c2 + c2' * c1 + c2' * c3 + c3' * c2)
        @test is_canonical(RotatingFrame(chain, t))
        @test occursin(
            "resonantly", (
                try
                    RotatingFrame(chain + η * (c1 - c3 + c1' - c3'), t)
                catch e
                    sprint(showerror, e)
                end
            )
        )
    end

    @testset "involutive blocks" begin
        # A detuned parametric amplifier has a nonzero, unshared diagonal, so no `2x2` form
        # reads off it, and its block is not Hermitian, so `eigen` will not take it either.
        # `C^2 = κ*I` closes it anyway, with the rate still symbolic in the general case.
        H0 = 1.0 * a' * a + 0.3 * (a * a + a' * a')
        U = RotatingFrame(H0, t)
        @test is_canonical(U)
        @test SecondQuantizedAlgebra._residual_iszero(transform(H0, U), 1.0e-10)
        # Every regime of the block closes, on both sides of threshold and on it. The root is
        # left unevaluated rather than taken, so `κ < 0` is a symbolic `√(-35)` instead of a
        # rounded `cosh` weight, and all three certify with no tolerance at all.
        for (lbl, Δv, gv) in (
                ("below", 2, 1 // 2), ("threshold", 2, 1), ("above", 1, 3),
            )
            Uk = RotatingFrame(Δv * a' * a + gv * (a * a + a' * a'), t)
            @test is_canonical(Uk; atol = 0.0)
            @test SecondQuantizedAlgebra._residual_iszero(
                transform(Δv * a' * a + gv * (a * a + a' * a'), Uk), 0.0
            )
        end
        # Symbolic couplings, where the sign of `κ = Δ^2 - 4g^2` is not even decidable. One
        # `cos`/`sin` form covers both branches, because both are even in `√κ`.
        @test is_canonical(
            RotatingFrame(Δ * a' * a + g * (a * a + a' * a'), t); atol = 0.0
        )
        # At threshold the block is nilpotent, so the series stops after one term and the
        # frame is a polynomial in `t` carrying no atom at all.
        @test isequal(
            conjugate(a, RotatingFrame(2 * a' * a + 1 * (a * a + a' * a'), t)),
            (1 - 2im * t) * a - 2im * t * a',
        )
        # The exact `2x2` path still goes first, so the forms it produces do not move. The
        # beamsplitter, harmonic and parametric spellings are pinned in "rotating frames";
        # the two-mode squeeze is only pinned here.
        @test isequal(
            conjugate(ma, RotatingFrame(g * (ma' * mb' + mb * ma), t)),
            cosh(g * t) * ma - im * sinh(g * t) * mb',
        )
    end

    @testset "rotating frames" begin
        # The per-site gauge terms sum to exactly `-H0`.
        for H0 in (
                ω * a' * a,
                ω * Sz,
                Δ * σ22,
                ω * ma' * ma + Ω * mb' * mb,
            )
            U = RotatingFrame(H0, t)
            @test isequal(gauge_term(U), -H0)
            @test iszero(transform(H0, U))
        end
        # Off-diagonal populations of a multilevel frame.
        h3 = NLevelSpace(:three, 3)
        t11 = Transition(h3, :τ, 1, 1)
        t22 = Transition(h3, :τ, 2, 2)
        H3 = ω * t11 + Δ * t22
        @test iszero(transform(H3, RotatingFrame(H3, t)))
        # A constant in `H0` is a global phase and lands in the gauge.
        Uc = RotatingFrame(ω * a' * a + 5, t)
        @test isequal(gauge_term(Uc), -ω * a' * a - 5)
        @test iszero(transform(ω * a' * a + 5, Uc))
        # The driven mode: the frame kills the number term and phases the drive.
        Hd = transform(ω * a' * a + η * (a + a'), RotatingFrame(ω * a' * a, t))
        @test length(Hd) == 2
        # The drive carries one phase atom per rate, so the two terms are mutual conjugates
        # at opposite exponents rather than four unrelated trigonometric factors.
        @test isequal(Hd, η * expim(-ω * t) * a + η * expim(ω * t) * a')
        # A symbolic level-diagonal frame is exact: one phase per level, so every residual
        # closes by exponent arithmetic with nothing left to bound.
        hsym = NLevelSpace(:sym3, 3)
        ς(i, j) = Transition(hsym, :ς, i, j)
        Symbolics.@variables Ea Eb Ec
        Usym = RotatingFrame(Ea * ς(1, 1) + Eb * ς(2, 2) + Ec * ς(3, 3), t)
        @test all(iszero, canonicality_residuals(Usym))
        @test isequal(conjugate(ς(1, 2), Usym), expim(Ea * t) * expim(-Eb * t) * ς(1, 2))
        # Opposite rates on one frequency share an atom, so the phases cancel outright.
        hsi = FockSpace(:si)
        asi, bsi = Destroy(hsi, :as), Destroy(hsi, :bs)
        Usi = RotatingFrame(ω * asi' * asi - ω * bsi' * bsi, t)
        @test isequal(conjugate(asi * bsi, Usi), 1 * (asi * bsi))
        # Mixing frames: the map comes from the `2x2` block of `ad_{H0}`, so a harmonic, a
        # beamsplitter and a parametric `H0` all work without being special-cased.
        # A rotating block is two phase atoms at opposite exponents; `_c`/`_s` spell the
        # `cos`/`sin` it amounts to. Scaling the operator rather than the bare phase keeps
        # the arithmetic in the coefficient algebra, where the phases stay one factor. The
        # parametric block is hyperbolic, not a phase, and keeps its `cosh`/`sinh`.
        _c(θ, o) = 0.5 * (expim(θ) * o + expim(-θ) * o)
        _s(θ, o) = -0.5im * (expim(θ) * o - expim(-θ) * o)
        for (H0, probe, want) in (
                (ω * (x * x + p * p) / 2, x, _c(ω * t, x) + _s(ω * t, p)),
                (
                    g * (ma' * mb + mb' * ma), ma,
                    _c(g * t, ma) - im * _s(g * t, mb),
                ),
                (
                    g * (ma * ma + ma' * ma'), ma,
                    cosh(2 * g * t) * ma - im * sinh(2 * g * t) * ma',
                ),
                (ω * Sz, Sx, _c(ω * t, Sx) - _s(ω * t, Sy)),
            )
            U = RotatingFrame(H0, t)
            @test isequal(conjugate(probe, U), want)
            @test is_canonical(U)
            @test isequal(gauge_term(U), -H0)
            @test iszero(transform(H0, U))
        end
        # An off-diagonal two-level `H0` factors through the dressed basis.
        for (Δv, Ωv) in ((1.0, 0.4), (-0.6, 1.3))
            H0 = (Δv / 2) * σ11 - (Δv / 2) * σ22 + Ωv * (σ12 + σ12')
            U = RotatingFrame(H0, t)
            @test is_canonical(U)
            @test SecondQuantizedAlgebra._residual_iszero(gauge_term(U) + H0, 1.0e-12)
            @test SecondQuantizedAlgebra._residual_iszero(transform(H0, U), 1.0e-12)
        end
        # Any number of coupled generators, as long as the block is one symbolic rate away
        # from numeric. One phase atom per eigenvalue is what makes these close: the trig
        # spelling would need the angle-addition identity between the eigenvalues.
        h3f = FockSpace(:c1) ⊗ FockSpace(:c2) ⊗ FockSpace(:c3)
        m1 = Destroy(h3f, :m1, 1)
        m2 = Destroy(h3f, :m2, 2)
        m3 = Destroy(h3f, :m3, 3)
        nn = m1' * m1 + m2' * m2 + m3' * m3
        chain = g * (m1' * m2 + m2' * m1 + m2' * m3 + m3' * m2)
        ring = chain + g * (m1' * m3 + m3' * m1)
        for H0 in (chain, ω * nn + chain, ω * nn + ring)
            U = RotatingFrame(H0, t)
            @test is_canonical(U)
            # The frame mixes modes, so cross-mode commutativity is a live condition and the
            # residual list has to carry it.
            @test SecondQuantizedAlgebra._mixes_sites(U)
            @test isequal(gauge_term(U), -H0)
        end
        # A uniform diagonal comes out as its own shared phase, so the mode frequency cancels
        # from a number-conserving observable instead of riding inside every eigenphase.
        @test SecondQuantizedAlgebra._residual_iszero(
            conjugate(nn, RotatingFrame(ω * nn + chain, t)) - nn, 1.0e-12
        )
        # Two modes with a nonzero diagonal: the exact `2x2` form does not apply (it needs a
        # traceless block), so this goes through the general path as well.
        U2c = RotatingFrame(ω * (ma' * ma + mb' * mb) + g * (ma' * mb + mb' * ma), t)
        @test is_canonical(U2c)
        # The atoms are exact; the weights in front of them come from `eigen`, so `0.5` lands
        # a rounding away and the comparison has to be the bounded one.
        @test SecondQuantizedAlgebra._residual_iszero(
            conjugate(ma, U2c) -
                (
                0.5 * (expim(g * t) + expim(-g * t)) * expim(-ω * t) * ma -
                    0.5 * (expim(g * t) - expim(-g * t)) * expim(-ω * t) * mb
            ),
            1.0e-12,
        )
        # A spin rotation about an arbitrary axis couples all three axes at once.
        @test is_canonical(RotatingFrame(ω * (Sx + Sz), t))
        @test is_canonical(RotatingFrame(ω * (σx + σy + σz), t))
        # Cross-site commutativity has teeth: flipping the sign of one cross term of the
        # two-mode squeeze leaves both same-site commutators right and breaks `[ã, b̃]`, so
        # only a residual across sites can reject it.
        bad2 = Dict{Op, QAdd}(
            ma => cosh(r) * ma + sinh(r) * mb', ma' => cosh(r) * ma' + sinh(r) * mb,
            mb => cosh(r) * mb - sinh(r) * ma', mb' => cosh(r) * mb' - sinh(r) * ma,
        )
        Ux = UnitaryTransform(
            bad2, bad2, _zero_qadd(), Num(0),
            ParamRelation[ParamRelation(cosh(r), sinh(r), 1)], Num[],
        )
        @test iszero(simplify(commutator(conjugate(ma, Ux), conjugate(ma', Ux)) - 1))
        @test iszero(simplify(commutator(conjugate(mb, Ux), conjugate(mb', Ux)) - 1))
        @test !iszero(simplify(commutator(conjugate(ma, Ux), conjugate(mb, Ux))))
        @test !is_canonical(Ux)
        # Unsupported frames throw with a specific message.
        @test_throws ArgumentError RotatingFrame(Ω * σ12, t)
        @test_throws ArgumentError RotatingFrame(_single_qadd(_CNUM_ONE, Op[]), t)
        @test_throws ArgumentError RotatingFrame(im * ω * a' * a, t)
        # Off the diagonal the eigenbasis is numeric, so symbolic couplings still throw.
        @test_throws ArgumentError RotatingFrame(Δ * σ11 + Ω * (σ12 + σ12'), t)
        # Beyond two levels the frame frequencies are commensurate. One phase per level
        # makes that exponent arithmetic, so the frame builds and certifies.
        h3o = NLevelSpace(:three_o, 3)
        τ(i, j) = Transition(h3o, :τ, i, j)
        H3o = 1.0 * τ(2, 2) + 2.5 * τ(3, 3) + 0.3 * (τ(1, 2) + τ(2, 1))
        U3o = RotatingFrame(H3o, t)
        @test is_canonical(U3o)
        @test SecondQuantizedAlgebra._residual_iszero(transform(H3o, U3o), 1.0e-12)
        # A four-level ladder: every difference is a sum of the ones below it.
        h4o = NLevelSpace(:four_o, 4)
        υ(i, j) = Transition(h4o, :υ, i, j)
        H4o = 1.0 * υ(2, 2) + 2.3 * υ(3, 3) + 3.9 * υ(4, 4) +
            0.4 * (υ(1, 2) + υ(2, 1)) + 0.7 * (υ(3, 4) + υ(4, 3))
        @test is_canonical(RotatingFrame(H4o, t))
        # A mixed atom-cavity `H0` leaves the linear span of the generators; the message
        # names the product that does it.
        hjc = FockSpace(:c) ⊗ NLevelSpace(:jc, 2)
        ajc = Destroy(hjc, :a, 1)
        sjc(i, j) = Transition(hjc, :σ, i, j, 2)
        msg = try
            RotatingFrame(Δ * sjc(2, 2) + g * (ajc' * sjc(1, 2) + ajc * sjc(2, 1)), t)
            ""
        catch e
            sprint(showerror, e)
        end
        @test occursin("exactly solvable", msg) && occursin("σ", msg)
    end

    @testset "dressed frames" begin
        h2 = NLevelSpace(:dressed, 2)
        sg(i, j) = Transition(h2, :σ, i, j)
        for (Δv, Ωv) in ((1.0, 0.4), (0.3, 2.0), (-1.5, 0.7))
            H = (Δv / 2) * sg(1, 1) - (Δv / 2) * sg(2, 2) + Ωv * (sg(1, 2) + sg(2, 1))
            U = DressedFrame(H)
            @test is_canonical(U)
            @test iszero(gauge_term(U))
            D = conjugate(H, U)
            off = [_to_complex(c) for (tm, c) in D if tm.ops[1].l1 != tm.ops[1].l2]
            dia = sort([real(_to_complex(c)) for (tm, c) in D if tm.ops[1].l1 == tm.ops[1].l2])
            @test maximum(abs, off; init = 0.0) < 1.0e-12
            @test dia ≈ [-sqrt(Δv^2 / 4 + Ωv^2), sqrt(Δv^2 / 4 + Ωv^2)]
        end
        # A symbolic frame is still exact even though the off-diagonal is not reduced away.
        @test is_canonical(
            DressedFrame((Δ / 2) * sg(1, 1) - (Δ / 2) * sg(2, 2) + Ω * (sg(1, 2) + sg(2, 1)))
        )
        # Numeric level energies diagonalize at any number of levels.
        h3d = NLevelSpace(:three_d, 3)
        td(i, j) = Transition(h3d, :τ, i, j)
        H3d = 1.0 * td(2, 2) + 2.5 * td(3, 3) + 0.3 * (td(1, 2) + td(2, 1)) +
            0.2 * (td(2, 3) + td(3, 2))
        U3 = DressedFrame(H3d)
        @test is_canonical(U3)
        D3 = conjugate(H3d, U3)
        @test maximum(
            abs, [_to_complex(c) for (tm, c) in D3 if tm.ops[1].l1 != tm.ops[1].l2];
            init = 0.0,
        ) < 1.0e-12
        # A symbolic `H0` diagonalizes block by block, so more than two levels is fine as
        # long as no block of three couples: here levels 1-2 mix, 3 is a spectator.
        Ub = DressedFrame(Δ * td(1, 1) + Ω * (td(1, 2) + td(2, 1)) + η * td(3, 3))
        @test is_canonical(Ub)
        # Three or more coupled levels diagonalize symbolically when the block is one real
        # rate away from a numeric Hermitian matrix: the eigenvectors are then numbers, which
        # is all a dressed basis is, and the energies `d + w*μ` stay symbolic. A uniformly
        # coupled chain qualifies with both its frequency and its coupling free.
        Hchain = ω * (td(1, 1) + td(2, 2) + td(3, 3)) +
            g * (td(1, 2) + td(2, 1) + td(2, 3) + td(3, 2))
        Uchain = DressedFrame(Hchain)
        @test is_canonical(Uchain)
        # The dressed matrix is diagonal. Its off-diagonal residue is rounding scaled by the
        # parameters, which no coefficient bound can wave through, so it is read at a point.
        chvals = Dict(ω => 0.7, g => 0.35)
        @test all(
            e -> only(e.first.ops).l1 == only(e.first.ops).l2 ||
                abs(_to_complex(_substitute_cnum(e.second, chvals))) < 1.0e-10,
            conjugate(Hchain, Uchain).arguments,
        )
        # The frame of the same `H0` follows, through the dressed factorization.
        @test is_canonical(RotatingFrame(Hchain, t))
        # A symbolic two-level block cannot supply dressed energies (they need a `sqrt` no
        # relation closes), so it stays a `DressedFrame` and is refused as a frame.
        Hs2 = (Δ / 2) * σ11 - (Δ / 2) * σ22 + Ω * (σ12 + σ12')
        @test is_canonical(DressedFrame(Hs2))
        @test occursin(
            "sqrt", (
                try
                    RotatingFrame(Hs2, t)
                catch e
                    sprint(showerror, e)
                end
            )
        )
        # A genuinely coupled symbolic three-level block with two independent couplings has
        # no such split, and is refused naming the levels.
        @test occursin(
            "levels 1, 2, 3", (
                try
                    DressedFrame(
                        Δ * td(1, 1) + Ω * (td(1, 2) + td(2, 1)) + η * (td(2, 3) + td(3, 2))
                    )
                catch e
                    sprint(showerror, e)
                end
            )
        )
        # Non-Hermitian input is refused.
        @test_throws ArgumentError DressedFrame(ω * a' * a)
        @test_throws ArgumentError DressedFrame(
            Ω * sg(1, 2) + 2 * Ω * sg(2, 1)           # not Hermitian
        )
        @test_throws ArgumentError DressedFrame(
            1.0 * sg(1, 2) + 2.0 * sg(2, 1)           # not Hermitian, numeric path
        )
    end

    @testset "silent-failure guards" begin
        # The level-diagonal N-level frame took `real(M[k,k])` unchecked, so a non-Hermitian
        # `H0` built fine, reported itself canonical, and carried a non-Hermitian gauge.
        @test_throws ArgumentError RotatingFrame(im * Δ * σ22, t)
        @test is_canonical(RotatingFrame(Δ * σ22, t))

        # A static transform contributes a zero gauge. Adopting a time variable its own rules
        # depend on makes that zero wrong rather than absent.
        @variables Amp
        @test_throws ArgumentError Displace(a, Amp * t) * RotatingFrame(ω * a' * a, t)
        @test Displace(a, Amp * t, t) * RotatingFrame(ω * a' * a, t) isa UnitaryTransform
        @test Displace(a, Amp) * RotatingFrame(ω * a' * a, t) isa UnitaryTransform
    end

    @testset "conjugate against a numeric U'AU oracle" begin
        # `is_canonical` certifies a transform through the same coefficient algebra that
        # produced its residuals, and the only independent oracle in this file was pointed at
        # three gauge terms. This one builds `U` as an actual matrix exponential and compares
        # `conjugate(A, U)` with `U'AU`, through the package's own `to_numeric`.
        N = 70
        keep = 1:25                       # displacing/squeezing spreads up the truncated ladder
        b = FockBasis(N - 1)
        mat(q) = Matrix(to_numeric(q, b).data)
        sub(M) = M[keep, keep]
        Am, Adm = mat(a), mat(a')

        for (name, U, G) in (
                # `U` built as a matrix exponential of its own generator, independent of
                # everything the transform layer computed.
                ("Displace", Displace(a, 0.35), 0.35 * (Adm - Am)),
                ("Rotation", Rotation(a, 0.7), -im * 0.7 * (Adm * Am)),
                ("Squeeze", Squeeze(a, 0.22), 0.22 * (Adm * Adm - Am * Am) / 2),
            )
            @testset "$name" begin
                Um = exp(G)
                for op in (a, a', a' * a)
                    @test norm(sub(mat(conjugate(op, U))) - sub(Um' * mat(op) * Um)) < 1.0e-8
                end
            end
        end

        # A rotating frame is where the phase atom has to survive all the way to a number.
        @testset "RotatingFrame at a concrete time" begin
            τ = 0.37
            Hd = transform(ω * a' * a + η * (a + a'), RotatingFrame(ω * a' * a, t))
            got = mat(substitute(Hd, Dict(ω => 1.3, η => 0.4, t => τ)))
            Um = exp(-im * 1.3 * τ * (Adm * Am))
            want = Um' * mat(0.4 * (a + a')) * Um
            @test norm(sub(got) - sub(want)) < 1.0e-8
        end

        # The operator-valued displacement gauge, against the same central-difference oracle
        # the scalar ones use. Here the piece that is a global phase for a number is an
        # operator, so it is not merely load-bearing but visible to the dynamics.
        @testset "operator-valued Displace gauge" begin
            hg = FockSpace(:cg) ⊗ FockSpace(:mg)
            acg, bmg = Destroy(hg, :acg, 1), Destroy(hg, :bmg, 2)
            Symbolics.@variables gq wq
            bg = FockBasis(7) ⊗ FockBasis(13)
            mg(q) = Matrix(to_numeric(q, bg).data)
            gv, wv, τ0, dτ = 0.37, 1.9, 0.41, 1.0e-5
            Amat = (-gv / wv) * mg(acg' * acg)
            Uτ(τ) = exp(τ * Amat * mg(bmg') - τ * Amat * mg(bmg))
            ref = 1im * ((Uτ(τ0 + dτ)' - Uτ(τ0 - dτ)') / (2dτ)) * Uτ(τ0)
            Ug = Displace(bmg, (-gq * t / wq) * (acg' * acg), t)
            got = mg(substitute(gauge_term(Ug), Dict(gq => gv, wq => wv, t => τ0)))
            kg = 1:64
            @test norm((got - ref)[kg, kg]) < 1.0e-7
            @test norm(ref[kg, kg]) > 1.0     # the gauge is not zero to begin with
        end

        # A symbolic dressed frame is built on a numeric eigenbasis, so `transform(H0, U)`
        # carries rounding scaled by the parameters of `H0` and no coefficient bound can
        # certify it. An oracle can.
        @testset "symbolic dressed frame at a point" begin
            hch = NLevelSpace(:chain, 3)
            τc(i, j) = Transition(hch, :τc, i, j)
            H0 = ω * (τc(1, 1) + τc(2, 2) + τc(3, 3)) +
                g * (τc(1, 2) + τc(2, 1) + τc(2, 3) + τc(3, 2))
            U = RotatingFrame(H0, t)
            vals = Dict(ω => 0.7, g => 0.35, t => 0.4)
            bch = NLevelBasis(3)
            mch(q) = Matrix(to_numeric(q, bch).data)
            Um = exp(-im * 0.4 * mch(substitute(H0, vals)))
            for i in 1:3, j in 1:3
                got = mch(substitute(conjugate(τc(i, j), U), vals))
                @test norm(got - Um' * mch(τc(i, j)) * Um) < 1.0e-10
            end
            @test norm(mch(substitute(transform(H0, U), vals))) < 1.0e-10
        end

        # An involutive block mixes `a` with `a'`, so it spreads up the ladder like a squeeze
        # and needs the deeper truncation.
        @testset "involutive block at a concrete time" begin
            Nin, kin, τ = 160, 1:40, 0.3
            bin = FockBasis(Nin - 1)
            min(q) = Matrix(to_numeric(q, bin).data)
            H0 = 1.0 * a' * a + 0.3 * (a * a + a' * a')
            Um = exp(-im * τ * min(H0))
            U = RotatingFrame(H0, t)
            for op in (a, a')
                got = min(substitute(conjugate(op, U), Dict(t => τ)))
                @test norm((got - Um' * min(op) * Um)[kin, kin]) < 1.0e-8
            end
        end

        # The affine shift is the one part of a frame `is_canonical` cannot pin: a constant
        # commutes out of the CCR residual, and the round trip only catches a forward and
        # inverse that disagree, not a shift wrong by the same factor in both. Only an
        # oracle fixes the magnitude.
        @testset "driven frame at a concrete time" begin
            τ, ωv, ηv = 0.37, 1.3, 0.4
            for (name, H0, U) in (
                    ("pure drive", ηv * (a + a'), RotatingFrame(η * (a + a'), t)),
                    (
                        "detuned drive", ωv * a' * a + ηv * (a + a'),
                        RotatingFrame(ω * a' * a + η * (a + a'), t),
                    ),
                )
                @testset "$name" begin
                    Um = exp(-im * τ * mat(H0))
                    for op in (a, a', a' * a)
                        got = mat(substitute(conjugate(op, U), Dict(ω => ωv, η => ηv, t => τ)))
                        @test norm(sub(got) - sub(Um' * mat(op) * Um)) < 1.0e-8
                    end
                end
            end
        end
    end

    @testset "conjugate is a *-algebra homomorphism" begin
        # A transform is stored rule by rule on the generators. That it extends to every sum
        # and product is what makes applying one to a Hamiltonian mean anything, and only the
        # closed forms on single generators were pinned above.
        # A bare product of two images never sees the transform's own reductions, so
        # `cosh(r)^2 - sinh(r)^2 = 1` cannot close symbolically and the comparison is at a
        # point rather than by `isequal`.
        agrees(q1, q2, pars) = all(
            e -> abs(_to_complex(_substitute_cnum(e.second, pars))) < 1.0e-12,
            (q1 - q2).arguments,
        )
        pars = Dict(
            θ => 0.37, r => 0.23, ϕ => 0.61, dx => 0.4, dp => -0.15,
            u => cosh(0.23), v => sinh(0.23),
        )
        for (name, U) in (
                "Rotation" => Rotation(a, θ),
                "Displace" => Displace(a, u),
                "Squeeze" => Squeeze(a, r),
                "Squeeze(r,ϕ)" => Squeeze(a, r, ϕ),
                "beamsplitter" => Rotation(ma, mb, θ),
                "two-mode squeeze" => Squeeze(ma, mb, r),
                "Bogoliubov" => Bogoliubov(ma, mb, u, v),
                "phase Rotation" => Rotation(x, p, θ),
                "phase Squeeze" => Squeeze(x, p, r),
                "phase Displace" => Displace(x, p, dx, dp),
                "spin Rotation" => Rotation(Sx, 3, θ),
                "NLevel Rotation" => Rotation(σ12, [cos(θ) -sin(θ); sin(θ) cos(θ)]),
            )
            @testset "$name" begin
                gs = generators(U)
                for g1 in gs, g2 in gs
                    @test agrees(
                        conjugate(g1 * g2, U), conjugate(g1, U) * conjugate(g2, U), pars
                    )
                    @test agrees(
                        conjugate(g1 + 2 * g2, U),
                        conjugate(g1, U) + 2 * conjugate(g2, U), pars,
                    )
                end
                for g in gs
                    @test agrees(conjugate(adjoint(g), U), adjoint(conjugate(g, U)), pars)
                end
                # A c-number is fixed by conjugation, including the zero and the unit.
                @test iszero(conjugate(_zero_qadd(), U))
                @test isequal(
                    conjugate(_single_qadd(_CNUM_ONE, Op[]), U),
                    _single_qadd(_CNUM_ONE, Op[]),
                )
            end
        end
    end

    @testset "composition and inversion laws" begin
        U1, U2, U3 = Rotation(a, θ), Squeeze(a, r, ϕ), Displace(a, u)
        # Composition is associative, and inversion reverses it. Both are properties of
        # `_merge_rules`, whose two calls run in opposite orders.
        for g in generators(U1)
            @test isequal(
                simplify(conjugate(g, (U1 * U2) * U3)),
                simplify(conjugate(g, U1 * (U2 * U3))),
            )
            @test isequal(
                simplify(conjugate(g, inv(U1 * U2))),
                simplify(conjugate(g, inv(U2) * inv(U1))),
            )
            @test isequal(simplify(conjugate(g, U1 * U2 * inv(U2) * inv(U1))), 1 * g)
        end
        # Constraints and blockers are properties of the factors and have to survive both.
        Ub = Bogoliubov(ma, mb, u, v)
        @test isequal(constraints(Ub * Rotation(ma, θ)), constraints(Ub))
        @test isequal(constraints(Rotation(ma, θ) * Ub), constraints(Ub))
        @test isequal(constraints(inv(Ub)), constraints(Ub))
        hpol = FockSpace(:cav2) ⊗ FockSpace(:mech2)
        ac, bm = Destroy(hpol, :ac2, 1), Destroy(hpol, :bm2, 2)
        Up = Displace(bm, (-g / ω) * (ac' * ac))
        for Uc in (Up * Rotation(bm, θ), Rotation(bm, θ) * Up, inv(Up))
            @test length(Uc.blockers) == 1
            @test isequal(conjugate(ac' * ac, Uc), ac' * ac)
            @test_throws ArgumentError conjugate(ac, Uc)
        end
    end

    @testset "a frame and its inverse cancel" begin
        # `transform` is not `conjugate`: undoing it has to undo the gauge as well, and
        # `inv` routes that through the inverse rules rather than negating it. Only a
        # `transform`-level round trip tests that, and for an operator-valued gauge the two
        # spellings genuinely differ.
        H = ω * a' * a + η * (a + a')
        for (name, U) in (
                "Rotation(ωt)" => Rotation(a, ω * t, t),
                "RotatingFrame" => RotatingFrame(ω * a' * a, t),
                "Displace(t)" => Displace(a, u * t, t),
                "Squeeze(ϕ(t))" => Squeeze(a, r, ω * t, t),
                "generic" => UnitaryTransform(a' * a, ω * sin(t), t),
            )
            @testset "$name" begin
                @test iszero(simplify(transform(transform(H, U), inv(U)) - H))
            end
        end
        # The polaron gauge is an operator, so `-U*gauge(U)*U'` does not collapse to
        # `-gauge(U)` and this is the case that pins the routing.
        hpg = FockSpace(:cavp) ⊗ FockSpace(:mechp)
        acp, bmp = Destroy(hpg, :acp, 1), Destroy(hpg, :bmp, 2)
        Ug = Displace(bmp, (-g * t / ω) * (acp' * acp), t)
        Hp = ω * bmp' * bmp + g * acp' * acp * (bmp + bmp') + Δ * acp' * acp
        @test iszero(simplify(transform(transform(Hp, Ug), inv(Ug)) - Hp))
        # Two moving factors: the composed gauge has to compose too.
        Ut1, Ut2 = Rotation(a, ω * t, t), Squeeze(a, r, ω * t, t)
        @test iszero(
            simplify(transform(H, Ut1 * Ut2) - transform(transform(H, Ut1), Ut2))
        )
        @test iszero(
            simplify(gauge_term(inv(Ut1 * Ut2)) - gauge_term(inv(Ut2) * inv(Ut1)))
        )
    end

    @testset "transform preserves Hermiticity" begin
        # A gauge term that is not Hermitian generates non-unitary motion. Nothing above
        # tests the sum, only the closed form of the gauge on its own.
        H = ω * a' * a + η * (a + a')
        for U in (
                Rotation(a, θ), Squeeze(a, r, ϕ), Displace(a, u),
                Rotation(a, ω * t, t), RotatingFrame(ω * a' * a, t),
                Squeeze(a, r, ω * t, t), UnitaryTransform(a' * a, ω * sin(t), t),
            )
            T = transform(H, U)
            @test iszero(simplify(T - adjoint(T)))
        end
    end

    @testset "spin and Pauli rotations against a matrix oracle" begin
        # The three-axis commutator is closed under any representation, so the same symbolic
        # rule has to reproduce `exp(-im*θ*S)` at every spin. A rule that silently assumed
        # spin-1/2 would pass at `1//2` and nowhere else.
        hsq = SpinSpace(:Sq)
        Sq(k) = Spin(hsq, :Sq, k)
        for s in (1 // 2, 1 // 1, 3 // 2)
            @testset "spin $s" begin
                bsp = SpinBasis(s)
                msp(q) = Matrix(to_numeric(q, bsp).data)
                for ax in 1:3
                    Um = exp(-im * 0.7 * msp(Sq(ax)))
                    Ur = Rotation(Sq(ax), ax, 0.7)
                    for k in 1:3
                        @test norm(
                            msp(conjugate(Sq(k), Ur)) - Um' * msp(Sq(k)) * Um
                        ) < 1.0e-12
                    end
                end
            end
        end
        # A Pauli rotation is generated by `σ/2`, which is where the factor of two that the
        # gauge term carries has to show up in the map as well.
        bpa = SpinBasis(1 // 2)
        mpa(q) = Matrix(to_numeric(q, bpa).data)
        for ax in 1:3
            Um = exp(-im * 0.7 / 2 * mpa(Pauli(hpa, :σ, ax)))
            Ur = Rotation(Pauli(hpa, :σ, ax), ax, 0.7)
            for k in 1:3
                pk = Pauli(hpa, :σ, k)
                @test norm(mpa(conjugate(pk, Ur)) - Um' * mpa(pk) * Um) < 1.0e-12
            end
        end
    end

    @testset "two-mode transforms against a matrix oracle" begin
        # The two-mode maps were pinned by their closed forms and by `is_canonical`, which
        # both run inside the same coefficient algebra. This builds `U` as a matrix
        # exponential of the generator instead.
        Nm = 14
        bm2 = FockBasis(Nm - 1) ⊗ FockBasis(Nm - 1)
        m2(q) = Matrix(to_numeric(q, bm2).data)
        # A two-mode squeeze moves quanta between the modes, so the comparison is read on the
        # block where neither word reaches the truncation.
        na, nb = real.(diag(m2(ma' * ma))), real.(diag(m2(mb' * mb)))
        kp = findall(i -> na[i] <= 3.5 && nb[i] <= 3.5, eachindex(na))
        for (name, U, G, val) in (
                ("beamsplitter", Rotation(ma, mb, 0.7), im * (ma' * mb - mb' * ma), 0.7),
                ("two-mode squeeze", Squeeze(ma, mb, 0.22), im * (ma' * mb' - mb * ma), 0.22),
            )
            @testset "$name" begin
                Um = exp(-im * val * m2(G))
                for op in (ma, ma', mb, mb', ma' * ma, ma' * mb)
                    @test norm(
                        (m2(conjugate(op, U)) - Um' * m2(op) * Um)[kp, kp]
                    ) < 1.0e-6
                end
            end
        end
    end

    @testset "indexed families: inverse and composition" begin
        # `U = ⊗ᵢUᵢ` has to survive both operations family-wide, not just at the wildcard the
        # rules are keyed on.
        hfi = FockSpace(:fi)
        afi = Destroy(hfi, :afi)
        ki, li = Index(hfi, :ki, 4, hfi), Index(hfi, :li, 4, hfi)
        ak, al = IndexedOperator(afi, ki), IndexedOperator(afi, li)
        Ur, Us = Rotation(ak, θ), Squeeze(ak, r)
        @test isequal(simplify(conjugate(conjugate(al, Ur), inv(Ur))), 1 * al)
        Uc = Ur * Us
        @test is_canonical(Uc)
        @test isequal(
            simplify(conjugate(al, Uc)), simplify(conjugate(conjugate(al, Ur), Us))
        )
        @test iszero(simplify(conjugate(al' * al, Us) - adjoint(conjugate(al' * al, Us))))
    end

    @testset "closed forms reproduce the Hadamard series" begin
        # Every closed form claims to be `exp(im*θ*G) X exp(-im*θ*G)` for a stated `G`, but
        # each was derived and then written down; the generator is never used again. A central
        # difference in the parameter recovers the first two terms of the series from the
        # rules alone, and compares them against the package's own commutator. Nothing else
        # here ties a rule map to its generator at all.
        maxc(q) = maximum(e -> abs(_to_complex(e.second)), q.arguments; init = 0.0)
        δ = 1.0e-5
        for (name, U, G, probes) in (
                ("Rotation", θv -> Rotation(a, θv), a' * a, (a, a', a' * a)),
                ("Displace", αv -> Displace(a, αv), im * (a' - a), (a, a', a' * a)),
                (
                    "Squeeze", rv -> Squeeze(a, rv),
                    im * (a' * a' - a * a) * (1 // 2), (a, a'),
                ),
                (
                    "beamsplitter", θv -> Rotation(ma, mb, θv),
                    im * (ma' * mb - mb' * ma), (ma, mb', ma' * ma),
                ),
                (
                    "quadrature squeeze", rv -> Squeeze(x, p, rv),
                    (x * p + p * x) * (1 // 2), (x, p),
                ),
                ("spin rotation", θv -> Rotation(Sx, 3, θv), 1 * Sz, (Sx, Sy, Sz)),
                ("Pauli rotation", θv -> Rotation(σx, 3, θv), (1 // 2) * σz, (σx, σy, σz)),
            )
            @testset "$name" begin
                for A in probes
                    d1 = (conjugate(A, U(δ)) - conjugate(A, U(-δ))) * (1 / (2δ))
                    @test maxc(d1 - im * commutator(G, A)) < 1.0e-8
                    d2 = (
                        conjugate(A, U(δ)) - 2 * conjugate(A, U(0.0)) +
                            conjugate(A, U(-δ))
                    ) * (1 / δ^2)
                    @test maxc(d2 + commutator(G, commutator(G, A))) < 1.0e-5
                end
            end
        end
    end

    @testset "one-parameter group law" begin
        # `U(θ₁)U(θ₂) = U(θ₁ + θ₂)` for every constructor with an additive parameter. The
        # composition path merges rules by substitution, so this is a statement about
        # `_merge_rules` on a site both factors cover, which nothing else exercises.
        @variables θ1 θ2
        pars = Dict(θ1 => 0.3, θ2 => 0.45)
        agrees(q1, q2) = all(
            e -> abs(_to_complex(_substitute_cnum(e.second, pars))) < 1.0e-12,
            (q1 - q2).arguments,
        )
        for (name, f, gens) in (
                ("Rotation", s -> Rotation(a, s), (a, a')),
                ("Displace", s -> Displace(a, s), (a, a')),
                ("Squeeze", s -> Squeeze(a, s), (a, a')),
                ("beamsplitter", s -> Rotation(ma, mb, s), (ma, ma', mb, mb')),
                ("two-mode squeeze", s -> Squeeze(ma, mb, s), (ma, ma', mb, mb')),
                ("phase Rotation", s -> Rotation(x, p, s), (x, p)),
                ("phase Squeeze", s -> Squeeze(x, p, s), (x, p)),
                ("spin Rotation", s -> Rotation(Sx, 3, s), (Sx, Sy, Sz)),
            )
            @testset "$name" begin
                U12, Usum = f(θ1) * f(θ2), f(θ1 + θ2)
                for gen in gens
                    @test agrees(conjugate(gen, U12), conjugate(gen, Usum))
                end
            end
        end
    end

    @testset "a frame change preserves the spectrum" begin
        # The one statement about `conjugate` that needs no reference to how the rules were
        # built: `U'HU` is similar to `H`. On an N-level space the representation is exact, so
        # the eigenvalues have to agree to rounding while the matrix itself moves.
        h4 = NLevelSpace(:sp4, 4)
        τ4(i, j) = Transition(h4, :τ, i, j)
        b4 = NLevelBasis(4)
        m4(q) = Matrix(to_numeric(q, b4).data)
        H4 = 0.4 * τ4(1, 1) + 1.1 * τ4(2, 2) - 0.7 * τ4(3, 3) +
            0.3 * (τ4(1, 2) + τ4(2, 1)) + 0.5 * (τ4(2, 3) + τ4(3, 2)) +
            0.2 * (τ4(1, 4) + τ4(4, 1))
        # Real, symmetric and orthogonal, so `U` is its own inverse and the round trip is a
        # second application rather than a call to `inv`.
        WH = [1.0 1 1 1; 1 -1 1 -1; 1 1 -1 -1; 1 -1 -1 1] / 2
        UH = Rotation(τ4(1, 2), WH)
        @test is_canonical(UH)
        C4 = conjugate(H4, UH)
        @test norm(m4(C4) - m4(H4)) > 1.0            # the frame really moved
        @test sort(real.(eigvals(Hermitian(m4(C4))))) ≈
            sort(real.(eigvals(Hermitian(m4(H4)))))
        @test norm(m4(conjugate(C4, UH)) - m4(H4)) < 1.0e-12
    end

    @testset "Bogoliubov diagonalization" begin
        # The acceptance test the squeezes exist for. Nothing is asserted about the rules:
        # the anomalous term is required to cancel at the textbook angle, and what is left
        # must be the textbook spectrum, zero-point constant included.
        ωv, gv = 1.0, 0.2
        D = conjugate(ωv * a' * a + gv * (a * a + a' * a'), Squeeze(a, atanh(-2gv / ωv) / 2))
        Ω = sqrt(ωv^2 - 4gv^2)
        @test _to_complex(D[Op[a', a]]) ≈ Ω
        @test _to_complex(D[Op[]]) ≈ (Ω - ωv) / 2
        @test abs(_to_complex(D[Op[a, a]])) < 1.0e-12
        @test abs(_to_complex(D[Op[a', a']])) < 1.0e-12
        # The two-mode pair-creation Hamiltonian, same story at half the angle.
        wv, gw = 1.3, 0.35
        D2 = conjugate(
            wv * (ma' * ma + mb' * mb) + gw * (ma' * mb' + ma * mb),
            Squeeze(ma, mb, atanh(-gw / wv) / 2),
        )
        Ω2 = sqrt(wv^2 - gw^2)
        @test _to_complex(D2[Op[ma', ma]]) ≈ Ω2
        @test _to_complex(D2[Op[mb', mb]]) ≈ Ω2
        @test _to_complex(D2[Op[]]) ≈ Ω2 - wv
        @test abs(_to_complex(D2[Op[ma', mb']])) < 1.0e-12
        @test abs(_to_complex(D2[Op[ma, mb]])) < 1.0e-12
    end

    @testset "random compositions" begin
        # The composed transforms above are pairs picked to exercise a named property. This
        # builds two- and three-factor products at random over three modes, so overlapping
        # sites, disjoint sites and repeated factors all turn up without being enumerated.
        hM = FockSpace(:m1) ⊗ FockSpace(:m2) ⊗ FockSpace(:m3)
        q1, q2, q3 = Destroy(hM, :m1, 1), Destroy(hM, :m2, 2), Destroy(hM, :m3, 3)
        @variables s1 s2 s3
        factories = [
            () -> Rotation(q1, s1), () -> Squeeze(q2, s2), () -> Displace(q3, s3),
            () -> Rotation(q1, q2, s1), () -> Squeeze(q2, q3, s2), () -> Rotation(q3, s3),
        ]
        # A composed round trip closes only through the factors' own `cosh`/`sinh` relations,
        # which a bare `simplify` on the difference never sees, so it is read at a point.
        pars = Dict(s1 => 0.4, s2 => 0.3, s3 => 0.55)
        agrees(qa, qb) = all(
            e -> abs(_to_complex(_substitute_cnum(e.second, pars))) < 1.0e-10,
            (qa - qb).arguments,
        )
        rng = MersenneTwister(17)
        for _ in 1:60
            U = foldl(*, [rand(rng, factories)() for _ in 1:rand(rng, 2:3)])
            @test is_canonical(U)
            for gen in generators(U)
                @test agrees(conjugate(conjugate(gen, U), inv(U)), 1 * gen)
                @test agrees(conjugate(adjoint(gen), U), adjoint(conjugate(gen, U)))
            end
        end
    end

    @testset "inference" begin
        U = Rotation(a, θ)
        @test @inferred(conjugate(a' * a, U)) isa QAdd
        @test @inferred(transform(a' * a, U)) isa QAdd
        @test @inferred(conjugate(a, U)) isa QAdd
    end
end
