using SecondQuantizedAlgebra
using Symbolics: Symbolics, @variables
using Test

import SecondQuantizedAlgebra: SuperOperator, sprepost, spre, spost, hamiltonian, channels,
    canonical_isless, independent_basis, support_basis, GKSLFormError, QAdd, to_num, jump,
    partner

side(x) = x * one(x)
numeric(c) = complex(Symbolics.value(real(to_num(c))), Symbolics.value(imag(to_num(c))))
same_map(A, B) = expand_completeness(SuperOperator(A)) == expand_completeness(SuperOperator(B))

function coordinates(S, basis)
    L = Lindbladian(S, basis)
    return hamiltonian(L), kossakowski(L, basis)
end

@testset "GKSL coordinates" begin
    h = FockSpace(:c)
    @qnumbers a::Destroy(h)
    ha = NLevelSpace(:atom, (:g, :e))
    σ(i, j) = Transition(ha, :σ, i, j)
    @variables γ::Real κ::Real Δ::Real g::Real

    @testset "bases are canonical modulo the identity" begin
        basis = JumpBasis(a + 2, a' * a)
        @test collect(basis) == side.([a, a' * a])
        @test basis[2] == side(a' * a)
        @test length(basis) == 2
        @test eltype(basis) == QAdd
        @test basis == JumpBasis([a, a' * a + 1])
        @test hash(basis) == hash(JumpBasis([a, a' * a + 1]))
        @test basis != JumpBasis(a' * a, a)
        @test isempty(JumpBasis())
        @test collect(JumpBasis(σ(:g, :g))) == [-σ(:e, :e)]
        @test_throws ArgumentError JumpBasis(a, 2a + 1)
        @test_throws ArgumentError JumpBasis(one(a))
    end

    @testset "independent bases keep the first of each direction" begin
        @test collect(independent_basis([a, 2a + 1, one(a), a' * a, a' * a + a, a'])) ==
            side.([a, a' * a, a'])
        @test collect(independent_basis([a + a', a - a', a])) == side.([a + a', a - a'])
        @test collect(independent_basis([g * a + a' * a, a' * a, a])) == [g * a + a' * a, side(a' * a)]
        @test isempty(independent_basis(QAdd[]))
    end

    @testset "bases display their canonical operators" begin
        basis = JumpBasis(a + 3, a^2)
        @test sprint(show, basis) == "JumpBasis(a, a * a)"
        @test sprint(show, MIME"text/plain"(), basis) ==
            "JumpBasis with 2 operators:\n  1: a\n  2: a * a"
    end

    @testset "Kossakowski matrix in diagonal and non-orthogonal bases" begin
        L = hamiltonian_action(Δ * a' * a) + κ * dissipator(a) + γ * dissipator(a' * a)
        @test isequal(kossakowski(L, JumpBasis(a, a' * a)), [κ 0; 0 γ])
        @test isequal(kossakowski(SuperOperator(L), JumpBasis(a' * a, a)), [γ 0; 0 κ])
        d = kossakowski(L, JumpBasis(a, a + a' * a))
        @test isequal(simplify.(d), simplify.([κ + γ -γ; -γ γ]))
        @test isequal(kossakowski(SuperOperator(L), JumpBasis(a, a + a' * a)), d)
        @test isequal(kossakowski(dissipator(a, a' * a) + dissipator(a' * a, a), JumpBasis(a, a' * a)), [0 1; 1 0])
        @test isequal(kossakowski(hamiltonian_action(Δ * a' * a), JumpBasis()), zeros(Int, 0, 0))
    end

    @testset "reading a Lindbladian agrees with reading its map" begin
        basis = JumpBasis(2a, a + a' * a)
        models = (
            κ * dissipator(a) + γ * dissipator(a' * a),
            κ * dissipator(3a + 1, a' * a) + κ * dissipator(a' * a, 3a + 1),
            0.5 * dissipator(a + 0.25 * a' * a),
        )
        for L in models
            @test isequal(kossakowski(L, basis), kossakowski(SuperOperator(L), basis))
        end
        @test_throws ArgumentError kossakowski(dissipator(a^2), JumpBasis(a))
        @test_throws ArgumentError kossakowski(dissipator(a), JumpBasis(a); rtol = -1)
        @test_throws GKSLFormError kossakowski(hamiltonian_action(im * a' * a) + dissipator(a), JumpBasis(a))
        floating = hamiltonian_action(0.4 * a' * a) + 0.7 * dissipator(0.3 * a + 0.7 * a' * a) +
            0.11 * dissipator(a' + 0.9 * a^2)
        @test_throws GKSLFormError kossakowski(floating, support_basis(floating))
        @test isequal(
            kossakowski(floating, support_basis(floating); rtol = 1.0e-12),
            kossakowski(SuperOperator(floating), support_basis(floating); rtol = 1.0e-12),
        )
    end

    @testset "a shifted jump is a Hamiltonian change" begin
        L = dissipator(a + 2)
        @test isequal(kossakowski(L, JumpBasis(a)), fill(1, 1, 1))
        recovered = Lindbladian(SuperOperator(L))
        @test recovered == L
        @test hamiltonian(recovered) == im * (a - a')
    end

    @testset "Lindbladian(SuperOperator) is a left inverse" begin
        models = (
            hamiltonian_action(Δ * a' * a) + κ * dissipator(a) + γ * dissipator(a' * a),
            hamiltonian_action(g * (a + a')) + κ * dissipator(a^2) + γ * dissipator(a, a' * a) +
                γ * dissipator(a' * a, a),
            hamiltonian_action(Δ * σ(:e, :e)) + γ * dissipator(σ(:g, :e)) +
                κ * dissipator(σ(:g, :g)),
            (1 // 10) * dissipator(a) + (1 // 4) * dissipator(a' * a),
            (1 // 10) * dissipator(a + (1 // 5) * a' * a) + (3 // 10) * dissipator(a' + (7 // 10) * a^2),
            hamiltonian_action(Δ * a' * a),
        )
        for L in models
            S = SuperOperator(L)
            @test same_map(Lindbladian(S), L)
            @test same_map(Lindbladian(S, support_basis(S)), L)
        end
    end

    @testset "support basis applies finite-level completeness" begin
        S = SuperOperator(γ * dissipator(σ(:g, :e)) + κ * dissipator(σ(:g, :g)))
        @test collect(support_basis(S)) == side.([σ(:g, :e), σ(:e, :e)])
        @test isequal(kossakowski(S, support_basis(S)), [γ 0; 0 κ])
        @test isempty(support_basis(hamiltonian_action(Δ * a' * a)))
    end

    @testset "maps outside GKSL form are rejected" begin
        basis = JumpBasis(a)
        @test_throws GKSLFormError Lindbladian(spre(a))
        @test_throws GKSLFormError Lindbladian(SuperOperator(hamiltonian_action(im * a' * a)))
        @test_throws GKSLFormError Lindbladian(spre(one(a)))
        @test_throws "does not preserve the trace" Lindbladian(sprepost(a, a'))
        @test_throws "not Hermitian modulo the identity" Lindbladian(sprepost(a, a') - spre(a' * a))
        hb = Index(h, :i, 10, h)
        summed = Σ(dissipator(IndexedOperator(Destroy(h, :b), hb)), hb)
        @test_throws ArgumentError kossakowski(summed, JumpBasis(a))
        @test_throws GKSLFormError kossakowski(sprepost(a, a), JumpBasis(a, a'))
        @test_throws ArgumentError kossakowski(SuperOperator(dissipator(a' * a)), basis)
        @test_throws ArgumentError kossakowski(SuperOperator(dissipator(a)), basis; rtol = -1)
        @test_throws ArgumentError JumpBasis(Σ(IndexedOperator(Destroy(h, :b), Index(h, :i, 10, h)), Index(h, :i, 10, h)))
    end

    @testset "rescaling the basis does not hide a term" begin
        defect = SuperOperator(dissipator(a)) + 0.3 * sprepost(a^2, (a')^2)
        outside = SuperOperator(dissipator(a) + 0.3 * dissipator(a^2))
        for s in (1.0e-6, 1.0, 1.0e6), rtol in (0, 1.0e-6)
            @test_throws "does not preserve the trace" Lindbladian(defect, JumpBasis(s * a); rtol)
            @test_throws ArgumentError Lindbladian(outside, JumpBasis(s * a); rtol)
        end
    end

    @testset "float maps are exact unless rtol opts in" begin
        J1 = 0.3 * a + 0.7 * a' * a
        J2 = a' + 0.9 * a^2
        S = SuperOperator(hamiltonian_action(0.4 * a' * a) + 0.7 * dissipator(J1) + 0.11 * dissipator(J2))
        err = try
            Lindbladian(S)
        catch e
            e
        end
        @test err isa GKSLFormError
        @test occursin("1//10", sprint(showerror, err))
        recovered = Lindbladian(S; rtol = 1.0e-12)
        difference = expand_completeness(SuperOperator(recovered)) - expand_completeness(S)
        @test all(((_, c),) -> abs(numeric(c)) < 1.0e-14, difference)
        exact_err = try
            Lindbladian(spre(a))
        catch e
            e
        end
        @test !occursin("1//10", sprint(showerror, exact_err))
    end

    @testset "a positive rtol is relative to the same parameter monomial" begin
        S = SuperOperator(2.0 * dissipator(a) + 1.0e-11 * hamiltonian_action(a' * a))
        @test hamiltonian(Lindbladian(S, JumpBasis(a))) == 1.0e-11 * a' * a
        @test iszero(hamiltonian(Lindbladian(S, JumpBasis(a); rtol = 1.0e-10)))
        T = SuperOperator(1.0 * dissipator(a) + 1.0e-11 * γ * dissipator(a' * a))
        @test isequal(kossakowski(T, JumpBasis(a, a' * a); rtol = 1.0e-10), [1.0 0; 0 1.0e-11 * γ])
        U = SuperOperator(g * dissipator(a, a' * a) + (1.0 + 1.0e-6) * g * dissipator(a' * a, a))
        @test_throws GKSLFormError kossakowski(U, JumpBasis(a, a' * a))
        @test_throws GKSLFormError kossakowski(U, JumpBasis(a, a' * a); rtol = 1.0e-10)
        @test kossakowski(U, JumpBasis(a, a' * a); rtol = 1.0e-5) isa Matrix
    end

    @testset "rescaled floating bases" begin
        L = 0.1 * dissipator(a) + 0.3 * dissipator(a' * a)
        for s in (1.0e-6, 1.0, 1.0e6)
            d = kossakowski(SuperOperator(L), JumpBasis(s * a, s * a' * a); rtol = 1.0e-12)
            @test map(c -> numeric(c) * s^2, d) ≈ [0.1 0; 0 0.3]
        end
    end

    @testset "channels are returned in channel order" begin
        L = γ * dissipator(a' * a) + κ * dissipator(a) + g * dissipator(a^2)
        entries = channels(L)
        @test entries isa Vector
        @test issorted(first.(entries); lt = canonical_isless)
        @test length(entries) == 3
        @test sum(rate * dissipator(jump(ch), partner(ch)) for (ch, rate) in entries) == L
    end

    @testset "complex coherences and identity shifts of H are exact" begin
        H = Δ * a' * a
        S = spre(-im * (H + 7)) + spost(im * (H + 7)) + SuperOperator(dissipator(a + im * a^2))
        H_extracted, d = coordinates(S, JumpBasis(a, a^2))
        @test isequal(d, [1 -im; im 1])
        @test iszero(simplify(H_extracted - H))
        @test iszero(simplify(hamiltonian(Lindbladian(S)) - H))
    end

    @testset "non-polynomial coefficients" begin
        L = hamiltonian_action(sin(g) * a' * a) + exp(κ) * dissipator(a) + cos(g)^2 * dissipator(a' * a)
        S = SuperOperator(L)
        @test same_map(Lindbladian(S), L)
        basis = JumpBasis(a, sin(g) * a + a' * a)
        @test same_map(Lindbladian(S, basis), L)
        @test isequal(kossakowski(L, basis), kossakowski(S, basis))
        @test length(independent_basis([a, sin(g) * a, a' * a + cos(g) * a])) == 2
    end

    @testset "inactive basis operators stay dark" begin
        @test isequal(kossakowski(dissipator(a), JumpBasis(a, a^2)), [1 0; 0 0])
    end

    @testset "adjoint thermal channels keep independent rates" begin
        H, d = coordinates(SuperOperator(κ * dissipator(a) + γ * dissipator(a')), JumpBasis(a, a'))
        @test isequal(d, [κ 0; 0 γ])
        @test iszero(simplify(H))
    end

    @testset "Pauli algebra closes inside the coordinates" begin
        hs = PauliSpace(:spin)
        σx, σy, σz = Pauli(hs, :σ, 1), Pauli(hs, :σ, 2), Pauli(hs, :σ, 3)
        L = hamiltonian_action(Δ * σz) + dissipator(σx + im * σy)
        H, d = coordinates(SuperOperator(L), JumpBasis(σx, σy))
        @test isequal(d, [1 -im; im 1])
        @test iszero(simplify(H - Δ * σz))
    end

    @testset "canonical phase-space operators" begin
        hp = PhaseSpace(:phase)
        x, p = Position(hp, :x), Momentum(hp, :p)
        H, d = coordinates(SuperOperator(dissipator(x + im * p)), JumpBasis(x, p))
        @test isequal(d, [1 -im; im 1])
        @test iszero(simplify(H))
        squeezing = hamiltonian_action(g * (x * p + p * x) / 2)
        @test same_map(Lindbladian(SuperOperator(squeezing)), squeezing)
        floating = hamiltonian_action(0.5 * (x * p + p * x))
        @test same_map(Lindbladian(SuperOperator(floating)), floating)
    end

    @testset "spin operators" begin
        hs = SpinSpace(:spin)
        Sx, Sy, Sz = Spin(hs, :S, 1), Spin(hs, :S, 2), Spin(hs, :S, 3)
        L = hamiltonian_action(Δ * Sz) + γ * dissipator(Sx + im * Sy)
        H, d = coordinates(SuperOperator(L), JumpBasis(Sx, Sy))
        @test isequal(d, [γ -im * γ; im * γ γ])
        @test iszero(simplify(H - Δ * Sz))
    end

    @testset "mixed product-space operators" begin
        hm = FockSpace(:boson) ⊗ PauliSpace(:spin)
        b = Destroy(hm, :b)
        σx = Pauli(hm, :σ, 1)
        H, d = coordinates(SuperOperator(dissipator(b + 2σx)), JumpBasis(b, σx))
        @test isequal(d, [1 2; 2 4])
        @test iszero(simplify(H))
    end

    @testset "bosonic monomial basis needs no cutoff" begin
        d = kossakowski(dissipator(a + 2a^2 - im * a' * a^2), JumpBasis(a, a^2, a' * a^2))
        @test isequal(d, [1 2 im; 2 4 2im; -im -2im 1])
    end
end
