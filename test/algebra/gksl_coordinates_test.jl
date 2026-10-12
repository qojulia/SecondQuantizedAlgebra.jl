using SecondQuantizedAlgebra
using Symbolics: @variables
using Test

import SecondQuantizedAlgebra: SuperOperator, sprepost, spre, spost, hamiltonian, channels,
    channel_order_key, directions, independent_frame, support_frame, gksl_coordinates, GKSLFormError

side(x) = x * one(x)
same_map(A, B) = expand_completeness(SuperOperator(A)) == expand_completeness(SuperOperator(B))

@testset "GKSL coordinates" begin
    h = FockSpace(:c)
    @qnumbers a::Destroy(h)
    ha = NLevelSpace(:atom, (:g, :e))
    σ(i, j) = Transition(ha, :σ, i, j)
    @variables γ::Real κ::Real Δ::Real g::Real

    @testset "frames are canonical modulo the identity" begin
        frame = DissipativeFrame(a + 2, a' * a)
        @test directions(frame) == side.([a, a' * a])
        @test frame == DissipativeFrame([a, a' * a + 1])
        @test hash(frame) == hash(DissipativeFrame([a, a' * a + 1]))
        @test frame != DissipativeFrame(a' * a, a)
        @test isempty(DissipativeFrame())
        @test directions(DissipativeFrame(σ(:g, :g))) == [-σ(:e, :e)]
        @test_throws ArgumentError DissipativeFrame(a, 2a + 1)
        @test_throws ArgumentError DissipativeFrame(one(a))
        @test directions(independent_frame([a, 2a + 1, one(a), a' * a, a' * a + a, a'])) ==
            side.([a, a' * a, a'])
    end

    @testset "frames display their canonical directions" begin
        frame = DissipativeFrame(a + 3, a^2)
        @test sprint(show, frame) == "DissipativeFrame(a, a * a)"
        @test sprint(show, MIME"text/plain"(), frame) ==
            "DissipativeFrame with 2 directions:\n  1: a\n  2: a * a"
    end

    @testset "Kossakowski matrix in diagonal and non-orthogonal frames" begin
        L = hamiltonian_action(Δ * a' * a) + κ * dissipator(a) + γ * dissipator(a' * a)
        @test isequal(kossakowski(L, DissipativeFrame(a, a' * a)), [κ 0; 0 γ])
        @test isequal(kossakowski(SuperOperator(L), DissipativeFrame(a' * a, a)), [γ 0; 0 κ])
        d = kossakowski(L, DissipativeFrame(a, a + a' * a))
        @test isequal(simplify.(d), simplify.([κ + γ -γ; -γ γ]))
        @test isequal(kossakowski(dissipator(a, a' * a) + dissipator(a' * a, a), DissipativeFrame(a, a' * a)), [0 1; 1 0])
        @test isequal(kossakowski(hamiltonian_action(Δ * a' * a), DissipativeFrame()), zeros(Int, 0, 0))
    end

    @testset "a shifted jump is a Hamiltonian change" begin
        L = dissipator(a + 2)
        @test isequal(kossakowski(L, DissipativeFrame(a)), fill(1, 1, 1))
        recovered = Lindbladian(SuperOperator(L))
        @test recovered == L
        @test hamiltonian(recovered) == im * (a - a')
        H, d = gksl_coordinates(L, DissipativeFrame(a))
        @test simplify(H - im * (a - a')) == simplify(H - hamiltonian(recovered))
        @test isequal(d, kossakowski(L, DissipativeFrame(a)))
    end

    @testset "Lindbladian(SuperOperator) is a left inverse" begin
        models = (
            hamiltonian_action(Δ * a' * a) + κ * dissipator(a) + γ * dissipator(a' * a),
            hamiltonian_action(g * (a + a')) + κ * dissipator(a^2) + γ * dissipator(a, a' * a) +
                γ * dissipator(a' * a, a),
            hamiltonian_action(Δ * σ(:e, :e)) + γ * dissipator(σ(:g, :e)) +
                κ * dissipator(σ(:g, :g)),
            (1 // 10) * dissipator(a) + 0.25 * dissipator(a' * a),
            hamiltonian_action(Δ * a' * a),
        )
        for L in models
            S = SuperOperator(L)
            recovered = Lindbladian(S)
            @test same_map(recovered, L)
            @test same_map(Lindbladian(S, support_frame(S)), L)
        end
    end

    @testset "support frame applies finite-level completeness" begin
        S = SuperOperator(γ * dissipator(σ(:g, :e)) + κ * dissipator(σ(:g, :g)))
        @test directions(support_frame(S)) == side.([σ(:g, :e), σ(:e, :e)])
        @test isequal(kossakowski(S, support_frame(S)), [γ 0; 0 κ])
    end

    @testset "maps outside GKSL form are rejected" begin
        frame = DissipativeFrame(a)
        @test_throws GKSLFormError Lindbladian(spre(a))
        @test_throws GKSLFormError Lindbladian(SuperOperator(hamiltonian_action(im * a' * a)))
        @test_throws GKSLFormError kossakowski(sprepost(a, a), DissipativeFrame(a, a'))
        @test_throws ArgumentError kossakowski(dissipator(a' * a), frame)
        @test_throws ArgumentError DissipativeFrame(Σ(IndexedOperator(Destroy(h, :b), Index(h, :i, 10, h)), Index(h, :i, 10, h)))
    end

    @testset "channels are returned in channel order" begin
        L = γ * dissipator(a' * a) + κ * dissipator(a) + g * dissipator(a^2)
        entries = channels(L)
        @test entries isa Vector
        @test issorted(entries; by = entry -> channel_order_key(entry.first))
        @test Set(entries) == Set(collect(L.channels))
    end

    @testset "complex coherences and identity shifts of H are exact" begin
        H = Δ * a' * a
        S = spre(-im * (H + 7)) + spost(im * (H + 7)) + SuperOperator(dissipator(a + im * a^2))
        H_extracted, d = gksl_coordinates(S, DissipativeFrame(a, a^2))
        @test isequal(d, [1 -im; im 1])
        @test iszero(simplify(H_extracted - H))
        @test iszero(simplify(hamiltonian(Lindbladian(S)) - H))
    end

    @testset "inactive frame directions stay dark" begin
        @test isequal(kossakowski(dissipator(a), DissipativeFrame(a, a^2)), [1 0; 0 0])
    end

    @testset "adjoint thermal channels keep independent rates" begin
        H, d = gksl_coordinates(κ * dissipator(a) + γ * dissipator(a'), DissipativeFrame(a, a'))
        @test isequal(d, [κ 0; 0 γ])
        @test iszero(simplify(H))
    end

    @testset "Pauli algebra closes inside the coordinates" begin
        hs = PauliSpace(:spin)
        σx, σy, σz = Pauli(hs, :σ, 1), Pauli(hs, :σ, 2), Pauli(hs, :σ, 3)
        L = hamiltonian_action(Δ * σz) + dissipator(σx + im * σy)
        H, d = gksl_coordinates(L, DissipativeFrame(σx, σy))
        @test isequal(d, [1 -im; im 1])
        @test iszero(simplify(H - Δ * σz))
    end

    @testset "canonical phase-space directions" begin
        hp = PhaseSpace(:phase)
        x, p = Position(hp, :x), Momentum(hp, :p)
        H, d = gksl_coordinates(dissipator(x + im * p), DissipativeFrame(x, p))
        @test isequal(d, [1 -im; im 1])
        @test iszero(simplify(H))
    end

    @testset "mixed product-space directions" begin
        hm = FockSpace(:boson) ⊗ PauliSpace(:spin)
        b = Destroy(hm, :b)
        σx = Pauli(hm, :σ, 1)
        H, d = gksl_coordinates(dissipator(b + 2σx), DissipativeFrame(b, σx))
        @test isequal(d, [1 2; 2 4])
        @test iszero(simplify(H))
    end

    @testset "bosonic monomial frame needs no cutoff" begin
        d = kossakowski(dissipator(a + 2a^2 - im * a' * a^2), DissipativeFrame(a, a^2, a' * a^2))
        @test isequal(d, [1 2 im; 2 4 2im; -im -2im 1])
    end
end
