using PrecompileTools: @setup_workload, @compile_workload

precompile(canonicalize!, (QTermDict, Vector{Op}, CNum, Vector{NonEqualPair}))
precompile(addto!, (QTermDict, Vector{Op}, CNum, Vector{NonEqualPair}))
precompile(iszero, (CNum,))
precompile(-, (CNum,))
precompile(merge_unique, (Vector{Index}, Vector{Index}))
precompile(merge_ne, (Vector{NonEqualPair}, Vector{NonEqualPair}))

for T in (Int, Float64, ComplexF64, Rational{Int}, Num)
    precompile(*, (T, Op))
    precompile(*, (Op, T))
end
precompile(get_prefactor, (QAdd,))
precompile(get_operators, (QAdd,))
precompile(get_variables, (QAdd,))
precompile(to_num, (CNum,))

@setup_workload begin
    @compile_workload begin
        let io = IOBuffer()
            h = FockSpace(:a); a = Destroy(h, :a); ad = Create(h, :a)
            show(io, a * ad); show(io, a + ad); show(io, a'); show(io, a - 2 * ad)
            show(io, commutator(a, ad)); show(io, (a + ad) * (a - ad))
            hn = NLevelSpace(:atom, 3); s12 = Transition(hn, :σ, 1, 2); s21 = Transition(hn, :σ, 2, 1)
            show(io, s12 * s21); show(io, normal_order(s21 * s12))
            show(io, expand_completeness(Transition(hn, :σ, 2, 2)))
            hs = SpinSpace(:s); show(io, commutator(Spin(hs, :S, 1), Spin(hs, :S, 2)))
            hp = PauliSpace(:p); show(io, Pauli(hp, :σ, 1) * Pauli(hp, :σ, 2))
            @variables g Δ
            H = g * ad * a + Δ * (a + ad); show(io, H); show(io, H')
            show(io, average(H)); show(io, undo_average(average(a * ad)))
            show(io, simplify(H)); substitute(H, Dict(g => 1.0))
            let i = Index(h, :i, 3, h)
                show(io, Σ(IndexedOperator(ad, i) * IndexedOperator(a, i), i))
            end
            L = hamiltonian_action(H) + g * dissipator(a)
            show(io, L); show(io, SuperOperator(L)); show(io, transpose(L)(ad * a))
            show(io, L ∘ L); show(io, sprepost(a, ad)(H))
        end
    end
end
