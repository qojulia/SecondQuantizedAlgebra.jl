using PrecompileTools: @setup_workload, @compile_workload

# === Explicit directives ===
# Internal entry points the workload does not fully exercise (indexed `_merge_*`).
precompile(_canonicalize!, (QTermDict, Vector{Op}, CNum, Vector{NonEqualPair}))
precompile(_addto!, (QTermDict, Vector{Op}, CNum, Vector{NonEqualPair}))
precompile(_iszero_cnum, (CNum,))
precompile(_neg_cnum, (CNum,))
precompile(_merge_unique, (Vector{Index}, Vector{Index}))
precompile(_merge_ne, (Vector{NonEqualPair}, Vector{NonEqualPair}))

# Scalar prefactor products across the numeric types users write literally, in both
# argument orders (a type loop; the workload only exercises `Int`/`Num` coefficients).
for T in (Int, Float64, ComplexF64, Rational{Int}, Num)
    precompile(*, (T, Op))
    precompile(*, (Op, T))
end
# Public single-term accessors, not reached by the display-driven workload.
precompile(prefactor, (QAdd,))
precompile(operators, (QAdd,))
precompile(to_num, (CNum,))

# === Representative workload ===
# Covers the common user surface: each operator kind, the full arithmetic/
# canonicalization pipeline (including QAdd*QAdd), normal_order, symbolic (parameter)
# coefficients, average/undo_average, simplify, substitute, expand_completeness, an
# indexed sum, and display of every result (the universal REPL path). Wrapped in
# `let` (Makie-style) so no workload local is serialized into the pkgimage.
#
# Balance (measured, Julia 1.12 / SymbolicUtils 4.40): precompilation ~6.1s -> ~11s
# (paid once per install), cutting first-call latency of a full build/manipulate/
# average/simplify/show workflow from ~3.5s to ~0.2s (paid every session). `average`
# and `simplify` share the SymbolicUtils machinery (together ~+3s), and the indexed
# and substitute lines ~+1s; drop those lines to trim precompilation if leaner is
# wanted.
@setup_workload begin
    @compile_workload begin
        let io = IOBuffer()
            h = FockSpace(:a); a = Destroy(h, :a); ad = Create(h, :a)
            show(io, a*ad); show(io, a + ad); show(io, a'); show(io, a - 2*ad)
            show(io, commutator(a, ad)); show(io, (a + ad)*(a - ad))
            hn = NLevelSpace(:atom, 3); s12 = Transition(hn, :σ, 1, 2); s21 = Transition(hn, :σ, 2, 1)
            show(io, s12*s21); show(io, normal_order(s21*s12))
            show(io, expand_completeness(Transition(hn, :σ, 2, 2)))
            hs = SpinSpace(:s); show(io, commutator(Spin(hs, :S, 1), Spin(hs, :S, 2)))
            hp = PauliSpace(:p); show(io, Pauli(hp, :σ, 1)*Pauli(hp, :σ, 2))
            @variables g Δ
            H = g*ad*a + Δ*(a + ad); show(io, H); show(io, H')
            show(io, average(H)); show(io, undo_average(average(a*ad)))
            show(io, simplify(H)); substitute(H, Dict(g => 1.0))
            let i = Index(h, :i, 3, h)
                show(io, Σ(IndexedOperator(ad, i)*IndexedOperator(a, i), i))
            end
        end
    end
end
