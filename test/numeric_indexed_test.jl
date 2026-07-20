using SecondQuantizedAlgebra
import SecondQuantizedAlgebra as S
using QuantumOpticsBase
using Symbolics: Symbolics, @variables, Num
using SymbolicUtils: SymbolicUtils
using Test

# Element-wise sup-norm distance between two complex matrices: avoids a hard
# `LinearAlgebra` dependency in the test environment.
_matnorm(A, B) = maximum(abs, A - B)

# Build a Tavis-Cummings Hamiltonian both symbolically (indexed sum) and
# explicitly (per-atom subspaces). Run both through `to_numeric` and verify
# the matrices agree element-by-element on small concrete N.
function _tc_indexed(N, cutoff)
    hc = FockSpace(:cavity)
    ha = NLevelSpace(:atom, 2)
    h = hc ⊗ ha
    a = Destroy(h, :a, 1)
    σ(α, β, k) = IndexedOperator(Transition(h, :σ, α, β, 2), k)
    @variables Δ
    gv(k) = IndexedVariable(:g, k)
    i = Index(h, :i, N, ha)
    H = -Δ * a' * a + Σ(gv(i) * (a' * σ(1, 2, i) + a * σ(2, 1, i)), i)

    bc = FockBasis(cutoff)
    bn = NLevelBasis(2)
    b = let acc = bc
        for _ in 1:N
            acc = acc ⊗ bn
        end
        acc
    end
    sites = Dict{Int, Vector{Int}}(1 => [1], 2 => collect(2:(N + 1)))
    return H, a, b, sites, Δ
end

function _tc_explicit(N, cutoff)
    h = let acc = FockSpace(:c)
        for k in 1:N
            acc = acc ⊗ NLevelSpace(Symbol("atom_", k), 2)
        end
        acc
    end
    a = Destroy(h, :a, 1)
    σ(α, β, k) = Transition(h, Symbol("σ_", k), α, β, k + 1)
    @variables Δ
    gsyms = Num[]
    for k in 1:N
        sym = SymbolicUtils.Sym{SymbolicUtils.SymReal}(Symbol(:gv_, k); type = Real)
        push!(gsyms, Num(sym))
    end
    H = -Δ * a' * a
    for k in 1:N
        H = H + gsyms[k] * (a' * σ(1, 2, k) + a * σ(2, 1, k))
    end
    return H, a, gsyms, Δ
end

# Build the (scalar_subs, subs_ex) pair from a shared (gnums, Δval) so both
# numeric pipelines see the same parameters.
function _build_subs(N, gnums, Δval, Δ_indexed, Δ_ex, gsyms_ex)
    ks_sym = SymbolicUtils.Sym{SymbolicUtils.SymReal}(
        :g;
        type = SymbolicUtils.FnType{Tuple{Int}, Real, Nothing},
        shape = UnitRange{Int}[],
    )
    scalar_subs = Dict{Num, ComplexF64}()
    scalar_subs[Δ_indexed] = Δval
    for k in 1:N
        scalar_subs[Num(ks_sym(k))] = gnums[k]
    end
    subs_ex = Dict{Num, ComplexF64}()
    subs_ex[Δ_ex] = Δval
    for k in 1:N
        subs_ex[gsyms_ex[k]] = gnums[k]
    end
    return scalar_subs, subs_ex
end

@testset "Indexed numeric: Tavis-Cummings" begin
    cutoff = 3
    for N in (2, 3, 4, 5)
        H_idx, a_idx, b, sites, Δ_idx = _tc_indexed(N, cutoff)
        H_ex, a_ex, gsyms_ex, Δ_ex = _tc_explicit(N, cutoff)

        gnums = ComplexF64[ComplexF64(0.1 + 0.05 * k, 0.02 + 0.03 * k) for k in 1:N]
        Δval = ComplexF64(0.5, 0.0)
        scalar_subs, subs_ex =
            _build_subs(N, gnums, Δval, Δ_idx, Δ_ex, gsyms_ex)

        @testset "N=$N: Hamiltonian" begin
            M_idx = to_numeric(
                H_idx, b, sites, Dict{S.QSym, Any}(), scalar_subs,
            )
            H_ex_subbed = substitute(H_ex, subs_ex)
            M_ex = to_numeric(H_ex_subbed, b)
            d1 = sparse(M_idx).data
            d2 = sparse(M_ex).data
            @test _matnorm(d1, d2) < 1.0e-10
        end

        @testset "N=$N: [H, a]" begin
            C_idx = commutator(H_idx, a_idx)
            C_ex = commutator(H_ex, a_ex)
            M_idx = to_numeric(
                C_idx, b, sites, Dict{S.QSym, Any}(), scalar_subs,
            )
            C_ex_subbed = substitute(C_ex, subs_ex)
            M_ex = to_numeric(C_ex_subbed, b)
            d1 = sparse(M_idx).data
            d2 = sparse(M_ex).data
            @test _matnorm(d1, d2) < 1.0e-10
        end
    end
end

@testset "Indexed numeric: d-override and resolved-site naming" begin
    hc = FockSpace(:cavity)
    ha = NLevelSpace(:atom, 2)
    h = hc ⊗ ha
    a = Destroy(h, :a, 1)
    sig(α, β, k) = IndexedOperator(Transition(h, :σ, α, β, 2), k)
    i = Index(h, :i, 2, ha)
    H = Σ(a' * sig(1, 2, i) + a * sig(2, 1, i), i)

    bc = FockBasis(2)
    bn = NLevelBasis(2)
    b = bc ⊗ bn ⊗ bn
    sites = Dict{Int, Vector{Int}}(1 => [1], 2 => [2, 3])

    # A resolved per-site index keeps its real name (not the `:_` sentinel).
    resolved = Index(i.name_id, i.range_id, i.space_index, Int32(2))
    @test index_name(resolved) === :i
    @test has_index(resolved)

    # A `d` override keyed by the abstract indexed op now applies on the sites
    # path; it was silently ignored when resolved ops carried name_id 0.
    M_plain = to_numeric(H, b, sites, Dict{S.QSym, Any}())
    custom = embed(b, 2, transition(bn, 1, 2))
    M_over = to_numeric(H, b, sites, Dict{S.QSym, Any}(sig(1, 2, i) => custom))
    @test sparse(M_plain).data != sparse(M_over).data
end

@testset "Indexed numeric: ne-constrained double sum" begin
    # A double sum `Σ_{i,j} σ⁺_i σ⁻_j` splits into a diagonal part (`i == j`) and an
    # off-diagonal part carrying an `ne = (i, j)` constraint that `_violates_ne` uses to
    # skip the `i == j` combinations during the unroll. The result must equal the full
    # explicit sum over all site pairs.
    ha = NLevelSpace(:atom, 2)
    σ(α, β, k) = IndexedOperator(Transition(ha, :σ, α, β), k)
    i = Index(ha, :i, 2, ha)
    j = Index(ha, :j, 2, ha)
    H = Σ(σ(1, 2, i) * σ(2, 1, j), i, j)

    bn = NLevelBasis(2)
    b = bn ⊗ bn
    sites = Dict{Int, Vector{Int}}(1 => [1, 2])
    M = to_numeric(H, b, sites, Dict{S.QSym, Any}())

    he = NLevelSpace(:a1, 2) ⊗ NLevelSpace(:a2, 2)
    t(α, β, k) = Transition(he, Symbol(:s, k), α, β, k)
    Hex = t(1, 2, 1) * t(2, 1, 1) + t(1, 2, 1) * t(2, 1, 2) +
        t(1, 2, 2) * t(2, 1, 1) + t(1, 2, 2) * t(2, 1, 2)
    Mex = to_numeric(Hex, b)
    @test _matnorm(sparse(M).data, sparse(Mex).data) < 1.0e-10
end

@testset "Indexed numeric: real-valued scalar_subs" begin
    # Real (non-`Complex`-typed) `scalar_subs` values exercise the real-only substitution
    # legs of `_split_scalar_subs`/`_apply_scalar_subs` (the `has_imag == false` path).
    hc = FockSpace(:cavity)
    ha = NLevelSpace(:atom, 2)
    h = hc ⊗ ha
    a = Destroy(h, :a, 1)
    σ(α, β, k) = IndexedOperator(Transition(h, :σ, α, β, 2), k)
    gv(k) = IndexedVariable(:g, k)
    i = Index(h, :i, 2, ha)
    H = Σ(gv(i) * (a' * σ(1, 2, i) + a * σ(2, 1, i)), i)

    bc = FockBasis(2)
    bn = NLevelBasis(2)
    b = bc ⊗ bn ⊗ bn
    sites = Dict{Int, Vector{Int}}(1 => [1], 2 => [2, 3])

    ks_sym = SymbolicUtils.Sym{SymbolicUtils.SymReal}(
        :g;
        type = SymbolicUtils.FnType{Tuple{Int}, Real, Nothing},
        shape = UnitRange{Int}[],
    )
    scalar_subs = Dict{Num, Float64}(Num(ks_sym(1)) => 0.3, Num(ks_sym(2)) => 0.7)
    M = to_numeric(H, b, sites, Dict{S.QSym, Any}(), scalar_subs)

    # Explicit reference on the same layout.
    he = FockSpace(:c) ⊗ NLevelSpace(:a1, 2) ⊗ NLevelSpace(:a2, 2)
    ae = Destroy(he, :a, 1)
    te(α, β, k) = Transition(he, Symbol(:s, k), α, β, k + 1)
    Hex = 0.3 * (ae' * te(1, 2, 1) + ae * te(2, 1, 1)) +
        0.7 * (ae' * te(1, 2, 2) + ae * te(2, 1, 2))
    Mex = to_numeric(Hex, b)
    @test _matnorm(sparse(M).data, sparse(Mex).data) < 1.0e-10
end

@testset "Indexed numeric: slot resolution edge cases" begin
    hc = FockSpace(:cavity)
    ha = NLevelSpace(:atom, 2)
    h = hc ⊗ ha
    a = Destroy(h, :a, 1)
    σ(α, β, k) = IndexedOperator(Transition(h, :σ, α, β, 2), k)
    i = Index(h, :i, 2, ha)
    H = Σ(a' * σ(1, 2, i) + a * σ(2, 1, i), i)

    bc = FockBasis(2)
    bn = NLevelBasis(2)
    b = bc ⊗ bn ⊗ bn
    sites = Dict{Int, Vector{Int}}(1 => [1], 2 => [2, 3])

    # A non-indexed subspace (the cavity) may be omitted from `sites`: its operator falls
    # back to its own `space_index` slot, so the result is unchanged.
    sites_no_cavity = Dict{Int, Vector{Int}}(2 => [2, 3])
    M_full = to_numeric(H, b, sites, Dict{S.QSym, Any}())
    M_omit = to_numeric(H, b, sites_no_cavity, Dict{S.QSym, Any}())
    @test _matnorm(sparse(M_full).data, sparse(M_omit).data) < 1.0e-10

    # An operator carrying an unresolved (non-integer) index on a multi-slot subspace is
    # ambiguous and must error. `σ_j` never binds to the summed index `i`.
    j = Index(h, :j, 2, ha)
    Hbad = Σ(σ(1, 2, i) * σ(1, 2, j), i)
    @test_throws ArgumentError to_numeric(Hbad, b, sites, Dict{S.QSym, Any}())
end

@testset "Indexed numeric: invalid site layouts are rejected" begin
    h = FockSpace(:site)
    a = Destroy(h, :a)
    i = Index(h, :i, 2, h)
    H = Σ(IndexedOperator(a, i), i)

    @test_throws ArgumentError to_numeric(H, FockBasis(2), Dict(1 => [1, 2]))
    b = FockBasis(2) ⊗ FockBasis(2)
    @test_throws ArgumentError to_numeric(H, b, Dict(1 => [1, 1]))
    @test_throws ArgumentError to_numeric(H, b, Dict(1 => [1, 3]))
end
