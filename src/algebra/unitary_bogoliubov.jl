# === General exact bosonic Bogoliubov transformations ===

function bogoliubov_modes(modes::AbstractVector{Op})
    isempty(modes) && unitary_error("`Bogoliubov` needs at least one Fock mode")
    lowering_modes = Op[]
    sizehint!(lowering_modes, length(modes))
    seen = Set{SiteKey}()
    for mode in modes
        lowering_mode = fock_or_throw(mode, "`Bogoliubov`")
        key = site_key(lowering_mode)
        key in seen && unitary_error(
            "`Bogoliubov` received the same Fock mode more than once: `$lowering_mode`",
        )
        push!(seen, key)
        push!(lowering_modes, lowering_mode)
    end
    return lowering_modes
end

bogoliubov_modes(mode::Op) = bogoliubov_modes(Op[mode])
bogoliubov_modes(modes::Tuple{Vararg{Op}}) = bogoliubov_modes(Op[modes...])

function bogoliubov_basis(modes::Vector{Op})
    n = length(modes)
    basis = Vector{Op}(undef, 2 * n)
    for i in 1:n
        basis[i] = modes[i]
        basis[n + i] = adjoint(modes[i])
    end
    return basis
end

function bogoliubov_matrix(S::AbstractMatrix, n::Int)
    dimension = 2 * n
    size(S) == (dimension, dimension) || unitary_error(
        "`Bogoliubov` needs a $dimension×$dimension Nambu matrix for $n modes; got $(size(S))",
    )
    matrix = Matrix{CNum}(undef, dimension, dimension)
    for j in 1:dimension, i in 1:dimension
        matrix[i, j] = to_cnum(S[i, j])
    end
    return matrix
end

function bogoliubov_matrix(U::AbstractMatrix, V::AbstractMatrix, n::Int)
    size(U) == (n, n) || unitary_error(
        "`Bogoliubov` needs an $n×$n `U` block; got $(size(U))",
    )
    size(V) == (n, n) || unitary_error(
        "`Bogoliubov` needs an $n×$n `V` block; got $(size(V))",
    )
    dimension = 2 * n
    matrix = Matrix{CNum}(undef, dimension, dimension)
    for j in 1:n, i in 1:n
        u = to_cnum(U[i, j])
        v = to_cnum(V[i, j])
        matrix[i, j] = u
        matrix[i, n + j] = v
        matrix[n + i, j] = conj_cnum(v)
        matrix[n + i, n + j] = conj_cnum(u)
    end
    return matrix
end

function bogoliubov_action(modes::Vector{Op}, matrix::Matrix{CNum})
    basis = bogoliubov_basis(modes)
    return AffineAction(
        AFFINE_BOSONIC_NAMBU, basis, matrix, fill(CNUM_ZERO, length(basis)),
    )
end

exact_bogoliubov(modes::Vector{Op}, matrix::Matrix{CNum}) =
    canonical_transform(bogoliubov_action(modes, matrix))

"""
    Bogoliubov(modes, S)

Construct a bosonic Bogoliubov transformation in Nambu ordering
`(a₁, …, aₙ, a₁', …, aₙ')`.

`S` is required by contract to preserve adjoints and the bosonic commutator form
`S*J*S' = J`. Satisfying these mathematical canonicality conditions is the caller's
responsibility. The constructor checks structural requirements such as modes and dimensions;
it introduces no hidden assumptions, tolerance-based checks, or projections.
"""
function Bogoliubov(modes::AbstractVector{Op}, S::AbstractMatrix)
    lowering_modes = bogoliubov_modes(modes)
    return exact_bogoliubov(
        lowering_modes, bogoliubov_matrix(S, length(lowering_modes)),
    )
end

Bogoliubov(mode::Op, S::AbstractMatrix) = Bogoliubov(Op[mode], S)
Bogoliubov(modes::Tuple{Vararg{Op}}, S::AbstractMatrix) = Bogoliubov(Op[modes...], S)

"""
    Bogoliubov(modes, U, V)

Construct the bosonic map `a ↦ U*a + V*a'`. The implied Nambu matrix is
`[U V; conj(V) conj(U)]`. The blocks are required by contract to satisfy the bosonic
canonical relations; satisfying those mathematical conditions is the caller's responsibility.
The constructor checks only structural requirements such as mode and block dimensions.
"""
function Bogoliubov(
        modes::AbstractVector{Op}, U::AbstractMatrix, V::AbstractMatrix,
    )
    lowering_modes = bogoliubov_modes(modes)
    return exact_bogoliubov(
        lowering_modes,
        bogoliubov_matrix(U, V, length(lowering_modes)),
    )
end

Bogoliubov(mode::Op, U::AbstractMatrix, V::AbstractMatrix) = Bogoliubov(Op[mode], U, V)
Bogoliubov(modes::Tuple{Vararg{Op}}, U::AbstractMatrix, V::AbstractMatrix) =
    Bogoliubov(Op[modes...], U, V)
