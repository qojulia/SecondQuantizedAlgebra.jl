# === General exact bosonic Bogoliubov transformations ===

function bogoliubov_modes(modes::AbstractVector{Op})
    isempty(modes) && unitary_error("`Bogoliubov` needs at least one Fock mode")
    lowering_modes = Op[]
    sizehint!(lowering_modes, length(modes))
    seen = Set{SiteKey}()
    for mode in modes
        lowering = fock_or_throw(mode, "`Bogoliubov`")
        key = site_key(lowering)
        key in seen && unitary_error(
            "`Bogoliubov` received the same Fock mode more than once: `$lowering`",
        )
        push!(seen, key)
        push!(lowering_modes, lowering)
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

function bogoliubov_zero(c::CNum, scratch::Vector{ParamRelation})
    reduced = reduce_all(c, ParamRelation[], true, scratch)
    return iszero_cnum(reduced)
end

function validate_bogoliubov_action(action::AffineAction)
    action.structure === AFFINE_BOSONIC_NAMBU || unitary_error(
        "internal Bogoliubov validation requires a bosonic Nambu affine action",
    )
    n = length(action.basis)
    half = n ÷ 2
    scratch = ParamRelation[]

    for i in 1:half
        for j in 1:half
            residual = add_cnum(
                action.linear[half + i, half + j],
                neg_cnum(conj_cnum(action.linear[i, j])),
            )
            bogoliubov_zero(residual, scratch) || unitary_error(
                "`Bogoliubov` Nambu matrix does not preserve adjoints",
            )
            residual = add_cnum(
                action.linear[half + i, j],
                neg_cnum(conj_cnum(action.linear[i, half + j])),
            )
            bogoliubov_zero(residual, scratch) || unitary_error(
                "`Bogoliubov` Nambu matrix does not preserve adjoints",
            )
        end
    end

    inverse = inverse_linear(action.linear, action.structure, action.relations)
    for (left, right) in ((action.linear, inverse), (inverse, action.linear))
        for j in 1:n, i in 1:n
            residual = i == j ? CNUM_NEG1 : CNUM_ZERO
            for k in 1:n
                residual = add_cnum(residual, mul_cnum(left[i, k], right[k, j]))
            end
            bogoliubov_zero(residual, scratch) || unitary_error(
                "`Bogoliubov` matrix is not canonical: residual ($i, $j) is " *
                    "`$(to_num(reduce_all(residual, ParamRelation[], true, scratch)))`",
            )
        end
    end
    return action
end

function exact_bogoliubov(modes::Vector{Op}, matrix::Matrix{CNum})
    basis = bogoliubov_basis(modes)
    action = AffineAction(
        BosonicNambu(), basis, matrix, fill(CNUM_ZERO, length(basis)),
    )
    validate_bogoliubov_action(action)
    return canonical_transform(action)
end

"""
    Bogoliubov(modes, S)

Construct an exact bosonic Bogoliubov transformation in Nambu ordering
`(a₁, …, aₙ, a₁', …, aₙ')`. `S` must preserve both adjoints and the bosonic
commutator form exactly. Symbolic matrices whose canonicality cannot be proven are rejected.
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

Construct the exact bosonic map `a ↦ U*a + V*a'`. The implied Nambu matrix is
`[U V; conj(V) conj(U)]` and must satisfy the bosonic canonical relations exactly.
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
