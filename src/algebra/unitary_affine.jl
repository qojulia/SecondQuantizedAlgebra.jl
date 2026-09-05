# Internal affine representation for exact total canonical transformations.
#
# `AffineAction` is the construction IR. `UnitaryTransform` remains the compiled execution
# representation used by `conjugate` and `transform`.
@enum AffineStructure::UInt8 begin
    AFFINE_GENERIC
    AFFINE_BOSONIC_NAMBU
    AFFINE_SYMPLECTIC_PHASE_SPACE
    AFFINE_ORTHOGONAL
    AFFINE_UNITARY_LINEAR
end

# Keep the descriptive constructor spellings at call sites while storing the structure as a
# value tag. This mirrors `OpKind`: algebra role is runtime data, not part of the Julia type.
GenericAffine() = AFFINE_GENERIC
BosonicNambu() = AFFINE_BOSONIC_NAMBU
SymplecticPhaseSpace() = AFFINE_SYMPLECTIC_PHASE_SPACE
OrthogonalAction() = AFFINE_ORTHOGONAL
UnitaryLinearAction() = AFFINE_UNITARY_LINEAR

struct AffineAction
    structure::AffineStructure
    basis::Vector{Op}
    linear::Matrix{CNum}
    shift::Vector{CNum}
    relations::Vector{ParamRelation}
end

function AffineAction(
        structure::AffineStructure, basis::Vector{Op}, linear::AbstractMatrix,
        shift::AbstractVector;
        relations::Vector{ParamRelation} = ParamRelation[],
    )
    n = length(basis)
    size(linear) == (n, n) || unitary_error(
        "an affine action on $n generators needs a $n×$n linear map; got $(size(linear))",
    )
    length(shift) == n || unitary_error(
        "an affine action on $n generators needs $n shifts; got $(length(shift))",
    )
    length(Set(basis)) == n || unitary_error("an affine action basis cannot contain duplicates")

    coefficients = Matrix{CNum}(undef, n, n)
    offsets = Vector{CNum}(undef, n)
    for j in 1:n, i in 1:n
        coefficients[i, j] = to_cnum(linear[i, j])
    end
    for i in 1:n
        offsets[i] = to_cnum(shift[i])
    end
    return AffineAction(
        structure, copy(basis), coefficients, offsets, copy(relations),
    )
end

function infer_affine_structure(basis::Vector{Op})
    n = length(basis)
    if n > 0 && all(is_fock, basis)
        iseven(n) || unitary_error("a bosonic Nambu basis needs an even number of generators")
        half = n ÷ 2
        for i in 1:half
            is_destroy(basis[i]) || unitary_error(
                "bosonic Nambu ordering requires annihilation operators first",
            )
            basis[half + i] == adjoint(basis[i]) || unitary_error(
                "bosonic Nambu ordering requires matching creation operators second",
            )
        end
        return BosonicNambu()
    elseif n > 0 && all(is_phase_space, basis)
        iseven(n) || unitary_error("a phase-space basis needs an even number of generators")
        half = n ÷ 2
        for i in 1:half
            is_position(basis[i]) || unitary_error(
                "phase-space ordering requires position operators first",
            )
            is_momentum(basis[half + i]) || unitary_error(
                "phase-space ordering requires momentum operators second",
            )
            site_key(basis[i]) == site_key(basis[half + i]) || unitary_error(
                "phase-space basis must pair position and momentum operators by site",
            )
        end
        return SymplecticPhaseSpace()
    elseif n > 0 && all(o -> is_pauli(o) || is_spin(o), basis)
        return OrthogonalAction()
    elseif n > 0 && all(is_transition, basis)
        return UnitaryLinearAction()
    end
    return GenericAffine()
end

AffineAction(basis::Vector{Op}, linear::AbstractMatrix, shift::AbstractVector; kwargs...) =
    AffineAction(infer_affine_structure(basis), basis, linear, shift; kwargs...)

function reduce_affine(c::CNum, relations::Vector{ParamRelation}, scratch::Vector{ParamRelation})
    isempty(relations) && return c
    return reduce_all(c, relations, true, scratch)
end

function dagger_linear(linear::Matrix{CNum})
    n, m = size(linear)
    out = Matrix{CNum}(undef, m, n)
    for j in 1:m, i in 1:n
        out[j, i] = conj_cnum(linear[i, j])
    end
    return out
end

function transpose_linear(linear::Matrix{CNum})
    n, m = size(linear)
    out = Matrix{CNum}(undef, m, n)
    for j in 1:m, i in 1:n
        out[j, i] = linear[i, j]
    end
    return out
end

function inverse_bosonic_nambu(linear::Matrix{CNum})
    n = size(linear, 1)
    iseven(n) || unitary_error("a bosonic Nambu action needs an even-dimensional basis")
    half = n ÷ 2
    out = dagger_linear(linear)
    for j in 1:n, i in 1:n
        left_sign = i <= half ? 1 : -1
        right_sign = j <= half ? 1 : -1
        left_sign == right_sign || (out[i, j] = neg_cnum(out[i, j]))
    end
    return out
end

function inverse_symplectic(linear::Matrix{CNum})
    n = size(linear, 1)
    iseven(n) || unitary_error("a phase-space action needs an even-dimensional basis")
    half = n ÷ 2
    transposed = transpose_linear(linear)
    out = Matrix{CNum}(undef, n, n)
    # Ω⁻¹ Aᵀ Ω for Ω = [0 I; -I 0].
    for j in 1:n, i in 1:n
        source_i = i <= half ? i + half : i - half
        source_j = j <= half ? j + half : j - half
        value = transposed[source_i, source_j]
        (i <= half) == (j <= half) || (value = neg_cnum(value))
        out[i, j] = value
    end
    return out
end

function inverse_generic(
        linear::Matrix{CNum}, relations::Vector{ParamRelation},
    )
    n = size(linear, 1)
    augmented = Matrix{CNum}(undef, n, 2n)
    for j in 1:n, i in 1:n
        augmented[i, j] = linear[i, j]
        augmented[i, n + j] = i == j ? CNUM_ONE : CNUM_ZERO
    end

    scratch = ParamRelation[]
    for column in 1:n
        pivot_offset = findfirst(
            row -> !iszero_cnum(reduce_affine(augmented[row, column], relations, scratch)),
            column:n,
        )
        pivot_offset === nothing && unitary_error("affine linear map is singular")
        pivot = column - 1 + pivot_offset
        if pivot != column
            for j in 1:(2n)
                augmented[column, j], augmented[pivot, j] =
                    augmented[pivot, j], augmented[column, j]
            end
        end

        pivot_coefficient = reduce_affine(augmented[column, column], relations, scratch)
        pivot_inverse = inv(pivot_coefficient)
        for j in 1:(2n)
            augmented[column, j] = reduce_affine(
                mul_cnum(augmented[column, j], pivot_inverse), relations, scratch,
            )
        end

        for row in 1:n
            row == column && continue
            factor = reduce_affine(augmented[row, column], relations, scratch)
            iszero_cnum(factor) && continue
            for j in 1:(2n)
                augmented[row, j] = reduce_affine(
                    add_cnum(
                        augmented[row, j],
                        neg_cnum(mul_cnum(factor, augmented[column, j])),
                    ),
                    relations,
                    scratch,
                )
            end
        end
    end

    inverse = Matrix{CNum}(undef, n, n)
    for j in 1:n, i in 1:n
        inverse[i, j] = reduce_affine(augmented[i, n + j], relations, scratch)
    end
    return inverse
end

function inverse_linear(
        linear::Matrix{CNum}, structure::AffineStructure,
        relations::Vector{ParamRelation},
    )
    structure === AFFINE_BOSONIC_NAMBU && return inverse_bosonic_nambu(linear)
    structure === AFFINE_SYMPLECTIC_PHASE_SPACE && return inverse_symplectic(linear)
    structure === AFFINE_ORTHOGONAL && return transpose_linear(linear)
    structure === AFFINE_UNITARY_LINEAR && return dagger_linear(linear)
    return inverse_generic(linear, relations)
end

function affine_rules(action::AffineAction)
    n = length(action.basis)
    rules = Dict{Op, QAdd}()
    sizehint!(rules, n)
    for i in 1:n
        pairs = Tuple{CNum, Vector{Op}}[]
        sizehint!(pairs, n + 1)
        for j in 1:n
            coefficient = action.linear[i, j]
            iszero_cnum(coefficient) || push!(pairs, (coefficient, Op[action.basis[j]]))
        end
        offset = action.shift[i]
        iszero_cnum(offset) || push!(pairs, (offset, Op[]))
        rules[action.basis[i]] = rule_qadd(pairs)
    end
    return rules
end

function affine_union_basis(first::AffineAction, second::AffineAction)
    raw = copy(first.basis)
    seen = Set(raw)
    for generator in second.basis
        generator in seen && continue
        push!(raw, generator)
        push!(seen, generator)
    end

    if all(is_fock, raw)
        lowerings = Op[]
        lowering_seen = Set{Op}()
        for generator in raw
            d = lowering(generator)
            d in lowering_seen && continue
            push!(lowerings, d)
            push!(lowering_seen, d)
        end
        sort!(lowerings)
        return vcat(lowerings, adjoint.(lowerings))
    elseif all(is_phase_space, raw)
        positions = sort!(Op[generator for generator in raw if is_position(generator)])
        momenta = Op[]
        for x in positions
            found = findfirst(
                p -> is_momentum(p) && site_key(p) == site_key(x), raw,
            )
            found === nothing && unitary_error(
                "phase-space affine composition lost the momentum paired with `$x`",
            )
            push!(momenta, raw[found])
        end
        return vcat(positions, momenta)
    end

    sort!(raw)
    return raw
end

function extend_affine(action::AffineAction, basis::Vector{Op})
    n = length(basis)
    locations = Dict{Op, Int}(generator => i for (i, generator) in enumerate(basis))
    linear = fill(CNUM_ZERO, n, n)
    shift = fill(CNUM_ZERO, n)
    for i in 1:n
        linear[i, i] = CNUM_ONE
    end

    for (source_row, generator) in enumerate(action.basis)
        target_row = locations[generator]
        for column in 1:n
            linear[target_row, column] = CNUM_ZERO
        end
        for (source_column, source_generator) in enumerate(action.basis)
            target_column = locations[source_generator]
            linear[target_row, target_column] = action.linear[source_row, source_column]
        end
        shift[target_row] = action.shift[source_row]
    end
    return linear, shift
end

function compose_affine_data(
        first_linear::Matrix{CNum}, first_shift::Vector{CNum},
        second_linear::Matrix{CNum}, second_shift::Vector{CNum},
        relations::Vector{ParamRelation},
    )
    n = length(first_shift)
    linear = Matrix{CNum}(undef, n, n)
    shift = Vector{CNum}(undef, n)
    scratch = ParamRelation[]

    # `first * second` must agree with sequential conjugation: apply the first
    # rule image and then substitute the second rules into it. For column vectors
    # this is A₁(A₂z + b₂) + b₁.
    for j in 1:n, i in 1:n
        value = CNUM_ZERO
        for k in 1:n
            value = add_cnum(
                value, mul_cnum(first_linear[i, k], second_linear[k, j]),
            )
        end
        linear[i, j] = reduce_affine(value, relations, scratch)
    end

    for i in 1:n
        value = first_shift[i]
        for k in 1:n
            value = add_cnum(value, mul_cnum(first_linear[i, k], second_shift[k]))
        end
        shift[i] = reduce_affine(value, relations, scratch)
    end
    return linear, shift
end

function compose_action_metadata(
        first::AffineAction, second::AffineAction, relations::Vector{ParamRelation},
    )
    basis = affine_union_basis(first, second)
    first_linear, first_shift = extend_affine(first, basis)
    second_linear, second_shift = extend_affine(second, basis)
    linear, shift = compose_affine_data(
        first_linear, first_shift, second_linear, second_shift, relations,
    )
    structure = infer_affine_structure(basis)
    return AffineAction(structure, basis, linear, shift; relations = relations)
end
