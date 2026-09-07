# Internal affine representation for exact total canonical transformations.
#
# Exact actions are stored as algebra-homogeneous blocks. `UnitaryTransform` keeps the
# compiled rule dictionaries used by `conjugate` and `transform`, but affine metadata is the
# semantic source for inversion and composition.

function AffineBlock(
        structure::AffineStructure, basis::Vector{Op}, linear::AbstractMatrix,
        shift::AbstractVector,
    )
    n = length(basis)
    size(linear) == (n, n) || unitary_error(
        "an affine block on $n generators needs a $n×$n linear map; got $(size(linear))",
    )
    length(shift) == n || unitary_error(
        "an affine block on $n generators needs $n shifts; got $(length(shift))",
    )
    length(Set(basis)) == n || unitary_error("an affine block basis cannot contain duplicates")

    coefficients = Matrix{CNum}(undef, n, n)
    offsets = Vector{CNum}(undef, n)
    for j in 1:n, i in 1:n
        coefficients[i, j] = to_cnum(linear[i, j])
    end
    for i in 1:n
        offsets[i] = to_cnum(shift[i])
    end
    return AffineBlock(structure, copy(basis), coefficients, offsets)
end

function validate_disjoint_blocks(blocks::Vector{AffineBlock})
    seen = Set{Op}()
    for block in blocks, generator in block.basis
        generator in seen && unitary_error(
            "affine blocks cannot overlap on generator `$generator`",
        )
        push!(seen, generator)
    end
    return blocks
end

function AffineAction(
        blocks::Vector{AffineBlock}; relations::Vector{ParamRelation} = ParamRelation[],
    )
    isempty(blocks) && unitary_error("an affine action needs at least one block")
    validate_disjoint_blocks(blocks)
    usable = all(is_usable_rel, relations) ? copy(relations) : filter(is_usable_rel, relations)
    return AffineAction(copy(blocks), usable)
end

function infer_affine_structure(basis::Vector{Op})
    n = length(basis)
    n > 0 || unitary_error("an affine action needs at least one generator")
    if all(is_fock, basis)
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
    elseif all(is_phase_space, basis)
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
    elseif all(o -> is_pauli(o) || is_spin(o), basis)
        return OrthogonalAction()
    elseif all(is_transition, basis)
        return UnitaryLinearAction()
    end
    unitary_error(
        "no single exact affine structure is registered for the supplied generator basis",
    )
end

function AffineAction(
        structure::AffineStructure, basis::Vector{Op}, linear::AbstractMatrix,
        shift::AbstractVector; relations::Vector{ParamRelation} = ParamRelation[],
    )
    return AffineAction(
        AffineBlock[AffineBlock(structure, basis, linear, shift)]; relations = relations,
    )
end

function AffineAction(
        basis::Vector{Op}, linear::AbstractMatrix, shift::AbstractVector;
        relations::Vector{ParamRelation} = ParamRelation[],
    )
    return AffineAction(
        infer_affine_structure(basis), basis, linear, shift; relations = relations,
    )
end

only_affine_block(action::AffineAction) = only(action.blocks)

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
    for j in 1:n, i in 1:n
        source_i = i <= half ? i + half : i - half
        source_j = j <= half ? j + half : j - half
        value = transposed[source_i, source_j]
        (i <= half) == (j <= half) || (value = neg_cnum(value))
        out[i, j] = value
    end
    return out
end

function inverse_linear(linear::Matrix{CNum}, structure::AffineStructure)
    structure === AFFINE_BOSONIC_NAMBU && return inverse_bosonic_nambu(linear)
    structure === AFFINE_SYMPLECTIC_PHASE_SPACE && return inverse_symplectic(linear)
    structure === AFFINE_ORTHOGONAL && return transpose_linear(linear)
    structure === AFFINE_UNITARY_LINEAR && return dagger_linear(linear)
    unitary_error("no exact inverse strategy is registered for affine structure `$structure`")
end

inverse_linear(
    linear::Matrix{CNum}, structure::AffineStructure, ::Vector{ParamRelation},
) = inverse_linear(linear, structure)

function affine_block_rules(block::AffineBlock)
    n = length(block.basis)
    rules = Dict{Op, QAdd}()
    sizehint!(rules, n)
    for i in 1:n
        pairs = Tuple{CNum, Vector{Op}}[]
        sizehint!(pairs, n + 1)
        for j in 1:n
            coefficient = block.linear[i, j]
            iszero_cnum(coefficient) || push!(pairs, (coefficient, Op[block.basis[j]]))
        end
        offset = block.shift[i]
        iszero_cnum(offset) || push!(pairs, (offset, Op[]))
        rules[block.basis[i]] = rule_qadd(pairs)
    end
    return rules
end

function affine_rules(action::AffineAction)
    rules = Dict{Op, QAdd}()
    for block in action.blocks
        merge!(rules, affine_block_rules(block))
    end
    return rules
end

function canonical_block_basis(structure::AffineStructure, raw::Vector{Op})
    unique!(raw)
    if structure === AFFINE_BOSONIC_NAMBU
        lowerings = sort!(unique(lowering.(raw)))
        return vcat(lowerings, adjoint.(lowerings))
    elseif structure === AFFINE_SYMPLECTIC_PHASE_SPACE
        positions = sort!(unique(Op[generator for generator in raw if is_position(generator)]))
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

function blocks_overlap(first::AffineBlock, second::AffineBlock)
    for generator in first.basis
        generator in second.basis && return true
    end
    return false
end

function extend_blocks(blocks::Vector{AffineBlock}, basis::Vector{Op})
    n = length(basis)
    locations = Dict{Op, Int}(generator => i for (i, generator) in enumerate(basis))
    linear = fill(CNUM_ZERO, n, n)
    shift = fill(CNUM_ZERO, n)
    for i in 1:n
        linear[i, i] = CNUM_ONE
    end
    for block in blocks
        for (source_row, generator) in enumerate(block.basis)
            target_row = locations[generator]
            for column in 1:n
                linear[target_row, column] = CNUM_ZERO
            end
            for (source_column, source_generator) in enumerate(block.basis)
                target_column = locations[source_generator]
                linear[target_row, target_column] = block.linear[source_row, source_column]
            end
            shift[target_row] = block.shift[source_row]
        end
    end
    return linear, shift
end

function subset_positions(subset::Vector{Op}, superset::Vector{Op})::Union{Nothing, Vector{Int}}
    positions = Vector{Int}(undef, length(subset))
    for i in eachindex(subset)
        position = findfirst(==(subset[i]), superset)
        position === nothing && return nothing
        positions[i] = position
    end
    return positions
end

function is_identity_affine_linear(linear::Matrix{CNum})
    n, m = size(linear)
    n == m || return false
    for j in 1:n, i in 1:n
        coefficient = linear[i, j]
        if i == j
            isequal(coefficient, CNUM_ONE) || return false
        else
            iszero_cnum(coefficient) || return false
        end
    end
    return true
end

function is_zero_affine_shift(shift::Vector{CNum})
    for coefficient in shift
        iszero_cnum(coefficient) || return false
    end
    return true
end

function compose_subset_block(
        first::AffineBlock, second::AffineBlock, positions::Vector{Int},
        relations::Vector{ParamRelation},
    )
    n = length(first.basis)
    m = length(second.basis)

    linear = if is_identity_affine_linear(second.linear)
        first.linear
    else
        out = copy(first.linear)
        scratch = ParamRelation[]
        for local_j in 1:m, i in 1:n
            value = CNUM_ZERO
            for local_k in 1:m
                value = add_cnum(
                    value,
                    mul_cnum(
                        first.linear[i, positions[local_k]], second.linear[local_k, local_j],
                    ),
                )
            end
            out[i, positions[local_j]] = reduce_affine(value, relations, scratch)
        end
        out
    end

    shift = if is_zero_affine_shift(second.shift)
        first.shift
    else
        out = Vector{CNum}(undef, n)
        scratch = ParamRelation[]
        for i in 1:n
            value = first.shift[i]
            for local_k in 1:m
                value = add_cnum(
                    value,
                    mul_cnum(first.linear[i, positions[local_k]], second.shift[local_k]),
                )
            end
            out[i] = reduce_affine(value, relations, scratch)
        end
        out
    end

    return AffineBlock(first.structure, first.basis, linear, shift)
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

    # `first * second` agrees with sequential conjugation: substitute the second map into
    # the first image, A₁(A₂z + b₂) + b₁.
    for j in 1:n, i in 1:n
        value = CNUM_ZERO
        for k in 1:n
            value = add_cnum(value, mul_cnum(first_linear[i, k], second_linear[k, j]))
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

function compose_overlapping_block(
        first::AffineBlock, second::AffineBlock, relations::Vector{ParamRelation},
    )
    structure = second.structure
    first.structure === structure || unitary_error(
        "overlapping affine blocks must describe the same canonical algebra",
    )

    positions = subset_positions(second.basis, first.basis)
    positions === nothing || return compose_subset_block(first, second, positions, relations)

    raw = copy(second.basis)
    append!(raw, first.basis)
    basis = canonical_block_basis(structure, raw)
    first_linear, first_shift = extend_blocks(AffineBlock[first], basis)
    second_linear, second_shift = extend_blocks(AffineBlock[second], basis)
    linear, shift = compose_affine_data(
        first_linear, first_shift, second_linear, second_shift, relations,
    )
    return AffineBlock(structure, basis, linear, shift)
end

function compose_overlapping_blocks(
        first_blocks::Vector{AffineBlock}, second::AffineBlock,
        relations::Vector{ParamRelation},
    )
    length(first_blocks) == 1 &&
        return compose_overlapping_block(only(first_blocks), second, relations)

    structure = second.structure
    all(block -> block.structure === structure, first_blocks) || unitary_error(
        "overlapping affine blocks must describe the same canonical algebra",
    )
    raw = copy(second.basis)
    for block in first_blocks
        append!(raw, block.basis)
    end
    basis = canonical_block_basis(structure, raw)
    first_linear, first_shift = extend_blocks(first_blocks, basis)
    second_linear, second_shift = extend_blocks(AffineBlock[second], basis)
    linear, shift = compose_affine_data(
        first_linear, first_shift, second_linear, second_shift, relations,
    )
    return AffineBlock(structure, basis, linear, shift)
end

function compose_action_metadata(
        first::AffineAction, second::AffineAction, relations::Vector{ParamRelation},
    )
    result = copy(first.blocks)
    for second_block in second.blocks
        overlapping = findall(block -> blocks_overlap(block, second_block), result)
        if isempty(overlapping)
            push!(result, second_block)
            continue
        end
        if length(overlapping) == 1
            index = only(overlapping)
            composed = compose_overlapping_block(result[index], second_block, relations)
            deleteat!(result, index)
            push!(result, composed)
            continue
        end
        first_blocks = result[overlapping]
        composed = compose_overlapping_blocks(first_blocks, second_block, relations)
        deleteat!(result, reverse(overlapping))
        push!(result, composed)
    end
    return AffineAction(result; relations = relations)
end