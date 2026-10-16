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

function validated_action(blocks::Vector{AffineBlock}, relations::Vector{ParamRelation})
    isempty(blocks) && unitary_error("an affine action needs at least one block")
    validate_disjoint_blocks(blocks)
    usable = all(is_usable_rel, relations) ? relations : filter(is_usable_rel, relations)
    return AffineAction(blocks, usable)
end

function operator_affine_structure(o::Op)
    is_fock(o) && return AFFINE_BOSONIC_NAMBU
    is_phase_space(o) && return AFFINE_SYMPLECTIC_PHASE_SPACE
    (is_pauli(o) || is_spin(o)) && return AFFINE_ORTHOGONAL
    is_transition(o) && return AFFINE_UNITARY_LINEAR
    return AFFINE_NONE
end

function infer_affine_structure(basis::Vector{Op})
    n = length(basis)
    n > 0 || unitary_error("an affine action needs at least one generator")
    structure = operator_affine_structure(first(basis))
    (structure !== AFFINE_NONE && all(o -> operator_affine_structure(o) === structure, basis)) ||
        unitary_error("no single exact affine structure is registered for the supplied generator basis")
    half = n ÷ 2
    if structure === AFFINE_BOSONIC_NAMBU
        iseven(n) || unitary_error("a bosonic Nambu basis needs an even number of generators")
        for i in 1:half
            is_destroy(basis[i]) || unitary_error(
                "bosonic Nambu ordering requires annihilation operators first",
            )
            basis[half + i] == adjoint(basis[i]) || unitary_error(
                "bosonic Nambu ordering requires matching creation operators second",
            )
        end
    elseif structure === AFFINE_SYMPLECTIC_PHASE_SPACE
        iseven(n) || unitary_error("a phase-space basis needs an even number of generators")
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
    end
    return structure
end

function AffineAction(
        basis::Vector{Op}, linear::Matrix{CNum}, shift::Vector{CNum},
        relations::Vector{ParamRelation} = ParamRelation[],
    )
    structure = infer_affine_structure(basis)
    n = length(basis)
    size(linear) == (n, n) || unitary_error(
        "an affine block on $n generators needs a $n×$n linear map; got $(size(linear))",
    )
    length(shift) == n || unitary_error(
        "an affine block on $n generators needs $n shifts; got $(length(shift))",
    )
    length(Set(basis)) == n || unitary_error("an affine block basis cannot contain duplicates")
    return validated_action(AffineBlock[AffineBlock(structure, basis, linear, shift)], relations)
end

function reduce_affine(c::CNum, relations::Vector{ParamRelation}, scratch::Vector{ParamRelation})
    isempty(relations) && return c
    return reduce_all(c, relations, true, scratch)
end

dagger_linear(linear::Matrix{CNum}) = Matrix{CNum}(linear')
transpose_linear(linear::Matrix{CNum}) = Matrix{CNum}(transpose(linear))

function inverse_bosonic_nambu(linear::Matrix{CNum})
    n = size(linear, 1)
    iseven(n) || unitary_error("a bosonic Nambu action needs an even-dimensional basis")
    half = n ÷ 2
    out = dagger_linear(linear)
    for j in 1:n, i in 1:n
        left_sign = i <= half ? 1 : -1
        right_sign = j <= half ? 1 : -1
        left_sign == right_sign || (out[i, j] = -out[i, j])
    end
    return out
end

function inverse_symplectic(linear::Matrix{CNum})
    n = size(linear, 1)
    iseven(n) || unitary_error("a phase-space action needs an even-dimensional basis")
    half = n ÷ 2
    out = Matrix{CNum}(undef, n, n)
    for j in 1:n, i in 1:n
        source_i = i <= half ? i + half : i - half
        source_j = j <= half ? j + half : j - half
        value = linear[source_j, source_i]
        (i <= half) == (j <= half) || (value = -value)
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

function affine_rules(action::AffineAction)
    rules = Dict{Op, QAdd}()
    sizehint!(rules, sum(length(block.basis) for block in action.blocks))
    for block in action.blocks
        n = length(block.basis)
        for i in 1:n
            terms = QTermDict()
            sizehint!(terms, n + 1)
            for j in 1:n
                coefficient = block.linear[i, j]
                iszero(coefficient) || addto!(terms, Op[block.basis[j]], coefficient)
            end
            offset = block.shift[i]
            iszero(offset) || addto!(terms, Op[], offset)
            rules[block.basis[i]] = QAdd(terms, EMPTY_INDICES)
        end
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
            found = position_of(
                p -> is_momentum(p) && site_key(p) == site_key(x), raw,
            )
            iszero(found) && unitary_error(
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

function identity_linear(n::Int)
    linear = fill(CNUM_ZERO, n, n)
    for i in 1:n
        linear[i, i] = CNUM_ONE
    end
    return linear
end

function extend_blocks(blocks::Vector{AffineBlock}, basis::Vector{Op})
    n = length(basis)
    locations = Dict{Op, Int}(generator => i for (i, generator) in enumerate(basis))
    linear = identity_linear(n)
    shift = fill(CNUM_ZERO, n)
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

subset_positions(subset::Vector{Op}, superset::Vector{Op}) =
    Int[position_of(==(operator), superset) for operator in subset]

function is_identity_affine_linear(linear::Matrix{CNum})
    n, m = size(linear)
    n == m || return false
    for j in 1:n, i in 1:n
        coefficient = linear[i, j]
        if i == j
            isequal(coefficient, CNUM_ONE) || return false
        else
            iszero(coefficient) || return false
        end
    end
    return true
end

function is_zero_affine_shift(shift::Vector{CNum})
    for coefficient in shift
        iszero(coefficient) || return false
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
                value = value + first.linear[i, positions[local_k]] * second.linear[local_k, local_j]
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
                value = value + first.linear[i, positions[local_k]] * second.shift[local_k]
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

    for j in 1:n, i in 1:n
        value = CNUM_ZERO
        for k in 1:n
            value = value + first_linear[i, k] * second_linear[k, j]
        end
        linear[i, j] = reduce_affine(value, relations, scratch)
    end
    for i in 1:n
        value = first_shift[i]
        for k in 1:n
            value = value + first_linear[i, k] * second_shift[k]
        end
        shift[i] = reduce_affine(value, relations, scratch)
    end
    return linear, shift
end

function compose_overlapping(
        first_blocks::Vector{AffineBlock}, second::AffineBlock, relations::Vector{ParamRelation},
    )
    structure = second.structure
    all(block -> block.structure === structure, first_blocks) || unitary_error(
        "overlapping affine blocks must describe the same canonical algebra",
    )
    if length(first_blocks) == 1
        first = only(first_blocks)
        positions = subset_positions(second.basis, first.basis)
        all(!iszero, positions) && return compose_subset_block(first, second, positions, relations)
    end
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

function compose_action(
        first::AffineAction, second::AffineAction, relations::Vector{ParamRelation},
    )
    result = copy(first.blocks)
    for block in second.blocks
        overlapping = findall(existing -> blocks_overlap(existing, block), result)
        if isempty(overlapping)
            push!(result, block)
        else
            composed = compose_overlapping(result[overlapping], block, relations)
            deleteat!(result, overlapping)
            push!(result, composed)
        end
    end
    return AffineAction(result, relations)
end

function canonical_block_inverse(block::AffineBlock, relations::Vector{ParamRelation})
    linear = inverse_linear(block.linear, block.structure)
    n = length(block.basis)
    shift = Vector{CNum}(undef, n)
    scratch = ParamRelation[]
    for i in 1:n
        value = CNUM_ZERO
        for j in 1:n
            value = value + linear[i, j] * block.shift[j]
        end
        shift[i] = reduce_affine(-value, relations, scratch)
    end
    return AffineBlock(block.structure, block.basis, linear, shift)
end

function canonical_affine_inverse(action::AffineAction)
    blocks = AffineBlock[
        canonical_block_inverse(block, action.relations) for block in action.blocks
    ]
    return AffineAction(blocks, action.relations)
end

canonical_transform(action::AffineAction) =
    validated_transform(action, zero_qadd(), StaticTime())

function substitute_relation(r::ParamRelation, rules::AbstractDict)
    hi = Symbolics.substitute(r.hi, rules)
    lo = Symbolics.substitute(r.lo, rules)
    return ParamRelation(hi, lo, r.sign)
end

function substitute_affine_block(block::AffineBlock, rules::AbstractDict)
    n = length(block.basis)
    linear = Matrix{CNum}(undef, n, n)
    shift = Vector{CNum}(undef, n)
    for j in 1:n, i in 1:n
        linear[i, j] = substitute_cnum(block.linear[i, j], rules)
    end
    for i in 1:n
        shift[i] = substitute_cnum(block.shift[i], rules)
    end
    return AffineBlock(block.structure, block.basis, linear, shift)
end

function substitute_affine_action(action::AffineAction, rules::AbstractDict)
    blocks = AffineBlock[substitute_affine_block(block, rules) for block in action.blocks]
    relations = ParamRelation[substitute_relation(r, rules) for r in action.relations]
    return validated_action(blocks, relations)
end

substitute_time(time::StaticTime, ::AbstractDict) = time
function substitute_time(time::DynamicTime, rules::AbstractDict)
    substituted = Symbolics.substitute(time.variable, rules)
    isequal(substituted, time.variable) || unitary_error(
        "substitution of the differentiation variable `$(time.variable)` is not supported",
    )
    return time
end

function validate_transform_substitution_rules(rules::AbstractDict)
    for key in keys(rules)
        key isa QSym && unitary_error(
            "unitary-transform substitution resolves scalar parameters only; " *
                "operator substitutions belong on the transformed expression",
        )
    end
    return rules
end

function conjugation_closed_substitutions(rules::AbstractDict)
    closed = ScalarRules()
    sizehint!(closed, 2 * length(rules))
    for (from, to) in rules
        add_scalar_rule!(closed, from, to)
    end
    for (from, to) in rules
        add_conjugate_rule!(closed, from, to)
    end
    return closed
end

add_conjugate_rule!(closed::ScalarRules, from, to) = closed
function add_conjugate_rule!(closed::ScalarRules, from::Union{Num, RawExpression}, to::Union{Number, RawExpression})
    raw_from = raw_rule(from)
    conjugated_from = raw_conj(raw_from)
    isequal(conjugated_from, raw_from) && return closed
    closed[conjugated_from] = raw_rule(conj(to))
    return closed
end

function SymbolicUtils.substitute(U::UnitaryTransform{T}, rules::AbstractDict) where {T}
    validate_transform_substitution_rules(rules)
    isempty(rules) && return U
    closed_rules = conjugation_closed_substitutions(rules)
    time = substitute_time(U.time, closed_rules)
    action = substitute_affine_action(U.action, closed_rules)
    gauge = substitute(U.gauge, closed_rules)
    return validated_transform(action, gauge, time)
end
