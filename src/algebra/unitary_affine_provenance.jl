# Proof-preserving helpers for affine actions that are canonical by construction or by API
# contract. No generic symbolic matrix inversion is used: each homogeneous affine block has a
# structural inverse formula determined by its algebra.

function canonical_block_inverse(block::AffineBlock, relations::Vector{ParamRelation})
    linear = inverse_linear(block.linear, block.structure)
    n = length(block.basis)
    shift = Vector{CNum}(undef, n)
    scratch = ParamRelation[]
    for i in 1:n
        value = CNUM_ZERO
        for j in 1:n
            value = add_cnum(value, mul_cnum(linear[i, j], block.shift[j]))
        end
        shift[i] = reduce_affine(neg_cnum(value), relations, scratch)
    end
    return AffineBlock(block.structure, block.basis, linear, shift)
end

function canonical_affine_inverse(action::AffineAction)
    blocks = AffineBlock[
        canonical_block_inverse(block, action.relations) for block in action.blocks
    ]
    # Structural inversion preserves block disjointness and does not mutate parameter relations.
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
    return AffineAction(blocks; relations = relations)
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
    return nothing
end

# Keep scalar substitutions closed under conjugation so exact rational/complex values replace
# both a symbolic coefficient and its conjugate before the affine map is recompiled.
function conjugation_closed_substitutions(rules::AbstractDict)
    isempty(rules) && return rules
    closed = Dict{Any, Any}(rules)
    sizehint!(closed, 2 * length(rules))
    for (from, to) in rules
        raw_from = if from isa Num
            SymbolicUtils.unwrap(from)
        elseif from isa SymbolicUtils.BasicSymbolic
            from
        else
            continue
        end
        conjugated_from = raw_conj(raw_from)
        isequal(conjugated_from, raw_from) && continue
        conjugated_to = if to isa Number || to isa Num || to isa SymbolicUtils.BasicSymbolic
            conj(to)
        else
            continue
        end
        closed[conjugated_from] = conjugated_to
    end
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
