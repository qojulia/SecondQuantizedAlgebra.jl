# Proof-preserving helpers for affine actions that are canonical by construction.
#
# Named constructors and exact generated flows establish canonicality structurally before
# reaching this layer. These helpers preserve that proof while compiling and inverting the
# affine metadata carried by `UnitaryTransform`.

function canonical_affine_inverse(action::AffineAction)
    linear = inverse_linear(action.linear, action.structure, action.relations)
    n = length(action.basis)
    shift = Vector{CNum}(undef, n)
    scratch = ParamRelation[]
    for i in 1:n
        value = CNUM_ZERO
        for j in 1:n
            value = add_cnum(value, mul_cnum(linear[i, j], action.shift[j]))
        end
        shift[i] = reduce_affine(neg_cnum(value), action.relations, scratch)
    end
    return AffineAction(
        action.structure, action.basis, linear, shift; relations = action.relations,
    )
end

function canonical_transform(action::AffineAction)
    inverse_action = canonical_affine_inverse(action)
    return validated_transform(
        affine_rules(action), affine_rules(inverse_action), zero_qadd(), StaticTime(),
        action.relations, action,
    )
end

compiled_inverse_action_metadata(action::AffineAction) = canonical_affine_inverse(action)

# Composition correctness is defined by the affine metadata product. The forward and inverse
# rule dictionaries already stored on the operands are compiled execution IR, so let
# `compose` reuse them through its rule-composition fallback instead of reconstructing the
# inverse affine map and recompiling both dictionaries on every `*`.
compile_composed_action_metadata(::AffineAction) = nothing
