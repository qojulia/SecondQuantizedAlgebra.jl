"""
    QField

Abstract type representing any expression involving operators.
"""
abstract type QField end

"""
    QSym <: QField

Abstract type representing fundamental operator types (leaves in the expression tree).
"""
abstract type QSym <: QField end

"""
    QTerm <: QField

Abstract type representing compound noncommutative expressions.
"""
abstract type QTerm <: QField end

"""
    OrderingConvention

Abstract type for operator ordering conventions used by `simplify`.
"""
abstract type OrderingConvention end

"""
    NormalOrder <: OrderingConvention

Normal ordering: creation operators to the left of annihilation operators.
Applies `[a, a†] = 1` for Fock and `[Sⱼ, Sₖ] = iϵⱼₖₗSₗ` for Spin.
"""
struct NormalOrder <: OrderingConvention end

# TODO: SymmetricOrder for TWA (Truncated Wigner Approximation)
# Symmetric (Weyl) ordering places operators in symmetrized form: (ab + ba)/2.
# To implement:
#   struct SymmetricOrder <: OrderingConvention end
# Then add a method for _apply_ordering_rule(a, b, same_space, i, arg_c, ops, ::SymmetricOrder)
# in simplify.jl. The ordering-independent reductions (Transition, Pauli) still apply;
# only the Fock and Spin commutation swap rules change.
