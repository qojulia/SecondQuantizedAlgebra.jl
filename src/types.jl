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

"""
    LazyOrder <: OrderingConvention

Lazy ordering: no rules applied during `*`. Operator algebra rules are deferred
until `simplify()` or `normal_order()` is called explicitly.
"""
struct LazyOrder <: OrderingConvention end
