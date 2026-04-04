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
