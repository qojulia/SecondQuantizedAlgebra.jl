module ExactResponseSolverPrototype

using SecondQuantizedAlgebra
using LinearAlgebra
using Symbolics

import SecondQuantizedAlgebra: CNum, _CNUM_ZERO, _CNUM_ONE, _add_cnum,
    _mul_cnum, _neg_cnum, _to_cnum, _iszero_cnum, to_num

export FLResult, faddeev_leverrier, solve_adjugate, solve_2x2,
    solve_structural_diagonal, multiply, residual_certificate,
    coefficient_size, coefficient_nodes, coefficient_terms, matrix_size

"""Disposable certificate for `adj(A)` and `det(A)` computed without symbolic pivots."""
struct FLResult
    determinant::CNum
    adjugate::Matrix{CNum}
end

function identity_cnum(n::Int)
    result = fill(_CNUM_ZERO, n, n)
    @inbounds for i in 1:n
        result[i, i] = _CNUM_ONE
    end
    return result
end

function multiply(a::Matrix{CNum}, b::Matrix{CNum})
    size(a, 2) == size(b, 1) || throw(DimensionMismatch())
    result = fill(_CNUM_ZERO, size(a, 1), size(b, 2))
    @inbounds for j in axes(b, 2), k in axes(a, 2), i in axes(a, 1)
        result[i, j] = _add_cnum(result[i, j], _mul_cnum(a[i, k], b[k, j]))
    end
    return result
end

function multiply(a::Matrix{CNum}, b::Vector{CNum})
    size(a, 2) == length(b) || throw(DimensionMismatch())
    result = fill(_CNUM_ZERO, size(a, 1))
    @inbounds for k in axes(a, 2), i in axes(a, 1)
        result[i] = _add_cnum(result[i], _mul_cnum(a[i, k], b[k]))
    end
    return result
end

function add_diagonal!(a::Matrix{CNum}, c::CNum)
    @inbounds for i in axes(a, 1)
        a[i, i] = _add_cnum(a[i, i], c)
    end
    return a
end

function trace_cnum(a::Matrix{CNum})
    result = _CNUM_ZERO
    @inbounds for i in axes(a, 1)
        result = _add_cnum(result, a[i, i])
    end
    return result
end

function scale(a::Matrix{CNum}, c::CNum)
    result = similar(a)
    @inbounds for i in eachindex(a)
        result[i] = _mul_cnum(c, a[i])
    end
    return result
end

"""
    faddeev_leverrier(A) -> FLResult

For `p(lambda) = lambda^n + c1 lambda^(n-1) + ... + cn`, iterate
`B0 = I`, `ck = -tr(A B(k-1))/k`, `Bk = A B(k-1) + ck I`.
Then `det(A) = (-1)^n cn` and `adj(A) = (-1)^(n+1) B(n-1)`.
Only exact integer divisions occur inside the recurrence.
"""
function faddeev_leverrier(a::Matrix{CNum})
    n, m = size(a)
    n == m || throw(DimensionMismatch("A must be square"))
    n > 0 || throw(ArgumentError("A must be nonempty"))
    b = identity_cnum(n)
    before_last = b
    c = _CNUM_ZERO
    @inbounds for k in 1:n
        ab = multiply(a, b)
        c = _neg_cnum(trace_cnum(ab)) / k
        k == n && break
        add_diagonal!(ab, c)
        b = ab
        k == n - 1 && (before_last = b)
    end
    determinant = isodd(n) ? _neg_cnum(c) : c
    adjugate = isodd(n + 1) ? scale(before_last, _to_cnum(-1)) : before_last
    return FLResult(determinant, adjugate)
end

function solve_adjugate(a::Matrix{CNum}, forcing::Vector{CNum})
    result = faddeev_leverrier(a)
    _iszero_cnum(result.determinant) && throw(ArgumentError("structurally singular system"))
    numerator = multiply(result.adjugate, forcing)
    return CNum[x / result.determinant for x in numerator], result, numerator
end

function solve_2x2(a::Matrix{CNum}, forcing::Vector{CNum})
    size(a) == (2, 2) || throw(DimensionMismatch())
    length(forcing) == 2 || throw(DimensionMismatch())
    determinant = _add_cnum(
        _mul_cnum(a[1, 1], a[2, 2]),
        _neg_cnum(_mul_cnum(a[1, 2], a[2, 1])),
    )
    _iszero_cnum(determinant) && throw(ArgumentError("structurally singular system"))
    n1 = _add_cnum(
        _mul_cnum(a[2, 2], forcing[1]),
        _neg_cnum(_mul_cnum(a[1, 2], forcing[2])),
    )
    n2 = _add_cnum(
        _mul_cnum(a[1, 1], forcing[2]),
        _neg_cnum(_mul_cnum(a[2, 1], forcing[1])),
    )
    return CNum[n1 / determinant, n2 / determinant], determinant, CNum[n1, n2]
end

"""
Solve only a structurally diagonal singular system. A zero diagonal with zero forcing gets
the deterministic zero coordinate; a zero diagonal with nonzero forcing is inconsistent.
This deliberately does not claim to solve connected rank-deficient blocks.
"""
function solve_structural_diagonal(a::Matrix{CNum}, forcing::Vector{CNum})
    n, m = size(a)
    n == m && m == length(forcing) || throw(DimensionMismatch())
    @inbounds for j in 1:n, i in 1:n
        i == j || _iszero_cnum(a[i, j]) || throw(
            ArgumentError("connected structurally singular blocks are unsupported"),
        )
    end
    result = fill(_CNUM_ZERO, n)
    @inbounds for i in 1:n
        if _iszero_cnum(a[i, i])
            _iszero_cnum(forcing[i]) || throw(
                ArgumentError("inconsistent forcing on structural null coordinate $i"),
            )
        else
            result[i] = forcing[i] / a[i, i]
        end
    end
    return result
end

"""Verify the division-free identity `A * numerator == determinant * forcing`."""
function residual_certificate(
        a::Matrix{CNum}, numerator::Vector{CNum}, determinant::CNum,
        forcing::Vector{CNum},
    )
    lhs = multiply(a, numerator)
    residual = CNum[
        _add_cnum(lhs[i], _neg_cnum(_mul_cnum(determinant, forcing[i])))
        for i in eachindex(lhs)
    ]
    return residual, all(_iszero_cnum, residual)
end

coefficient_size(c::CNum) = ncodeunits(string(to_num(c)))
matrix_size(a::AbstractArray{CNum}) = sum(coefficient_size, a; init = 0)

function _nodes(x)
    x isa Complex && return _nodes(real(x)) + _nodes(imag(x))
    x isa Num && return _nodes(Symbolics.unwrap(x))
    x isa Number && return 1
    SymbolicUtils.iscall(x) || return 1
    return 1 + sum(_nodes, SymbolicUtils.arguments(x); init = 0)
end

coefficient_nodes(c::CNum) = _nodes(to_num(c))

function coefficient_terms(c::CNum)
    tail = getfield(c, :tail)
    hasfield(typeof(tail), :terms) && return length(getfield(tail, :terms))
    return 1
end

end

