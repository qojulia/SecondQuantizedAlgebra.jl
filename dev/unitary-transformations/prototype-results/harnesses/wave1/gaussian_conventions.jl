module GaussianConventionsPrototype

using SecondQuantizedAlgebra
using Symbolics

import SecondQuantizedAlgebra: Coeff, Op, QAdd, UnitaryTransform,
    _CNUM_IM, _CNUM_NEG_IM, _CNUM_ONE, _CNUM_ZERO,
    _add_cnum, _conj_cnum, _dt, _iszero_cnum, _mul_cnum, _neg_cnum,
    _reduce_params, _rule_qadd, _substitute_cnum, _to_cnum, _to_complex,
    _zero_qadd, to_num

export affine_coordinates, commutator_matrix, direct_affine_action,
    interleaved_order, split_order, scanned_affine_action,
    reconstructed_gauge, scalar_coefficient, operator_part,
    velocity_residual, materialize_rows

interleaved_order(lowering::Vector{Op}) =
    reduce(vcat, (Op[a, adjoint(a)] for a in lowering); init = Op[])

split_order(lowering::Vector{Op}) = Op[lowering; adjoint.(lowering)]

function commutator_matrix(basis::Vector{Op})
    n = length(basis)
    result = fill(_CNUM_ZERO, n, n)
    for j in 1:n, k in 1:n
        q = commutator(basis[j], basis[k])
        for (term, coefficient) in q
            isempty(term.ops) || error("basis is not canonical: non-scalar commutator")
            result[j, k] = _add_cnum(result[j, k], coefficient)
        end
    end
    return result
end

function adjoint_permutation(basis::Vector{Op})
    [something(findfirst(isequal(adjoint(op)), basis), 0) for op in basis]
end

function affine_coordinates(q::QAdd, basis::Vector{Op})
    linear = fill(_CNUM_ZERO, length(basis))
    scalar = _CNUM_ZERO
    for (term, coefficient) in q
        if isempty(term.ops)
            scalar = _add_cnum(scalar, coefficient)
        elseif length(term.ops) == 1
            index = findfirst(isequal(first(term.ops)), basis)
            index === nothing && error("operator outside affine basis: $(first(term.ops))")
            linear[index] = _add_cnum(linear[index], coefficient)
        else
            error("non-affine term in action: $(term.ops)")
        end
    end
    return linear, scalar
end

affine_coordinates(op::Op, basis::Vector{Op}) =
    affine_coordinates(1 * op, basis)

function direct_affine_action(H::QAdd, basis::Vector{Op})
    n = length(basis)
    action = fill(_CNUM_ZERO, n, n)
    forcing = fill(_CNUM_ZERO, n)
    for row in 1:n
        action[row, :], forcing[row] = affine_coordinates(
            im * commutator(H, basis[row]), basis,
        )
    end
    return action, forcing
end

"""Scan a degree-two Hamiltonian once using [xy,z]=x[y,z]+[x,z]y."""
function scanned_affine_action(H::QAdd, basis::Vector{Op})
    n = length(basis)
    C = commutator_matrix(basis)
    predicted = [_zero_qadd() for _ in 1:n]
    for (term, coefficient) in H
        degree = length(term.ops)
        degree == 0 && continue
        degree <= 2 || error("Hamiltonian is not quadratic")
        indices = map(term.ops) do op
            index = findfirst(isequal(op), basis)
            index === nothing && error("operator outside canonical basis: $op")
            index
        end
        if degree == 1
            i = only(indices)
            for k in 1:n
                factor = _mul_cnum(_CNUM_IM, _mul_cnum(coefficient, C[i, k]))
                _iszero_cnum(factor) ||
                    (predicted[k] = predicted[k] + to_num(factor))
            end
        else
            i, j = indices
            for k in 1:n
                left = _mul_cnum(_CNUM_IM, _mul_cnum(coefficient, C[j, k]))
                right = _mul_cnum(_CNUM_IM, _mul_cnum(coefficient, C[i, k]))
                _iszero_cnum(left) ||
                    (predicted[k] = predicted[k] + to_num(left) * basis[i])
                _iszero_cnum(right) ||
                    (predicted[k] = predicted[k] + to_num(right) * basis[j])
            end
        end
    end
    action = fill(_CNUM_ZERO, n, n)
    forcing = fill(_CNUM_ZERO, n)
    for row in 1:n
        action[row, :], forcing[row] = affine_coordinates(predicted[row], basis)
    end
    return action, forcing
end

function rule_affine_coordinates(U::UnitaryTransform, basis::Vector{Op}; inverse = false)
    rules = inverse ? U.inverse_rules : U.rules
    n = length(basis)
    matrix = fill(_CNUM_ZERO, n, n)
    shift = fill(_CNUM_ZERO, n)
    for row in 1:n
        matrix[row, :], shift[row] = affine_coordinates(rules[basis[row]], basis)
    end
    return matrix, shift
end

function _matmul(A::Matrix{Coeff}, B::Matrix{Coeff})
    size(A, 2) == size(B, 1) || throw(DimensionMismatch())
    C = fill(_CNUM_ZERO, size(A, 1), size(B, 2))
    for i in axes(A, 1), j in axes(B, 2), k in axes(A, 2)
        C[i, j] = _add_cnum(C[i, j], _mul_cnum(A[i, k], B[k, j]))
    end
    return C
end

function _matvec(A::Matrix{Coeff}, x::Vector{Coeff})
    length(x) == size(A, 2) || throw(DimensionMismatch())
    y = fill(_CNUM_ZERO, size(A, 1))
    for i in axes(A, 1), k in axes(A, 2)
        y[i] = _add_cnum(y[i], _mul_cnum(A[i, k], x[k]))
    end
    return y
end

_differentiate(A::Matrix{Coeff}, t::Num) = map(c -> _dt(c, t), A)
_differentiate(x::Vector{Coeff}, t::Num) = map(c -> _dt(c, t), x)

function _inverse_iC(C::Matrix{Coeff})
    # Both candidate bosonic commutator matrices are monomial matrices.  Invert
    # iC exactly without relying on `zero(::Coeff)` in generic LinearAlgebra.
    n = size(C, 1)
    R = fill(_CNUM_ZERO, n, n)
    for row in 1:n
        column = findfirst(j -> !_iszero_cnum(C[row, j]), 1:n)
        column === nothing && error("singular commutator matrix")
        R[column, row] = inv(_mul_cnum(_CNUM_IM, C[row, column]))
    end
    return R
end

"""
Reconstruct the non-central gauge from z' = Mz+d and its derivative.

With G = 1/2 z'Kz + l'z + s and package convention
G = i*(dU'/dt)*U, dot(z') = -i[G,z'].  The stored gauge is expressed in
the untransformed generator labels, so M^-1*dot(M)=i*C*K and
M^-1*dot(d)=i*C*l.  `lift_scalar` is intentionally explicit because the
affine action cannot determine it.
"""
function reconstructed_gauge(
        U::UnitaryTransform, basis::Vector{Op}, t::Num;
        lift_scalar = _CNUM_ZERO,
    )
    M, d = rule_affine_coordinates(U, basis)
    Minv, _ = rule_affine_coordinates(U, basis; inverse = true)
    body_action = _matmul(Minv, _differentiate(M, t))
    body_forcing = _matvec(Minv, _differentiate(d, t))
    inv_iC = _inverse_iC(commutator_matrix(basis))
    K = _matmul(inv_iC, body_action)
    l = _matvec(inv_iC, body_forcing)

    gauge = _zero_qadd()
    for i in eachindex(basis), j in eachindex(basis)
        coefficient = _mul_cnum(_to_cnum(1 // 2), K[i, j])
        _iszero_cnum(coefficient) ||
            (gauge = gauge + to_num(coefficient) * basis[i] * basis[j])
    end
    for i in eachindex(basis)
        _iszero_cnum(l[i]) || (gauge = gauge + to_num(l[i]) * basis[i])
    end
    _iszero_cnum(lift_scalar) || (gauge = gauge + to_num(lift_scalar))
    return _reduce_params(gauge, U.relations, true)
end

function scalar_coefficient(q::QAdd)
    scalar = _CNUM_ZERO
    for (term, coefficient) in q
        isempty(term.ops) && (scalar = _add_cnum(scalar, coefficient))
    end
    return scalar
end

function operator_part(q::QAdd)
    terms = Tuple{Coeff, Vector{Op}}[]
    for (term, coefficient) in q
        isempty(term.ops) || push!(terms, (coefficient, copy(term.ops)))
    end
    return _rule_qadd(terms)
end

function velocity_residual(U::UnitaryTransform, basis::Vector{Op}, t::Num)
    residuals = QAdd[]
    for z in basis
        image = conjugate(z, U)
        derivative = SecondQuantizedAlgebra._map_coefficients(c -> _dt(c, t), image)
        residual = simplify(derivative + im * commutator(U.gauge, image))
        push!(residuals, _reduce_params(residual, U.relations, true))
    end
    return residuals
end

"""Minimal ordering-neutral coefficient-row to rule materializer."""
function materialize_rows(matrix::Matrix{Coeff}, basis::Vector{Op})
    size(matrix) == (length(basis), length(basis)) || throw(DimensionMismatch())
    rules = Dict{Op, QAdd}()
    for row in eachindex(basis)
        pairs = Tuple{Coeff, Vector{Op}}[]
        for column in eachindex(basis)
            coefficient = matrix[row, column]
            _iszero_cnum(coefficient) ||
                push!(pairs, (coefficient, Op[basis[column]]))
        end
        rules[basis[row]] = _rule_qadd(pairs)
    end
    return rules
end

end

