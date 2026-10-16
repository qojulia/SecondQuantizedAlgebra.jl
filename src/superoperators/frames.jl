function reduce_map_params(S::SuperOperator, relations::Vector{ParamRelation})
    isempty(relations) && return S
    scratch = ParamRelation[]
    return map_coefficients(c -> reduce_all(c, relations, true, scratch), S)
end

function conjugate(S::SuperOperator, U::UnitaryTransform)
    conjugated = map_sides((ops, ne) -> conjugate(side_qadd(ops, ne), U), S)
    return reduce_map_params(conjugated, U.action.relations)
end

conjugate(L::Lindbladian, U::UnitaryTransform) =
    map_lindbladian(H -> conjugate(H, U), J -> conjugate(J, U), identity, L)

transform(S::SuperOperator, U::UnitaryTransform{StaticTime}) = conjugate(S, U)
transform(S::SuperOperator, U::UnitaryTransform{DynamicTime}) =
    conjugate(S, U) + SuperOperator(hamiltonian_action(U.gauge))

transform(L::Lindbladian, U::UnitaryTransform) =
    map_lindbladian(H -> transform(H, U), J -> conjugate(J, U), identity, L)

struct FockDrift
    rate::CNum
    offset::CNum
end

function linear_fock_drift(L::Lindbladian{QAdd}, d::Op)
    drift = transpose(L)(d)
    isempty(drift.indices) || unitary_error(
        "`DisplacementFrame` requires a drift without bound sums; got `$(sprint(show, drift))`",
    )
    rate = CNUM_ZERO
    offset = CNUM_ZERO
    for (term, c) in drift
        if isempty(term.ops)
            offset = offset + c
        elseif term.ops == Op[d]
            rate = rate + c
        else
            unitary_error(
                "`DisplacementFrame` requires a drift affine in `$d`; " *
                    "got `$(sprint(show, drift))`",
            )
        end
    end
    iszero(rate) && unitary_error(
        "`DisplacementFrame` found a drift with no linear response in `$d`",
    )
    return FockDrift(rate, offset)
end

function DisplacementFrame(a::Op, L::Lindbladian{QAdd})
    d = fock_or_throw(a, "`DisplacementFrame`")
    drift = linear_fock_drift(L, d)
    α = -drift.offset / drift.rate
    return fock_displacement(d, α)
end

function DisplacementFrame(a::Op, L::Lindbladian{QAdd}, t::Num)
    d = fock_or_throw(a, "`DisplacementFrame`")
    tt = time_or_throw(t)
    drift = linear_fock_drift(L, d)
    coefficient_depends_on(drift.rate, SymbolicUtils.unwrap(tt)) && unitary_error(
        "`DisplacementFrame` requires a time-independent drift rate; " *
            "got `$(to_num(drift.rate))`",
    )
    offset = exponential_form(drift.offset)
    amplitude = bounded_harmonic_displacement(
        CNUM_IM * offset, CNUM_IM * drift.rate, tt,
    )
    derivative = drift.rate * amplitude + offset
    derivative_adjoint = conj(derivative)
    gauge = fock_displacement_gauge(
        d, amplitude, derivative, derivative_adjoint,
        CNUM_NEG_IM * derivative, CNUM_IM * derivative_adjoint,
    )
    return timed_transform(fock_displacement(d, amplitude), gauge, tt)
end
