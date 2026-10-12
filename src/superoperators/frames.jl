function reduce_map_params(S::SuperOperator, relations::Vector{ParamRelation})
    isempty(relations) && return S
    scratch = ParamRelation[]
    return map_coefficients(c -> reduce_all(c, relations, true, scratch), S)
end

"""
    conjugate(S, U::UnitaryTransform)

Conjugate every side of a map, `A ↦ U' * A * U`. For a [`Lindbladian`](@ref) the
Hamiltonian and every jump are conjugated, so `𝒟[J] ↦ 𝒟[U' J U]`. Conjugation is a
homomorphism: `conjugate(A ∘ B, U) == conjugate(A, U) ∘ conjugate(B, U)`.
"""
function conjugate(S::SuperOperator, U::UnitaryTransform)
    conjugated = map_sides((ops, ne) -> conjugate(side_qadd(ops, ne), U), S)
    return reduce_map_params(conjugated, U.action.relations)
end

conjugate(L::Lindbladian, U::UnitaryTransform) =
    map_lindbladian(H -> conjugate(H, U), J -> conjugate(J, U), identity, L)

"""
    transform(L, U::UnitaryTransform)

Move a Schrödinger-picture generator into the frame of `U`: `conjugate(L, U)` plus the
coherent gauge term `hamiltonian_action(gauge_term(U))`. For a [`Lindbladian`](@ref) the
Hamiltonian is transformed with [`transform`](@ref) and the jumps are conjugated; jump
constants created by a displacement are folded back into the Hamiltonian. `transform` is
affine, so it does not distribute over composition.
"""
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
            offset = add_cnum(offset, c)
        elseif term.ops == Op[d]
            rate = add_cnum(rate, c)
        else
            unitary_error(
                "`DisplacementFrame` requires a drift affine in `$d`; " *
                    "got `$(sprint(show, drift))`",
            )
        end
    end
    iszero_cnum(rate) && unitary_error(
        "`DisplacementFrame` found a drift with no linear response in `$d`",
    )
    return FockDrift(rate, offset)
end

"""
    DisplacementFrame(a, Llin::Lindbladian)

Construct the static displacement `a ↦ a + α` that removes the coherent drive from the
damped linear reference generator `Llin`. The drift `transpose(Llin)(a) = λ a + c` must be
affine in `a`, and `α = -c/λ` is its fixed point, so heating and loss enter through the
exact drift rather than through a no-jump Hamiltonian.
"""
function DisplacementFrame(a::Op, L::Lindbladian{QAdd})
    d = fock_or_throw(a, "`DisplacementFrame`")
    drift = linear_fock_drift(L, d)
    α = neg_cnum(drift.offset) / drift.rate
    return fock_displacement(d, α)
end

"""
    DisplacementFrame(a, Llin::Lindbladian, t)

Construct the bounded moving displacement of the damped linear generator `Llin`, whose
drift `transpose(Llin)(a) = λ a + c(t)` has a time-independent `λ` and a drive `c(t)` that
is a finite sum of exact harmonic phases. The displacement solves `α̇ = λ α + c(t)` with the
free homogeneous solution set to zero.
"""
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
        mul_cnum(CNUM_IM, offset), mul_cnum(CNUM_IM, drift.rate), tt,
    )
    derivative = add_cnum(mul_cnum(drift.rate, amplitude), offset)
    derivative_adjoint = conj_cnum(derivative)
    gauge = fock_displacement_gauge(
        d, amplitude, derivative, derivative_adjoint,
        mul_cnum(CNUM_NEG_IM, derivative), mul_cnum(CNUM_IM, derivative_adjoint),
    )
    return timed_transform(fock_displacement(d, amplitude), gauge, tt)
end
