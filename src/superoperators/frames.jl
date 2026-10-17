function reduce_params(S::SuperOperator, relations::Vector{ParamRelation})
    isempty(relations) && return S
    scratch = ParamRelation[]
    return map_coefficients(c -> reduce_all(c, relations, true, scratch), S)
end

function conjugate(S::SuperOperator, U::UnitaryTransform)
    conjugated = map_sides((ops, ne) -> conjugate(side_qadd(ops, ne), U), S)
    return reduce_params(conjugated, U.action.relations)
end

conjugate(L::Lindbladian, U::UnitaryTransform) =
    map_lindbladian(H -> conjugate(H, U), J -> conjugate(J, U), identity, L)

transform(S::SuperOperator, U::UnitaryTransform{StaticTime}) = conjugate(S, U)
transform(S::SuperOperator, U::UnitaryTransform{DynamicTime}) =
    conjugate(S, U) + SuperOperator(hamiltonian_action(U.gauge))

transform(L::Lindbladian, U::UnitaryTransform) =
    map_lindbladian(H -> transform(H, U), J -> conjugate(J, U), identity, L)

function lindblad_drift(basis::Vector{Op}, primary::Int, L::Lindbladian{QAdd})
    n = length(basis)
    locations = Dict{Op, Int}(operator => i for (i, operator) in enumerate(basis))
    rate = fill(CNUM_ZERO, n, n)
    drive = fill(CNUM_ZERO, n)
    heisenberg = transpose(L)
    for (row, operator) in enumerate(basis)
        drift = heisenberg(operator)
        isempty(drift.indices) || unitary_error(
            "`DisplacementFrame` requires a drift without bound sums; got `$(sprint(show, drift))`",
        )
        for (term, c) in drift
            add_adjoint_entry!(rate, drive, row, locations, term.ops, c) || unitary_error(
                "`DisplacementFrame` requires a drift affine in the displaced modes; " *
                    "got `$(sprint(show, drift))`",
            )
        end
    end
    return LinearDrift(basis, primary, rate, drive, CNUM_IM)
end

DisplacementFrame(modes::FockModes, L::Lindbladian{QAdd}) =
    displacement_frame(lindblad_drift(fock_frame_basis(modes)..., L))
DisplacementFrame(modes::FockModes, L::Lindbladian{QAdd}, t::Num) =
    displacement_frame(lindblad_drift(fock_frame_basis(modes)..., L), t)
