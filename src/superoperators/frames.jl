function reduce_params(S::SuperOperator, relations::Vector{ParamRelation})
    isempty(relations) && return S
    scratch = ParamRelation[]
    return map_coefficients(c -> reduce_all(c, relations, true, scratch), S)
end

function conjugate(S::SuperOperator, U::UnitaryTransform)
    conjugated = map_sides((ops, ne) -> conjugate(side_qadd(ops, ne), U), S)
    return reduce_params(conjugated, U.action.relations)
end

transform(S::SuperOperator, U::UnitaryTransform{StaticTime}) = conjugate(S, U)
transform(S::SuperOperator, U::UnitaryTransform{DynamicTime}) =
    conjugate(S, U) + SuperOperator(hamiltonian_action(U.gauge))

transform(L::Lindbladian{O}, U::UnitaryTransform) where {O} =
    map_lindbladian(H -> transform(H, U), J -> conjugate(J, U), identity, same_scope, O, L)

"""
    local_lindbladian(L, space)

The Hamiltonian terms and channels of `L` that act on `space`. Every other part commutes with
an operator on `space`, so its Heisenberg action on that operator vanishes.
"""
function local_lindbladian(L::Lindbladian{QAdd}, space::Int32)
    H = QTermDict()
    for (t, c) in L.hamiltonian.arguments
        acts_on_space(t.ops, space) && (H[copy_key(t)] = c)
    end
    chans = Dict{Channel{QAdd}, CNum}()
    for (ch, γ) in L.channels
        (acts_on_space(ch.jump, space) || acts_on_space(ch.partner, space)) && (chans[ch] = γ)
    end
    return Lindbladian{QAdd}(QAdd(H, drop_unused_indices(H, copy(L.hamiltonian.indices))), chans)
end

acts_on_space(ops::Vector{Op}, space::Int32) = any(o -> o.space_index == space, ops)
acts_on_space(q::QAdd, space::Int32) = any(t -> acts_on_space(t.ops, space), keys(q.arguments))

function lindblad_drift(basis::Vector{Op}, primary::Int, L::Lindbladian{QAdd})
    n = length(basis)
    locations = Dict{Op, Int}(operator => i for (i, operator) in enumerate(basis))
    rate = fill(CNUM_ZERO, n, n)
    drive = fill(CNUM_ZERO, n)
    heisenberg = Dict{Int32, SuperOperator}()
    for (row, operator) in enumerate(basis)
        space = operator.space_index
        drift = get!(() -> transpose(local_lindbladian(L, space)), heisenberg, space)(operator)
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
