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

struct AffineDrift
    linear::CNum
    mixing::CNum
    offset::CNum
end

function affine_drift(heisenberg::SuperOperator, u::Op, v::Op)
    drift = heisenberg(u)
    isempty(drift.indices) || unitary_error(
        "`DisplacementFrame` requires a drift without bound sums; got `$(sprint(show, drift))`",
    )
    linear = CNUM_ZERO
    mixing = CNUM_ZERO
    offset = CNUM_ZERO
    for (term, c) in drift
        if isempty(term.ops)
            offset = add_cnum(offset, c)
        elseif term.ops == Op[u]
            linear = add_cnum(linear, c)
        elseif term.ops == Op[v]
            mixing = add_cnum(mixing, c)
        else
            unitary_error(
                "`DisplacementFrame` requires a drift affine in `$u` and `$v`; " *
                    "got `$(sprint(show, drift))`",
            )
        end
    end
    return AffineDrift(linear, mixing, offset)
end

struct LinearDriftPair
    m11::CNum
    m12::CNum
    m21::CNum
    m22::CNum
    offset1::CNum
    offset2::CNum
end

function fock_drift_pair(heisenberg::SuperOperator, d::Op, row::AffineDrift)
    raising = adjoint(d)
    adjoint_row = affine_drift(heisenberg, raising, d)
    offset = exponential_form(row.offset)
    adjoint_offset = exponential_form(adjoint_row.offset)
    (
        isequal(adjoint_row.linear, conj_cnum(row.linear)) &&
            isequal(adjoint_row.mixing, conj_cnum(row.mixing)) &&
            isequal(adjoint_offset, conj_cnum(offset))
    ) || unitary_error(
        "`DisplacementFrame` requires a Hermiticity-preserving drift, " *
            "`transpose(L)($raising) == transpose(L)($d)'`",
    )
    return LinearDriftPair(
        row.linear, row.mixing, adjoint_row.mixing, adjoint_row.linear,
        offset, adjoint_offset,
    )
end

function quadrature_drift_pair(heisenberg::SuperOperator, x::Op, p::Op)
    row_x = affine_drift(heisenberg, x, p)
    row_p = affine_drift(heisenberg, p, x)
    pair = LinearDriftPair(
        row_x.linear, row_x.mixing, row_p.mixing, row_p.linear,
        exponential_form(row_x.offset), exponential_form(row_p.offset),
    )
    for c in (pair.m11, pair.m12, pair.m21, pair.m22, pair.offset1, pair.offset2)
        isequal(c, conj_cnum(c)) || unitary_error(
            "`DisplacementFrame` requires a real quadrature drift; got coefficient " *
                "`$(to_num(c))`",
        )
    end
    return pair
end

function validate_drift_time(pair::LinearDriftPair, t::Num)
    variable = SymbolicUtils.unwrap(t)
    for c in (pair.m11, pair.m12, pair.m21, pair.m22)
        coefficient_depends_on(c, variable) && unitary_error(
            "`DisplacementFrame` requires a time-independent drift matrix; " *
                "got entry `$(to_num(c))`",
        )
    end
    return pair
end

@noinline function drift_resonance_error(drive1::CNum, drive2::CNum, harmonic::CNum)
    unitary_error(
        "`DisplacementFrame` found a resonant drive component " *
            "`($(to_num(drive1)), $(to_num(drive2)))` at frequency `$(to_num(harmonic))`",
    )
end

function drift_response(
        pair::LinearDriftPair, drive1::CNum, drive2::CNum, harmonic::CNum,
    )
    s = mul_cnum(CNUM_IM, harmonic)
    a11 = add_cnum(s, neg_cnum(pair.m11))
    a22 = add_cnum(s, neg_cnum(pair.m22))
    determinant = add_cnum(mul_cnum(a11, a22), neg_cnum(mul_cnum(pair.m12, pair.m21)))
    iszero_cnum(determinant) && drift_resonance_error(drive1, drive2, harmonic)
    first = add_cnum(mul_cnum(a22, drive1), mul_cnum(pair.m12, drive2)) / determinant
    second = add_cnum(mul_cnum(pair.m21, drive1), mul_cnum(a11, drive2)) / determinant
    return (first = first, second = second)
end

function static_drift_displacement(pair::LinearDriftPair)
    (iszero_cnum(pair.offset1) && iszero_cnum(pair.offset2)) &&
        return (first = CNUM_ZERO, second = CNUM_ZERO)
    return drift_response(pair, pair.offset1, pair.offset2, CNUM_ZERO)
end

struct HarmonicDrive
    harmonic::CNum
    drive1::CNum
    drive2::CNum
end

function add_harmonic!(
        groups::Vector{HarmonicDrive}, harmonic::CNum, component::CNum, slot::Int,
    )
    for i in eachindex(groups)
        group = groups[i]
        isequal(group.harmonic, harmonic) || continue
        groups[i] = slot == 1 ?
            HarmonicDrive(harmonic, add_cnum(group.drive1, component), group.drive2) :
            HarmonicDrive(harmonic, group.drive1, add_cnum(group.drive2, component))
        return groups
    end
    push!(
        groups, slot == 1 ? HarmonicDrive(harmonic, component, CNUM_ZERO) :
            HarmonicDrive(harmonic, CNUM_ZERO, component),
    )
    return groups
end

function collect_harmonics!(
        groups::Vector{HarmonicDrive}, drive::CNum, t::Num, slot::Int,
    )
    iszero_cnum(drive) && return groups
    tail = drive.tail
    if tail isa Native
        return add_harmonic!(groups, CNUM_ZERO, drive, slot)
    elseif !(tail isa Poly)
        unitary_error(
            "`DisplacementFrame` supports only finite harmonic drive coefficients; " *
                "got `$(to_num(drive))`",
        )
    end
    for monomial in tail.terms
        add_harmonic!(groups, harmonic_frequency(monomial, t), from_poly([monomial]), slot)
    end
    return groups
end

function bounded_drift_displacement(pair::LinearDriftPair, t::Num)
    groups = HarmonicDrive[]
    collect_harmonics!(groups, pair.offset1, t, 1)
    collect_harmonics!(groups, pair.offset2, t, 2)
    first = CNUM_ZERO
    second = CNUM_ZERO
    for group in groups
        response = drift_response(pair, group.drive1, group.drive2, group.harmonic)
        first = add_cnum(first, response.first)
        second = add_cnum(second, response.second)
    end
    return (first = first, second = second)
end

function drift_derivative(pair::LinearDriftPair, first::CNum, second::CNum)
    return (
        first = add_cnum(
            add_cnum(mul_cnum(pair.m11, first), mul_cnum(pair.m12, second)), pair.offset1,
        ),
        second = add_cnum(
            add_cnum(mul_cnum(pair.m21, first), mul_cnum(pair.m22, second)), pair.offset2,
        ),
    )
end

function moving_fock_displacement(d::Op, α::CNum, derivative::CNum, t::Num)
    derivative_adjoint = conj_cnum(derivative)
    gauge = fock_displacement_gauge(
        d, α, derivative, derivative_adjoint,
        mul_cnum(CNUM_NEG_IM, derivative), mul_cnum(CNUM_IM, derivative_adjoint),
    )
    return timed_transform(fock_displacement(d, α), gauge, t)
end

function scalar_drift_rate(row::AffineDrift)
    iszero_cnum(row.linear) && unitary_error(
        "`DisplacementFrame` found a drift with no linear response",
    )
    return row.linear
end

"""
    DisplacementFrame(a, L::Union{Lindbladian, SuperOperator})

Construct the static displacement `a ↦ a + α` that removes the coherent drive from the
linear generator `L`, a [`Lindbladian`](@ref) or any [`SuperOperator`](@ref) such as a
Caldeira–Leggett master equation. The Heisenberg drift must be affine in the mode,
`d/dt (a, a') = M (a, a') + c`, and `(α, conj(α)) = -M⁻¹ c` is its fixed point, so heating
and loss enter through the exact drift rather than through a no-jump Hamiltonian. The drift
must preserve Hermiticity, including when it does not mix `a` with `a'`:
`transpose(L)(a') == transpose(L)(a)'`. A structurally singular `M` is rejected.
"""
function DisplacementFrame(a::Op, L::Union{Lindbladian{QAdd}, SuperOperator})
    d = fock_or_throw(a, "`DisplacementFrame`")
    heisenberg = transpose(L)
    row = affine_drift(heisenberg, d, adjoint(d))
    pair = fock_drift_pair(heisenberg, d, row)
    if iszero_cnum(row.mixing)
        α = neg_cnum(row.offset) / scalar_drift_rate(row)
        return fock_displacement(d, α)
    end
    return fock_displacement(d, static_drift_displacement(pair).first)
end

"""
    DisplacementFrame(a, L::Union{Lindbladian, SuperOperator}, t)

Construct the bounded moving displacement of the linear generator `L`, whose drift
`d/dt (a, a') = M (a, a') + c(t)` has a time-independent matrix `M` and a drive `c(t)` that
is a finite sum of exact harmonic phases. The displacement solves `α̇ = M₁₁ α + M₁₂ conj(α)
+ c₁(t)` with the free homogeneous solution set to zero: each harmonic `expim(ν*t)` is
answered by the exact solution of `(im*ν - M) v = c_ν`.

A structurally resonant harmonic is rejected. Symbolic determinants are retained as exact
quotients and therefore describe the response away from their resonance surfaces.
"""
function DisplacementFrame(a::Op, L::Union{Lindbladian{QAdd}, SuperOperator}, t::Num)
    d = fock_or_throw(a, "`DisplacementFrame`")
    tt = time_or_throw(t)
    heisenberg = transpose(L)
    row = affine_drift(heisenberg, d, adjoint(d))
    pair = fock_drift_pair(heisenberg, d, row)
    if iszero_cnum(row.mixing)
        rate = scalar_drift_rate(row)
        coefficient_depends_on(rate, SymbolicUtils.unwrap(tt)) && unitary_error(
            "`DisplacementFrame` requires a time-independent drift rate; " *
                "got `$(to_num(rate))`",
        )
        offset = exponential_form(row.offset)
        amplitude = bounded_harmonic_displacement(
            mul_cnum(CNUM_IM, offset), mul_cnum(CNUM_IM, rate), tt,
        )
        derivative = add_cnum(mul_cnum(rate, amplitude), offset)
        return moving_fock_displacement(d, amplitude, derivative, tt)
    end
    pair = validate_drift_time(pair, tt)
    α = bounded_drift_displacement(pair, tt).first
    derivative = drift_derivative(pair, α, conj_cnum(α)).first
    return moving_fock_displacement(d, α, derivative, tt)
end

"""
    DisplacementFrame(x, p, L::Union{Lindbladian, SuperOperator})

Construct the static translation `x ↦ x + X`, `p ↦ p + P` that removes the drive from the
linear generator `L` of a canonical `(Position, Momentum)` pair. The real Heisenberg drift
`d/dt (x, p) = M (x, p) + c` must be affine in the pair, and `(X, P) = -M⁻¹ c` is its
fixed point. A structurally singular driven response is rejected.
"""
function DisplacementFrame(x::Op, p::Op, L::Union{Lindbladian{QAdd}, SuperOperator})
    phase_pair(x, p, "`DisplacementFrame`")
    pair = quadrature_drift_pair(transpose(L), x, p)
    shift = static_drift_displacement(pair)
    return quadrature_displacement(x, p, shift.first, shift.second)
end

"""
    DisplacementFrame(x, p, L::Union{Lindbladian, SuperOperator}, t)

Construct the bounded moving translation of a canonical `(Position, Momentum)` pair driven
by the linear generator `L`, whose real drift `d/dt (x, p) = M (x, p) + c(t)` has a
time-independent `M` and finite harmonic drives. The homogeneous solution is set to zero,
the moving-frame gauge is the exact Weyl gauge of [`Displace`](@ref), and structural
resonances are rejected.
"""
function DisplacementFrame(
        x::Op, p::Op, L::Union{Lindbladian{QAdd}, SuperOperator}, t::Num,
    )
    phase_pair(x, p, "`DisplacementFrame`")
    tt = time_or_throw(t)
    pair = validate_drift_time(quadrature_drift_pair(transpose(L), x, p), tt)
    shift = bounded_drift_displacement(pair, tt)
    derivative = drift_derivative(pair, shift.first, shift.second)
    gauge = quadrature_displacement_gauge(
        x, p, shift.first, shift.second, derivative.first, derivative.second,
        derivative.second, neg_cnum(derivative.first),
    )
    static = quadrature_displacement(x, p, shift.first, shift.second)
    return timed_transform(static, gauge, tt)
end
