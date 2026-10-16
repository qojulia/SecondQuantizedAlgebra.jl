function hermitian_generator_or_throw(G::QAdd)
    iszero(simplify(G - adjoint(G))) || unitary_error(
        "`UnitaryTransform(G, θ)` requires `G` to be provably Hermitian exactly",
    )
    return G
end

function phase_space_site_basis(operators::Vector{Op})
    position_index = position_of(is_position, operators)
    momentum_index = position_of(is_momentum, operators)
    iszero(position_index) && unitary_error(
        "`UnitaryTransform(G, θ)` cannot infer a missing position rule for a phase-space " *
            "site; include both canonical quadratures in the generator",
    )
    iszero(momentum_index) && unitary_error(
        "`UnitaryTransform(G, θ)` cannot infer a missing momentum rule for a phase-space " *
            "site; include both canonical quadratures in the generator",
    )
    return Op[operators[position_index], operators[momentum_index]]
end

function closed_adjoint_basis(G::QAdd)
    isempty(G.indices) || unitary_error(
        "`UnitaryTransform(G, θ)` does not yet support summed generators",
    )

    by_site = Dict{SiteKey, Vector{Op}}()
    for (term, _) in G
        isempty(term.ne) || unitary_error(
            "`UnitaryTransform(G, θ)` does not support constrained generator terms",
        )
        for operator in term.ops
            is_free_family(operator) && unitary_error(
                "`UnitaryTransform(G, θ)` does not yet support free indexed families",
            )
            push!(get!(by_site, site_key(operator), Op[]), operator)
        end
    end
    isempty(by_site) && unitary_error(
        "a scalar generator has a trivial adjoint action and no transformed operator basis",
    )

    raw = Op[]
    for operators in values(by_site)
        first_operator = first(operators)
        local_basis = if is_phase_space(first_operator)
            phase_space_site_basis(operators)
        else
            site_generators(first_operator)
        end
        isempty(local_basis) && unitary_error(
            "no finite exact generator basis is registered for `$first_operator`",
        )
        append!(raw, local_basis)
    end

    structure = all(is_fock, raw) ? AFFINE_BOSONIC_NAMBU :
        all(is_phase_space, raw) ? AFFINE_SYMPLECTIC_PHASE_SPACE : AFFINE_UNITARY_LINEAR
    return canonical_block_basis(structure, raw)
end

function affine_commutator_data(G::QAdd, basis::Vector{Op})
    n = length(basis)
    locations = Dict{Op, Int}(operator => i for (i, operator) in enumerate(basis))
    linear = fill(CNUM_ZERO, n, n)
    shift = fill(CNUM_ZERO, n)

    for (row, generator) in enumerate(basis)
        image = simplify(commutator(G, generator))
        isempty(image.indices) || unitary_error(
            "the adjoint action on `$generator` produced a summed expression",
        )
        for (term, coefficient) in image
            isempty(term.ne) || unitary_error(
                "the adjoint action on `$generator` produced constrained terms",
            )
            if isempty(term.ops)
                shift[row] = shift[row] + coefficient
            elseif length(term.ops) == 1
                target = first(term.ops)
                column = get(locations, target, 0)
                column == 0 && unitary_error(
                    "the adjoint action leaves the finite generator basis: " *
                        "`[$G, $generator]` contains `$target`",
                )
                linear[row, column] = linear[row, column] + coefficient
            else
                unsupported = rule_qadd((coefficient, copy(term.ops)))
                unitary_error(
                    "the adjoint action on `$generator` is not affine-linear; got " *
                        "`$(sprint(show, unsupported))`",
                )
            end
        end
    end
    return linear, shift
end

function coefficient_equal(left::CNum, right::CNum)
    return iszero(left - right)
end

function exact_real_argument(rate::CNum, θ::Real)
    coefficient_equal(rate, CNUM_ONE) && return θ
    coefficient_equal(rate, CNUM_NEG1) && return -θ
    return real(to_num(rate)) * θ
end

function generated_scalar_exponential(flow::CNum, θ::Real)
    iszero(flow) && return CNUM_ONE
    if coefficient_is_real(flow)
        coefficient_equal(flow, CNUM_ONE) && return to_cnum(exp(θ))
        coefficient_equal(flow, CNUM_NEG1) && return to_cnum(inv(exp(θ)))
        return to_cnum(exp(exact_real_argument(flow, θ)))
    end

    phase_rate = CNUM_NEG_IM * flow
    if coefficient_is_real(phase_rate)
        return phase(exact_real_argument(phase_rate, θ))
    end

    unitary_error(
        "no exact scalar exponential strategy supports the adjoint rate `$(to_num(flow))`",
    )
end

function connected_components(linear::Matrix{CNum})
    n = size(linear, 1)
    visited = falses(n)
    components = Vector{Vector{Int}}()
    for start in 1:n
        visited[start] && continue
        stack = Int[start]
        component = Int[]
        visited[start] = true
        while !isempty(stack)
            current = pop!(stack)
            push!(component, current)
            for other in 1:n
                other == current && continue
                connected = !iszero(linear[current, other]) ||
                    !iszero(linear[other, current])
                connected || continue
                visited[other] && continue
                visited[other] = true
                push!(stack, other)
            end
        end
        sort!(component)
        push!(components, component)
    end
    return components
end

function exact_rotation_pair(rate::CNum, θ::Real)
    if coefficient_equal(rate, CNUM_ONE)
        return to_cnum(cos(θ)), to_cnum(sin(θ))
    elseif coefficient_equal(rate, CNUM_NEG1)
        return to_cnum(cos(θ)), -to_cnum(sin(θ))
    end
    angle = exact_real_argument(rate, θ)
    return to_cnum(cos(angle)), to_cnum(sin(angle))
end

function rational_square_root(square::CNum)
    none = (0, 0 // 1)
    is_native(square) || return none
    z = native_scalar(square)
    iszero(imag(z)) || return none
    x = real(z)
    (x isa Rational{Int} && !iszero(x)) || return none
    magnitude = Rational{BigInt}(abs(x))
    is_rational_power(magnitude, 1 // 2) || return none
    return Int(sign(x)), Rational{Int}(rational_power_value(magnitude, 1 // 2))
end

is_real_scaled_phase(terms::Vector{PhaseTerm}) =
    length(terms) == 1 && coefficient_is_real(only(terms).amplitude)

function trig_pair!(relations::Vector{ParamRelation}, angle::Real)
    push_relation!(relations, trig_rel(angle))
    return to_cnum(cos(angle)), to_cnum(sin(angle))
end

function hyperbolic_pair!(relations::Vector{ParamRelation}, angle::Real)
    push_relation!(relations, hyp_rel(angle))
    return to_cnum(cosh(angle)), to_cnum(sinh(angle))
end

function set_two_by_two!(result::Matrix{CNum}, i::Int, j::Int, diagonal::CNum, upper::CNum, lower::CNum)
    result[i, i] = diagonal
    result[i, j] = upper
    result[j, i] = lower
    result[j, j] = diagonal
    return result
end

function exact_two_by_two_flow!(
        result::Matrix{CNum}, flow::Matrix{CNum}, i::Int, j::Int, θ::Real,
        relations::Vector{ParamRelation},
    )
    iszero(flow[i, i]) && iszero(flow[j, j]) || unitary_error(
        "exact two-dimensional adjoint blocks currently require zero diagonal entries",
    )
    upper = flow[i, j]
    lower = flow[j, i]

    if coefficient_is_real(upper) && coefficient_is_real(lower) &&
            coefficient_equal(lower, -upper)
        cosine, sine = exact_rotation_pair(upper, θ)
        push_relation!(relations, trig_rel(exact_real_argument(upper, θ)))
        return set_two_by_two!(result, i, j, cosine, sine, -sine)
    end

    if coefficient_equal(lower, conj(upper))
        if coefficient_equal(upper * lower, CNUM_ONE)
            ch, sh = hyperbolic_pair!(relations, θ)
            return set_two_by_two!(result, i, j, ch, upper * sh, lower * sh)
        elseif coefficient_is_real(upper) && coefficient_equal(lower, upper)
            ch, sh = hyperbolic_pair!(relations, exact_real_argument(upper, θ))
            return set_two_by_two!(result, i, j, ch, sh, sh)
        end
        terms = phase_terms(upper)
        if is_real_scaled_phase(terms)
            direction = phase(only(terms).phase)
            ch, sh = hyperbolic_pair!(relations, exact_real_argument(only(terms).amplitude, θ))
            return set_two_by_two!(result, i, j, ch, direction * sh, conj(direction) * sh)
        end
    end

    sign_of_square, rate = rational_square_root(upper * lower)
    if !iszero(sign_of_square)
        angle = exact_real_argument(to_cnum(rate), θ)
        diagonal, off_diagonal = sign_of_square < 0 ? trig_pair!(relations, angle) :
            hyperbolic_pair!(relations, angle)
        scaled = to_cnum(inv(rate)) * off_diagonal
        return set_two_by_two!(result, i, j, diagonal, upper * scaled, lower * scaled)
    end

    unitary_error(
        "no exact two-dimensional exponential strategy supports the closed adjoint block " *
            "[$(to_num(upper)), $(to_num(lower))]",
    )
end

function exact_closed_adjoint_action(
        basis::Vector{Op}, linear::Matrix{CNum}, shift::Vector{CNum}, θ::Real,
    )
    n = length(basis)
    flow = Matrix{CNum}(undef, n, n)
    affine_flow = Vector{CNum}(undef, n)
    for j in 1:n, i in 1:n
        flow[i, j] = CNUM_IM * linear[i, j]
    end
    for i in 1:n
        affine_flow[i] = CNUM_IM * shift[i]
    end

    has_linear = any(value -> !iszero(value), flow)
    has_shift = any(value -> !iszero(value), affine_flow)
    if has_shift
        has_linear && unitary_error(
            "closed affine adjoint actions with simultaneous linear and translation parts " *
                "do not yet have an exact exponential strategy",
        )
        parameter = to_cnum(θ)
        offsets = CNum[parameter * value for value in affine_flow]
        return AffineAction(basis, identity_linear(n), offsets)
    end

    result = fill(CNUM_ZERO, n, n)
    relations = ParamRelation[]
    for component in connected_components(flow)
        if length(component) == 1
            i = only(component)
            result[i, i] = generated_scalar_exponential(flow[i, i], θ)
        elseif length(component) == 2
            exact_two_by_two_flow!(
                result, flow, component[1], component[2], θ, relations,
            )
        else
            unitary_error(
                "the adjoint action closes on a $(length(component))-dimensional coupled " *
                    "block, but no exact exponential strategy is registered for it",
            )
        end
    end
    return AffineAction(
        basis, result, fill(CNUM_ZERO, n), relations,
    )
end

function generator_structure_parts(G::QAdd)
    isempty(G.indices) || return QAdd[G]
    groups = Dict{AffineStructure, QTermDict}()
    for (term, coefficient) in G
        isempty(term.ops) && continue
        structure = operator_affine_structure(first(term.ops))
        structure === AFFINE_NONE && return QAdd[G]
        for operator in term.ops
            operator_affine_structure(operator) === structure || return QAdd[G]
        end
        terms = get!(groups, structure, QTermDict())
        addto_key!(terms, copy_key(term), coefficient)
    end
    length(groups) <= 1 && return QAdd[G]

    parts = QAdd[]
    for structure in instances(AffineStructure)
        haskey(groups, structure) && push!(parts, QAdd(groups[structure], EMPTY_INDICES))
    end
    return parts
end

function exact_generator_transform(G::QAdd, θ::Real)
    basis = closed_adjoint_basis(G)
    linear, shift = affine_commutator_data(G, basis)
    return canonical_transform(exact_closed_adjoint_action(basis, linear, shift, θ))
end

function UnitaryTransform(G::QAdd, θ::Real)
    hermitian_generator_or_throw(G)
    parts = generator_structure_parts(G)
    transform = exact_generator_transform(first(parts), θ)
    for i in 2:length(parts)
        transform = transform * exact_generator_transform(parts[i], θ)
    end
    return transform
end

UnitaryTransform(G::QSym, θ::Real) =
    UnitaryTransform(to_qadd(G), θ)

"""
    RotatingFrame(H0, t)

Frame rotating with a time-independent Hamiltonian `H0`, the interaction picture.

The transformation is ``U = e^{-i H_0 t}``. With [`transform`](@ref), a Hamiltonian
``H_0 + V`` becomes ``U^\\dagger V U``, since the motion of the frame cancels `H0`. `H0` must
be Hermitian and of a form that [`UnitaryTransform`](@ref) can exponentiate.

# Examples
```jldoctest
julia> a = Destroy(FockSpace(:cavity), :a);

julia> @variables ω g t;

julia> transform(ω * a' * a + g * (a + a'), RotatingFrame(ω * a' * a, t))
exp(-im*t*ω)*g * a + exp(im*t*ω)*g * a'
```

See also [`transform`](@ref), [`Rotation`](@ref).
"""
function RotatingFrame(H0::QAdd, t::Num)
    tt = time_or_throw(t)
    qadd_depends_on(H0, SymbolicUtils.unwrap(tt)) && unitary_error(
        "`RotatingFrame(H0, t)` requires a time-independent `H0`",
    )
    U = UnitaryTransform(H0, tt)
    return timed_transform(U, -H0, tt)
end

RotatingFrame(H0::QSym, t::Num) =
    RotatingFrame(to_qadd(H0), t)
