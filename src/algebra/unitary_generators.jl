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

    leading = operator_affine_structure(first(raw))
    shared = all(o -> operator_affine_structure(o) === leading, raw)
    return canonical_block_basis(shared ? leading : AFFINE_NONE, raw)
end

function add_adjoint_entry!(
        linear::Matrix{CNum}, shift::Vector{CNum}, row::Int, locations::Dict{Op, Int},
        ops::Vector{Op}, weight::CNum,
    )
    if isempty(ops)
        shift[row] = shift[row] + weight
        return true
    end
    length(ops) == 1 || return false
    column = get(locations, only(ops), 0)
    column == 0 && return false
    linear[row, column] = linear[row, column] + weight
    return true
end

function affine_commutator_data(G::QAdd, basis::Vector{Op})
    n = length(basis)
    locations = Dict{Op, Int}(operator => i for (i, operator) in enumerate(basis))
    site_rows = Dict{SiteKey, Vector{Int}}()
    for (row, generator) in enumerate(basis)
        push!(get!(Vector{Int}, site_rows, site_key(generator)), row)
    end
    all_rows = collect(1:n)
    indexed = any(generator -> has_index(generator.index), basis)
    no_rows = Int[]
    linear = fill(CNUM_ZERO, n, n)
    shift = fill(CNUM_ZERO, n)
    residuals = Dict{Int, QAdd}()
    for (term, coefficient) in G, (position, factor) in enumerate(term.ops)
        rows = (indexed || has_index(factor.index)) ? all_rows : get(site_rows, site_key(factor), no_rows)
        for row in rows, (inner, scale) in commutator(factor, basis[row])
            ops = vcat(view(term.ops, 1:(position - 1)), inner.ops, view(term.ops, (position + 1):length(term.ops)))
            weight = coefficient * scale
            add_adjoint_entry!(linear, shift, row, locations, ops, weight) ||
                (residuals[row] = get(residuals, row, zero_qadd()) + weight * prod(ops))
        end
    end
    for row in sort!(collect(keys(residuals)))
        generator = basis[row]
        for (term, coefficient) in simplify(residuals[row])
            add_adjoint_entry!(linear, shift, row, locations, term.ops, coefficient) && continue
            length(term.ops) == 1 && unitary_error(
                "the adjoint action leaves the finite generator basis: " *
                    "`[$G, $generator]` contains `$(only(term.ops))`",
            )
            unitary_error(
                "the adjoint action on `$generator` is not affine-linear; got " *
                    "`$(sprint(show, single_qadd(coefficient, copy(term.ops))))`",
            )
        end
    end
    for entries in (linear, shift), k in eachindex(entries)
        is_native(entries[k]) || (entries[k] = simplify(entries[k]))
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
        coefficient_equal(flow, CNUM_NEG1) && return inv(to_cnum(exp(θ)))
        return to_cnum(exp(exact_real_argument(flow, θ)))
    end

    phase_rate = CNUM_NEG_IM * flow
    if coefficient_is_real(phase_rate)
        coefficient_equal(phase_rate, CNUM_ONE) && return phase(θ)
        coefficient_equal(phase_rate, CNUM_NEG1) && return conj_phase(θ)
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

function exact_rotation_pair!(relations::Vector{ParamRelation}, rate::CNum, θ::Real)
    coefficient_equal(rate, CNUM_ONE) && return trig_pair!(relations, θ)
    if coefficient_equal(rate, CNUM_NEG1)
        cosine, sine = trig_pair!(relations, θ)
        return cosine, -sine
    end
    return trig_pair!(relations, exact_real_argument(rate, θ))
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

function relation_of(relation::ParamRelation, sign::Int8, angle)
    relation.sign == sign || return false
    hi = SymbolicUtils.unwrap(relation.hi)
    return SymbolicUtils.iscall(hi) && isequal(only(SymbolicUtils.arguments(hi)), angle)
end

function relation_pair!(relations::Vector{ParamRelation}, angle::Real, build::F, sign::Int8) where {F}
    raw = SymbolicUtils.unwrap(angle)
    found = position_of(relation -> relation_of(relation, sign, raw), relations)
    relation = iszero(found) ? build(angle) : relations[found]
    iszero(found) && push_relation!(relations, relation)
    return to_cnum(relation.hi), to_cnum(relation.lo)
end

trig_pair!(relations::Vector{ParamRelation}, angle::Real) = relation_pair!(relations, angle, trig_rel, Int8(-1))

hyperbolic_pair!(relations::Vector{ParamRelation}, angle::Real) = relation_pair!(relations, angle, hyp_rel, Int8(1))

function exact_two_by_two_flow(upper::CNum, lower::CNum, θ::Real, relations::Vector{ParamRelation})
    if coefficient_is_real(upper) && coefficient_is_real(lower) &&
            coefficient_equal(lower, -upper)
        cosine, sine = exact_rotation_pair!(relations, upper, θ)
        return cosine, sine, -sine
    end

    if coefficient_equal(lower, conj(upper))
        if coefficient_equal(upper * lower, CNUM_ONE)
            ch, sh = hyperbolic_pair!(relations, θ)
            return ch, upper * sh, lower * sh
        elseif coefficient_is_real(upper) && coefficient_equal(lower, upper)
            ch, sh = hyperbolic_pair!(relations, exact_real_argument(upper, θ))
            return ch, sh, sh
        end
        terms = phase_terms(upper)
        if is_real_scaled_phase(terms)
            direction = phase(only(terms).phase)
            ch, sh = hyperbolic_pair!(relations, exact_real_argument(only(terms).amplitude, θ))
            return ch, direction * sh, conj(direction) * sh
        end
    end

    sign_of_square, rate = rational_square_root(upper * lower)
    if !iszero(sign_of_square)
        angle = exact_real_argument(to_cnum(rate), θ)
        diagonal, off_diagonal = sign_of_square < 0 ? trig_pair!(relations, angle) :
            hyperbolic_pair!(relations, angle)
        scaled = to_cnum(inv(rate)) * off_diagonal
        return diagonal, upper * scaled, lower * scaled
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
    exponentials = Dict{CNum, CNum}()
    rotations = Dict{Tuple{CNum, CNum}, Tuple{CNum, CNum, CNum}}()
    for component in connected_components(flow)
        if length(component) == 1
            k = only(component)
            rate = flow[k, k]
            result[k, k] = get(exponentials, rate) do
                haskey(exponentials, -rate) ? inv(exponentials[-rate]) : generated_scalar_exponential(rate, θ)
            end
            exponentials[rate] = result[k, k]
        elseif length(component) == 2
            i, j = component[1], component[2]
            iszero(flow[i, i]) && iszero(flow[j, j]) || unitary_error(
                "exact two-dimensional adjoint blocks currently require zero diagonal entries",
            )
            rising, falling = flow[i, j], flow[j, i]
            diagonal, upper, lower = get!(rotations, (rising, falling)) do
                exact_two_by_two_flow(rising, falling, θ, relations)
            end
            result[i, i] = diagonal
            result[i, j] = upper
            result[j, i] = lower
            result[j, j] = diagonal
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
