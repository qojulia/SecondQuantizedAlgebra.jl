const SiteKey = Tuple{Int32, Index, Int32}

"""Marker stored by a unitary transformation whose rules are time independent."""
struct StaticTime end

"""Marker stored by a unitary transformation differentiated with respect to `variable`."""
struct DynamicTime
    variable::Num
end

# Exact actions are partitioned into algebra-homogeneous affine blocks. These types live next
# to `UnitaryTransform` so every transform can carry concrete affine metadata directly. Their
# container fields are immutable by convention after construction and may be shared safely.
@enum AffineStructure::UInt8 begin
    AFFINE_BOSONIC_NAMBU
    AFFINE_SYMPLECTIC_PHASE_SPACE
    AFFINE_ORTHOGONAL
    AFFINE_UNITARY_LINEAR
end

struct AffineBlock
    structure::AffineStructure
    basis::Vector{Op}
    linear::Matrix{CNum}
    shift::Vector{CNum}
end

struct AffineAction
    blocks::Vector{AffineBlock}
    relations::Vector{ParamRelation}
end

"""
    UnitaryTransform

An exact change of frame compiled from a canonical affine action on a complete set of site
generators. Construct transforms with [`Displace`](@ref), [`Rotation`](@ref),
[`Squeeze`](@ref), [`Bogoliubov`](@ref), or the frame constructors; apply them with
[`conjugate`](@ref) or [`transform`](@ref).
"""
struct UnitaryTransform{T}
    action::AffineAction
    rules::Dict{Op, QAdd}
    generators::Vector{Op}
    gauge::QAdd
    time::T

    function UnitaryTransform{T}(
            action::AffineAction, rules::Dict{Op, QAdd}, generators::Vector{Op},
            gauge::QAdd, time::T,
        ) where {T}
        (T === StaticTime || T === DynamicTime) ||
            throw(ArgumentError("invalid unitary-transform time marker `$T`"))
        return new{T}(action, rules, generators, gauge, time)
    end
end

@noinline unitary_error(message::AbstractString) = throw(ArgumentError(message))

is_fock(o::Op) = is_destroy(o) || is_create(o)
is_phase_space(o::Op) = is_position(o) || is_momentum(o)
lowering(o::Op) = is_create(o) ? adjoint(o) : o

function fock_or_throw(o::Op, what::AbstractString)
    is_fock(o) || unitary_error("$what expects a Fock ladder operator, got $(o.kind)")
    return lowering(o)
end

function site_generators(o::Op)
    if is_fock(o)
        d = lowering(o)
        return Op[d, adjoint(d)]
    elseif is_pauli(o) || is_spin(o)
        return Op[
            Op(o.kind, o.name_id, o.space_index, o.index, Int32(axis), 0, 0, 0)
                for axis in 1:3
        ]
    elseif is_transition(o)
        n = Int(o.nlev)
        return Op[
            Op(OP_TRANSITION, o.name_id, o.space_index, o.index, Int32(i), Int32(j), o.g, o.nlev)
                for i in 1:n for j in 1:n
        ]
    end
    return Op[]
end

function validate_complete(generators::Vector{Op})
    available = Set(generators)
    checked = Set{SiteKey}()
    for first_generator in generators
        key = site_key(first_generator)
        key in checked && continue
        push!(checked, key)

        expected = site_generators(first_generator)
        if isempty(expected)
            if is_phase_space(first_generator)
                has_x = false
                has_p = false
                for generator in generators
                    site_key(generator) == key || continue
                    has_x |= is_position(generator)
                    has_p |= is_momentum(generator)
                end
                (has_x && has_p) || unitary_error(
                    "incomplete rule set: `$first_generator` has no rule for its conjugate variable",
                )
            end
            continue
        end
        for generator in expected
            generator in available || unitary_error(
                "incomplete rule set: `$first_generator` is covered but `$generator` is not",
            )
        end
    end
    return nothing
end

function validated_transform(
        action::AffineAction, gauge::QAdd, time::T,
    ) where {T <: Union{StaticTime, DynamicTime}}
    rules = affine_rules(action)
    generators = sort!(collect(keys(rules)))
    for generator in generators
        (has_index(generator.index) && index_slot(generator.index) === nothing) &&
            unitary_error(
            "unitary transforms of free indexed-operator families are not part of " *
                "the exact closed-form API; resolve the index to one site first",
        )
    end
    validate_complete(generators)
    return UnitaryTransform{T}(action, rules, generators, gauge, time)
end

function time_or_throw(t::Num)
    raw = SymbolicUtils.unwrap(t)
    SymbolicUtils.issym(raw) || unitary_error(
        "time must be a real symbolic variable, got `$t`",
    )
    return t
end

function timed_transform(U::UnitaryTransform{StaticTime}, gauge::QAdd, t::Num)
    time = DynamicTime(time_or_throw(t))
    reduced = reduce_params(gauge, U.action.relations, true)
    return UnitaryTransform{DynamicTime}(U.action, U.rules, U.generators, reduced, time)
end

function covered_site(U::UnitaryTransform, key::SiteKey)
    for generator in U.generators
        site_key(generator) == key && return true
    end
    return false
end

function validate_coverage(q::QAdd, U::UnitaryTransform)
    for (term, _) in q, operator in term.ops
        haskey(U.rules, operator) && continue
        covered_site(U, site_key(operator)) || continue
        unitary_error(
            "`$operator` acts on a site covered by this transform but has no rule; " *
                "the constructor must cover every generator of a transformed site",
        )
    end
    return nothing
end

apply_rules(q::QAdd, rules::Dict{Op, QAdd}) = substitute_op_rules(q, rules)

function reduce_params(q::QAdd, relations::Vector{ParamRelation}, gated::Bool)
    isempty(relations) && return q
    scratch = ParamRelation[]
    return map_coefficients(c -> reduce_all(c, relations, gated, scratch), q)
end

"""
    conjugate(A, U::UnitaryTransform)

Return the observable change of frame `U' * A * U`. For a Hamiltonian in a moving frame,
use [`transform`](@ref) to include the time-dependent gauge term.
"""
function conjugate(q::QAdd, U::UnitaryTransform)
    validate_coverage(q, U)
    return reduce_params(apply_rules(q, U.rules), U.action.relations, true)
end

conjugate(o::QSym, U::UnitaryTransform) =
    conjugate(single_qadd(CNUM_ONE, Op[o]), U)

"""
    transform(H, U::UnitaryTransform)

Return `U' * H * U + im*(∂ₜU')*U`. For a static transform this is exactly
[`conjugate`](@ref).
"""
transform(q::QAdd, U::UnitaryTransform{StaticTime}) = conjugate(q, U)
transform(q::QAdd, U::UnitaryTransform{DynamicTime}) = conjugate(q, U) + U.gauge
transform(o::QSym, U::UnitaryTransform) =
    transform(single_qadd(CNUM_ONE, Op[o]), U)

"""Return the Hamiltonian gauge term stored by `U`."""
gauge_term(U::UnitaryTransform) = U.gauge

"""Return the fundamental operators transformed by `U`, in canonical order."""
generators(U::UnitaryTransform) = copy(U.generators)

function Base.inv(U::UnitaryTransform{T}) where {T}
    inverse_action = canonical_affine_inverse(U.action)
    inverse_rules = affine_rules(inverse_action)
    relations = U.action.relations
    gauge = if T === StaticTime || iszero(U.gauge)
        U.gauge
    else
        -reduce_params(apply_rules(U.gauge, inverse_rules), relations, true)
    end
    return UnitaryTransform{T}(
        inverse_action, inverse_rules, U.generators, gauge, U.time,
    )
end

Base.adjoint(U::UnitaryTransform) = inv(U)

function contains_relation(relations::Vector{ParamRelation}, relation::ParamRelation)
    for existing in relations
        isequal(existing.hi, relation.hi) && isequal(existing.lo, relation.lo) &&
            existing.sign == relation.sign && return true
    end
    return false
end

function merge_relations(a::Vector{ParamRelation}, b::Vector{ParamRelation})
    isempty(a) && return b
    isempty(b) && return a

    first_new = 0
    for i in eachindex(b)
        if !contains_relation(a, b[i])
            first_new = i
            break
        end
    end
    iszero(first_new) && return a

    out = copy(a)
    for i in first_new:lastindex(b)
        relation = b[i]
        contains_relation(out, relation) || push!(out, relation)
    end
    return out
end

function compose_rule_image(image::QAdd, rules::Dict{Op, QAdd})
    isempty(image.indices) || return apply_rules(image, rules)
    out = QTermDict()
    for (term, coefficient) in image
        if isempty(term.ops)
            addto_key!(out, copy_key(term), coefficient)
        elseif length(term.ops) == 1 && haskey(rules, first(term.ops))
            replacement = rules[first(term.ops)]
            for (replacement_term, replacement_coefficient) in replacement
                addto_key!(
                    out, copy_key(replacement_term),
                    mul_cnum(coefficient, replacement_coefficient),
                )
            end
        elseif length(term.ops) == 1
            addto_key!(out, copy_key(term), coefficient)
        else
            return apply_rules(image, rules)
        end
    end
    return QAdd(out, EMPTY_INDICES)
end

function compose_rules(first::Dict{Op, QAdd}, second::Dict{Op, QAdd})
    out = Dict{Op, QAdd}()
    sizehint!(out, length(first) + length(second))
    for (generator, image) in first
        out[generator] = compose_rule_image(image, second)
    end
    for (generator, image) in second
        haskey(out, generator) || (out[generator] = image)
    end
    return out
end

function add_gauges(left::QAdd, right::QAdd)::QAdd
    if isempty(left.indices) && isempty(right.indices) &&
            length(left.arguments) == length(right.arguments) == 1
        first_term, first_coefficient = first(left.arguments)
        second_term, second_coefficient = first(right.arguments)
        if isequal(first_term, second_term)
            coefficient = add_cnum(first_coefficient, second_coefficient)
            iszero_cnum(coefficient) && return zero_qadd()
            return QAdd(QTermDict(first_term => coefficient), EMPTY_INDICES)
        end
    end
    return left + right
end

function coefficient_depends_on(c::CNum, variable)
    tail = c.tail
    tail isa Native && return false
    if tail isa Poly
        for monomial in tail.terms, factor in monomial.syms
            raw_depends_on(factor, variable) && return true
        end
        return false
    end
    return raw_depends_on(tail.expr, variable)
end

function rules_depend_on(U::UnitaryTransform{StaticTime}, t::Num)
    variable = SymbolicUtils.unwrap(t)
    for image in values(U.rules), (_, coefficient) in image
        coefficient_depends_on(coefficient, variable) && return true
    end
    return false
end

function check_adopted_time(U::UnitaryTransform{StaticTime}, t::Num)
    rules_depend_on(U, t) && unitary_error(
        "a static transform whose rules depend on `$t` cannot be composed with a timed " *
            "transform; construct the moving transform with its timed constructor",
    )
    return nothing
end

function compose_action(
        first::AffineAction, second::AffineAction, relations::Vector{ParamRelation},
    )
    if length(first.blocks) == 1 && length(second.blocks) == 1
        first_block = only(first.blocks)
        second_block = only(second.blocks)
        if blocks_overlap(first_block, second_block)
            return AffineAction(
                AffineBlock[compose_overlapping_block(first_block, second_block, relations)],
                relations,
            )
        end
        return AffineAction(AffineBlock[first_block, second_block], relations)
    end
    return compose_action_metadata(first, second, relations)
end

function compose(
        first::UnitaryTransform, second::UnitaryTransform, time::T,
    ) where {T <: Union{StaticTime, DynamicTime}}
    relations = merge_relations(first.action.relations, second.action.relations)
    action = compose_action(first.action, second.action, relations)
    rules = compose_rules(first.rules, second.rules)

    gauge = if iszero(first.gauge)
        second.gauge
    else
        transported = reduce_params(
            apply_rules(first.gauge, second.rules), relations, true,
        )
        iszero(second.gauge) ? transported : add_gauges(transported, second.gauge)
    end

    if length(rules) == length(first.rules)
        same_layout = all(generator -> haskey(rules, generator), first.generators)
        generators = same_layout ? first.generators : sort!(collect(keys(rules)))
    else
        generators = sort!(collect(keys(rules)))
    end
    return UnitaryTransform{T}(action, rules, generators, gauge, time)
end

Base.:*(
    first::UnitaryTransform{StaticTime}, second::UnitaryTransform{StaticTime},
) = compose(first, second, StaticTime())

function Base.:*(
        first::UnitaryTransform{StaticTime}, second::UnitaryTransform{DynamicTime},
    )
    t = second.time.variable
    check_adopted_time(first, t)
    return compose(first, second, second.time)
end

function Base.:*(
        first::UnitaryTransform{DynamicTime}, second::UnitaryTransform{StaticTime},
    )
    t = first.time.variable
    check_adopted_time(second, t)
    return compose(first, second, first.time)
end

function Base.:*(
        first::UnitaryTransform{DynamicTime}, second::UnitaryTransform{DynamicTime},
    )
    isequal(first.time.variable, second.time.variable) || unitary_error(
        "cannot compose transforms with different time variables " *
            "(`$(first.time.variable)` and `$(second.time.variable)`)",
    )
    return compose(first, second, first.time)
end

# Rule-building primitives used by affine compilation and gauge construction.
function rule_qadd(pairs::Vector{Tuple{CNum, Vector{Op}}})
    out = QTermDict()
    for (coefficient, operators) in pairs
        iszero_cnum(coefficient) || addto!(out, operators, coefficient)
    end
    return QAdd(out, EMPTY_INDICES)
end

rule_qadd(pairs::Vararg{Tuple{CNum, Vector{Op}}}) =
    rule_qadd(Tuple{CNum, Vector{Op}}[pairs...])
scaled(coefficient::CNum, operator::Op) = rule_qadd((coefficient, Op[operator]))

phase(ϕ::Real) = phase_coeff(as_num(ϕ))
conj_phase(ϕ::Real) = conj_cnum(phase(ϕ))
trig_rel(θ::Real) = ParamRelation(cos(θ), sin(θ), -1)
hyp_rel(r::Real) = ParamRelation(cosh(r), sinh(r), 1)

dt(c::CNum, t::Num) = Symbolics.derivative(c, t)
dt(x::Coefficient, t::Num) = dt(to_cnum(x), t)

gauge(generator::QAdd, θ::Real, t::Num) = generator * neg_cnum(dt(θ, t))
