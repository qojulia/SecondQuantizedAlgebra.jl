const SiteKey = Tuple{Int32, Index, Int32}

"""Marker stored by a unitary transformation whose rules are time independent."""
struct StaticTime end

"""Marker stored by a unitary transformation differentiated with respect to `variable`."""
struct DynamicTime
    variable::Num
end

struct SiteInfo
    key::SiteKey
    generators::Vector{Op}
end

"""
    UnitaryTransform

An exact change of frame compiled from a canonical action on a complete set of site
generators. Construct transforms with [`Displace`](@ref), [`Rotation`](@ref),
[`Squeeze`](@ref), [`Bogoliubov`](@ref), or the frame constructors; apply them with
[`conjugate`](@ref) or [`transform`](@ref).
"""
struct UnitaryTransform{T, A}
    action::A
    rules::Dict{Op, QAdd}
    inverse_rules::Dict{Op, QAdd}
    generators::Vector{Op}
    sites::Vector{SiteInfo}
    gauge::QAdd
    time::T
    relations::Vector{ParamRelation}

    function UnitaryTransform{T, A}(
            action::A, rules::Dict{Op, QAdd}, inverse_rules::Dict{Op, QAdd},
            generators::Vector{Op}, sites::Vector{SiteInfo}, gauge::QAdd, time::T,
            relations::Vector{ParamRelation}, ::Val{:validated},
        ) where {T, A}
        (T === StaticTime || T === DynamicTime) ||
            throw(ArgumentError("invalid unitary-transform time marker `$T`"))
        return new{T, A}(
            action, rules, inverse_rules, generators, sites, gauge, time, relations,
        )
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

function site_infos(generators::Vector{Op})
    sites = SiteInfo[]
    for g in generators
        key = site_key(g)
        found = findfirst(site -> site.key == key, sites)
        if found === nothing
            push!(sites, SiteInfo(key, Op[g]))
        else
            push!(sites[found].generators, g)
        end
    end
    sort!(sites; by = site -> (site.key[1], index_key(site.key[2]), name_rank(site.key[3])))
    return sites
end

function validate_complete(sites::Vector{SiteInfo})
    for site in sites
        first_generator = first(site.generators)
        expected = site_generators(first_generator)
        if isempty(expected)
            if is_phase_space(first_generator)
                has_x = any(is_position, site.generators)
                has_p = any(is_momentum, site.generators)
                (has_x && has_p) || unitary_error(
                    "incomplete rule set: `$first_generator` has no rule for its conjugate variable",
                )
            end
            continue
        end
        available = Set(site.generators)
        for generator in expected
            generator in available || unitary_error(
                "incomplete rule set: `$first_generator` is covered but `$generator` is not",
            )
        end
    end
    return nothing
end

function validated_transform(
        rules::Dict{Op, QAdd}, inverse_rules::Dict{Op, QAdd}, gauge::QAdd, time::T,
        relations::Vector{ParamRelation} = ParamRelation[], action = nothing,
    ) where {T <: Union{StaticTime, DynamicTime}}
    isempty(rules) && unitary_error("a `UnitaryTransform` needs at least one rule")
    length(rules) == length(inverse_rules) || unitary_error(
        "forward and inverse rules must cover the same generators",
    )
    generators = sort!(collect(keys(rules)))
    for generator in generators
        (has_index(generator.index) && index_slot(generator.index) === nothing) &&
            unitary_error(
            "unitary transforms of free indexed-operator families are not part of " *
                "the exact closed-form API; resolve the index to one site first",
        )
        haskey(inverse_rules, generator) || unitary_error(
            "inverse rules are missing the generator `$generator`",
        )
    end
    sites = site_infos(generators)
    validate_complete(sites)
    validate_complete(site_infos(sort!(collect(keys(inverse_rules)))))
    usable = all(is_usable_rel, relations) ? relations : filter(is_usable_rel, relations)
    return UnitaryTransform{T, typeof(action)}(
        action, rules, inverse_rules, generators, sites, gauge, time, usable,
        Val(:validated),
    )
end

static_transform(
    rules::Dict{Op, QAdd}, inverse_rules::Dict{Op, QAdd},
    relations::Vector{ParamRelation} = ParamRelation[],
) = validated_transform(rules, inverse_rules, zero_qadd(), StaticTime(), relations)

function time_or_throw(t::Num)
    raw = SymbolicUtils.unwrap(t)
    SymbolicUtils.issym(raw) || unitary_error(
        "time must be a real symbolic variable, got `$t`",
    )
    return t
end

function timed_transform(U::UnitaryTransform{StaticTime, A}, gauge::QAdd, t::Num) where {A}
    time = DynamicTime(time_or_throw(t))
    reduced = reduce_params(gauge, U.relations, true)
    return UnitaryTransform{DynamicTime, A}(
        U.action, U.rules, U.inverse_rules, U.generators, U.sites, reduced, time,
        U.relations, Val(:validated),
    )
end

function covered_site(U::UnitaryTransform, key::SiteKey)
    for site in U.sites
        site.key == key && return true
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
    return reduce_params(apply_rules(q, U.rules), U.relations, true)
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

compiled_inverse_action_metadata(::Nothing) = nothing
compiled_inverse_action_metadata(action) = nothing
compile_composed_action_metadata(action) = nothing

function Base.inv(U::UnitaryTransform{T, A}) where {T, A}
    gauge = if T === StaticTime || iszero(U.gauge)
        U.gauge
    else
        -reduce_params(apply_rules(U.gauge, U.inverse_rules), U.relations, true)
    end
    inverse_action = compiled_inverse_action_metadata(U.action)::A
    return UnitaryTransform{T, A}(
        inverse_action, copy(U.inverse_rules), copy(U.rules), U.generators, U.sites,
        gauge, U.time, copy(U.relations), Val(:validated),
    )
end

Base.adjoint(U::UnitaryTransform) = inv(U)

function merge_relations(a::Vector{ParamRelation}, b::Vector{ParamRelation})
    isempty(a) && return b
    isempty(b) && return a
    out = copy(a)
    for relation in b
        any(
            r -> isequal(r.hi, relation.hi) && isequal(r.lo, relation.lo) &&
                r.sign == relation.sign, out
        ) || push!(out, relation)
    end
    return out
end

# Legacy rule composition is retained only for transforms that do not carry affine metadata.
# Every named exact total transform is expected to use the affine path.
function compose_rule_image(image::QAdd, rules::Dict{Op, QAdd})
    isempty(image.indices) || return apply_rules(image, rules)
    out = QTermDict()
    for (term, coefficient) in image
        if isempty(term.ops)
            addto_key!(out, copy_key(term), coefficient)
        elseif length(term.ops) == 1
            replacement = get(rules, first(term.ops), nothing)
            if replacement === nothing
                addto_key!(out, copy_key(term), coefficient)
            else
                for (replacement_term, replacement_coefficient) in replacement
                    addto_key!(
                        out, copy_key(replacement_term),
                        mul_cnum(coefficient, replacement_coefficient),
                    )
                end
            end
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

compose_action_metadata(first, second, relations::Vector{ParamRelation}) = nothing

function compose(
        first::UnitaryTransform, second::UnitaryTransform, time::T,
    ) where {T <: Union{StaticTime, DynamicTime}}
    relations = merge_relations(first.relations, second.relations)
    action = compose_action_metadata(first.action, second.action, relations)
    compiled = compile_composed_action_metadata(action)
    rules, inverse_rules = if compiled === nothing
        (
            compose_rules(first.rules, second.rules),
            compose_rules(second.inverse_rules, first.inverse_rules),
        )
    else
        compiled
    end

    gauge = if iszero(first.gauge)
        second.gauge
    else
        transported = reduce_params(
            apply_rules(first.gauge, second.rules), relations, true,
        )
        iszero(second.gauge) ? transported : add_gauges(transported, second.gauge)
    end

    generators = sort!(collect(keys(rules)))
    sites = site_infos(generators)
    return UnitaryTransform{T, typeof(action)}(
        action, rules, inverse_rules, generators, sites, gauge, time, relations,
        Val(:validated),
    )
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

# Rule-building primitives. Named constructors produce linear images, so these avoid a
# general expression canonicalization pass during construction.
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
with_adjoint(generator::Op, image::QAdd) =
    Dict{Op, QAdd}(generator => image, adjoint(generator) => adjoint(image))
pair_rules(x::Op, p::Op, fx::QAdd, fp::QAdd) = Dict{Op, QAdd}(x => fx, p => fp)

phase(ϕ::Real) = phase_coeff(as_num(ϕ))
conj_phase(ϕ::Real) = conj_cnum(phase(ϕ))
trig_rel(θ::Real) = ParamRelation(cos(θ), sin(θ), -1)
hyp_rel(r::Real) = ParamRelation(cosh(r), sinh(r), 1)

dt(c::CNum, t::Num) = Symbolics.derivative(c, t)
dt(x::Coefficient, t::Num) = dt(to_cnum(x), t)

gauge(generator::QAdd, θ::Real, t::Num) = generator * neg_cnum(dt(θ, t))

include("unitary_constructors.jl")
