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

An exact change of frame represented by its forward and inverse action on a complete set of
site generators. Construct transforms with [`Displace`](@ref), [`Rotation`](@ref), or
[`Squeeze`](@ref); apply them with [`conjugate`](@ref) or [`transform`](@ref).
"""
struct UnitaryTransform{T}
    rules::Dict{Op, QAdd}
    inverse_rules::Dict{Op, QAdd}
    generators::Vector{Op}
    sites::Vector{SiteInfo}
    gauge::QAdd
    time::T
    relations::Vector{ParamRelation}

    function UnitaryTransform{T}(
            rules::Dict{Op, QAdd}, inverse_rules::Dict{Op, QAdd},
            generators::Vector{Op}, sites::Vector{SiteInfo}, gauge::QAdd, time::T,
            relations::Vector{ParamRelation}, ::Val{:validated},
        ) where {T}
        (T === StaticTime || T === DynamicTime) ||
            throw(ArgumentError("invalid unitary-transform time marker `$T`"))
        return new{T}(rules, inverse_rules, generators, sites, gauge, time, relations)
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
        relations::Vector{ParamRelation} = ParamRelation[],
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
    return UnitaryTransform{T}(
        rules, inverse_rules, generators, sites, gauge, time, usable, Val(:validated),
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

function timed_transform(U::UnitaryTransform{StaticTime}, gauge::QAdd, t::Num)
    time = DynamicTime(time_or_throw(t))
    reduced = reduce_params(gauge, U.relations, true)
    return UnitaryTransform{DynamicTime}(
        U.rules, U.inverse_rules, U.generators, U.sites, reduced, time, U.relations,
        Val(:validated),
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

function Base.inv(U::UnitaryTransform{T}) where {T}
    gauge = if T === StaticTime || iszero(U.gauge)
        U.gauge
    else
        -reduce_params(apply_rules(U.gauge, U.inverse_rules), U.relations, true)
    end
    return UnitaryTransform{T}(
        copy(U.inverse_rules), copy(U.rules), U.generators, U.sites, gauge, U.time,
        copy(U.relations), Val(:validated),
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

# The first transform is applied first. Images belonging only to the second transform are
# preserved, while overlapping images are transported through the second map.
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

function diagonal_entry(
        image::QAdd, generator::Op,
    )::Union{Nothing, Tuple{QTerm, Coeff}}
    isempty(image.indices) || return nothing
    length(image.arguments) == 1 || return nothing
    term, coefficient = first(image.arguments)
    length(term.ops) == 1 || return nothing
    only(term.ops) == generator || return nothing
    return (term, coefficient)
end

function diagonal_coefficient(image::QAdd, generator::Op)::Union{Coeff, Nothing}
    entry = diagonal_entry(image, generator)
    return entry === nothing ? nothing : entry[2]
end

diagonal_rule(term::QTerm, coefficient::Coeff) =
    QAdd(QTermDict(term => coefficient), EMPTY_INDICES)

function coefficients_are_inverse(left::Coeff, right::Coeff)::Bool
    if left.tail isa Native && right.tail isa Native
        return isone(left.z * right.z)
    end
    (left.tail isa Poly && right.tail isa Poly) || return false
    length(left.tail.terms) == length(right.tail.terms) == 1 || return false
    left_term = only(left.tail.terms)
    right_term = only(right.tail.terms)
    isone(left_term.scalar * right_term.scalar) || return false
    length(left_term.syms) == length(right_term.syms) || return false
    @inbounds for i in eachindex(left_term.syms)
        left_term.syms[i] === right_term.syms[i] || return false
        left_term.exps[i] == -right_term.exps[i] || return false
    end
    return true
end

function compose_diagonal_rules(
        first::UnitaryTransform, second::UnitaryTransform,
    )::Union{Nothing, Tuple{Dict{Op, QAdd}, Dict{Op, QAdd}}}
    length(first.rules) == length(second.rules) || return nothing
    rules = Dict{Op, QAdd}()
    inverse_rules = Dict{Op, QAdd}()
    sizehint!(rules, length(first.rules))
    sizehint!(inverse_rules, length(first.rules))
    for generator in first.generators
        first_entry = diagonal_entry(first.rules[generator], generator)
        first_entry === nothing && return nothing
        term, first_coefficient = first_entry
        first_inverse_coefficient = diagonal_coefficient(
            first.inverse_rules[generator], generator,
        )
        first_inverse_coefficient === nothing && return nothing
        coefficients_are_inverse(first_coefficient, first_inverse_coefficient) ||
            return nothing
        second_image = get(second.rules, generator, nothing)
        second_image === nothing && return nothing
        second_coefficient = diagonal_coefficient(second_image, generator)
        second_coefficient === nothing && return nothing
        second_inverse_coefficient = diagonal_coefficient(
            second.inverse_rules[generator], generator,
        )
        second_inverse_coefficient === nothing && return nothing
        coefficients_are_inverse(second_coefficient, second_inverse_coefficient) ||
            return nothing
        paired_generator = adjoint(generator)
        paired_rule = get(rules, paired_generator, nothing)
        if paired_rule !== nothing
            paired_first_coefficient = diagonal_coefficient(
                first.rules[paired_generator], paired_generator,
            )
            paired_second_coefficient = diagonal_coefficient(
                second.rules[paired_generator], paired_generator,
            )
            paired_first_coefficient === nothing && return nothing
            paired_second_coefficient === nothing && return nothing
            coefficients_are_inverse(first_coefficient, paired_first_coefficient) ||
                return nothing
            coefficients_are_inverse(second_coefficient, paired_second_coefficient) ||
                return nothing
            paired_inverse_rule = inverse_rules[paired_generator]
            coefficient = diagonal_coefficient(paired_inverse_rule, paired_generator)
            inverse_coefficient = diagonal_coefficient(paired_rule, paired_generator)
            coefficient === nothing && return nothing
            inverse_coefficient === nothing && return nothing
            rules[generator] = diagonal_rule(term, coefficient)
            inverse_rules[generator] = diagonal_rule(term, inverse_coefficient)
            continue
        end
        coefficient = mul_cnum(first_coefficient, second_coefficient)
        rules[generator] = diagonal_rule(term, coefficient)
        inverse_rules[generator] = diagonal_rule(term, inv(coefficient))
    end
    return (rules, inverse_rules)
end

function invariant_under_diagonal_rules(gauge::QAdd, rules::Dict{Op, QAdd})::Bool
    isempty(gauge.indices) || return false
    for (term, _) in gauge
        for operator in term.ops
            image = get(rules, operator, nothing)
            image === nothing && continue
            coefficient = diagonal_coefficient(image, operator)
            coefficient === nothing && return false
            paired = adjoint(operator)
            if paired == operator
                isequal(coefficient, CNUM_ONE) || return false
            else
                paired_image = get(rules, paired, nothing)
                paired_image === nothing && return false
                paired_coefficient = diagonal_coefficient(paired_image, paired)
                paired_coefficient === nothing && return false
                coefficients_are_inverse(coefficient, paired_coefficient) || return false
                count(==(operator), term.ops) == count(==(paired), term.ops) || return false
            end
        end
    end
    return true
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

function compose(
        first::UnitaryTransform, second::UnitaryTransform, time::T,
    ) where {T <: Union{StaticTime, DynamicTime}}
    relations = merge_relations(first.relations, second.relations)
    diagonal_rules = compose_diagonal_rules(first, second)
    rules, inverse_rules = diagonal_rules === nothing ?
        (
            compose_rules(first.rules, second.rules),
            compose_rules(second.inverse_rules, first.inverse_rules),
        ) : diagonal_rules
    gauge = if iszero(first.gauge)
        second.gauge
    elseif invariant_under_diagonal_rules(first.gauge, second.rules)
        iszero(second.gauge) ? first.gauge : add_gauges(first.gauge, second.gauge)
    else
        transported = reduce_params(
            apply_rules(first.gauge, second.rules), relations, true,
        )
        iszero(second.gauge) ? transported : add_gauges(transported, second.gauge)
    end
    if length(rules) == length(first.rules)
        same_layout = true
        for generator in first.generators
            haskey(rules, generator) || (same_layout = false; break)
        end
        generators = same_layout ? first.generators : sort!(collect(keys(rules)))
        sites = same_layout ? first.sites : site_infos(generators)
    else
        generators = sort!(collect(keys(rules)))
        sites = site_infos(generators)
    end
    return UnitaryTransform{T}(
        rules, inverse_rules, generators, sites, gauge, time, relations, Val(:validated),
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

# Rule-building primitives. All named constructors produce linear images, so these avoid a
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
