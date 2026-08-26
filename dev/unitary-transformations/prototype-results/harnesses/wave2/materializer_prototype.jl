module MaterializerPrototype

using SecondQuantizedAlgebra
const SQA = SecondQuantizedAlgebra

const CNum = SQA.CNum
const Coeff = SQA.Coeff
const Op = SQA.Op
const QAdd = SQA.QAdd
const QTerm = SQA.QTerm
const QTermDict = SQA.QTermDict
const ParamRelation = SQA.ParamRelation

"""Dense closure-coordinate images for the public generator rows."""
struct CoefficientRuleMap
    generators::Vector{Op}
    coordinates::Vector{QTerm}
    coefficients::Matrix{Coeff}
    offsets::Vector{Coeff}
end

struct AffineRulePair
    forward::CoefficientRuleMap
    inverse::CoefficientRuleMap
end

"""Candidate nested storage. This is deliberately not retained by the materializer."""
struct RuleTransformData
    rules::Dict{Op, QAdd}
    inverse_rules::Dict{Op, QAdd}
    generators::Vector{Op}
    sites::Vector{SQA.SiteInfo}
    relations::Vector{ParamRelation}
end

struct NestedTransform{T}
    data::RuleTransformData
    gauge::QAdd
    time::T
end

function _validate_map(map::CoefficientRuleMap)
    nrows, ncols = size(map.coefficients)
    nrows == length(map.generators) || throw(DimensionMismatch("generator rows"))
    ncols == length(map.coordinates) || throw(DimensionMismatch("coordinate columns"))
    nrows == length(map.offsets) || throw(DimensionMismatch("offset rows"))
    return map
end

function _materialize_rules(map::CoefficientRuleMap)::Dict{Op, QAdd}
    _validate_map(map)
    rules = Dict{Op, QAdd}()
    sizehint!(rules, length(map.generators))
    @inbounds for i in eachindex(map.generators)
        terms = QTermDict()
        sizehint!(terms, length(map.coordinates) + 1)
        for j in eachindex(map.coordinates)
            coefficient = map.coefficients[i, j]
            SQA._iszero_cnum(coefficient) || SQA._addto_key!(
                terms, SQA._copy_key(map.coordinates[j]), coefficient,
            )
        end
        offset = map.offsets[i]
        SQA._iszero_cnum(offset) || SQA._addto!(terms, Op[], offset)
        rules[map.generators[i]] = QAdd(terms, SQA._EMPTY_INDICES)
    end
    return rules
end

function _materialize_gauge(operator_part::QAdd, scalar::CNum)::QAdd
    SQA._iszero_cnum(scalar) && return operator_part
    terms = copy(operator_part.arguments)
    SQA._addto!(terms, Op[], scalar)
    return QAdd(terms, SQA._EMPTY_INDICES)
end

function materialize_direct(
        pair::AffineRulePair,
        relations::Vector{ParamRelation} = ParamRelation[],
    )::SQA.UnitaryTransform{SQA.StaticTime}
    return SQA._static_transform(
        _materialize_rules(pair.forward), _materialize_rules(pair.inverse), relations,
    )
end

function materialize_direct(
        pair::AffineRulePair, operator_gauge::QAdd, scalar_gauge::CNum,
        t::SQA.Num, relations::Vector{ParamRelation} = ParamRelation[],
    )::SQA.UnitaryTransform{SQA.DynamicTime}
    static = materialize_direct(pair, relations)
    return SQA._timed_transform(
        static, _materialize_gauge(operator_gauge, scalar_gauge), t,
    )
end

function _rule_data(
        pair::AffineRulePair,
        relations::Vector{ParamRelation} = ParamRelation[],
    )::RuleTransformData
    rules = _materialize_rules(pair.forward)
    inverse_rules = _materialize_rules(pair.inverse)
    generators = sort!(collect(keys(rules)))
    return RuleTransformData(
        rules, inverse_rules, generators, SQA._site_infos(generators), relations,
    )
end

function materialize_via_wrapper(
        pair::AffineRulePair,
        relations::Vector{ParamRelation} = ParamRelation[],
    )::SQA.UnitaryTransform{SQA.StaticTime}
    data = _rule_data(pair, relations)
    return SQA._validated_transform(
        data.rules, data.inverse_rules, SQA._zero_qadd(), SQA.StaticTime(),
        data.relations,
    )
end

function nested(U::SQA.UnitaryTransform{T}) where {T}
    data = RuleTransformData(
        U.rules, U.inverse_rules, U.generators, U.sites, U.relations,
    )
    return NestedTransform{T}(data, U.gauge, U.time)
end

function _covered_site(U::NestedTransform, key)
    for site in U.data.sites
        site.key == key && return true
    end
    return false
end

function _validate_coverage(q::QAdd, U::NestedTransform)
    for (term, _) in q, operator in term.ops
        haskey(U.data.rules, operator) && continue
        _covered_site(U, SQA.site_key(operator)) || continue
        throw(ArgumentError("incomplete nested transform coverage"))
    end
    return nothing
end

function nested_conjugate(q::QAdd, U::NestedTransform)::QAdd
    _validate_coverage(q, U)
    return SQA._reduce_params(
        SQA._apply_rules(q, U.data.rules), U.data.relations, true,
    )
end

nested_transform(q::QAdd, U::NestedTransform{SQA.StaticTime}) =
    nested_conjugate(q, U)
nested_transform(q::QAdd, U::NestedTransform{SQA.DynamicTime}) =
    nested_conjugate(q, U) + U.gauge

function nested_inv(U::NestedTransform{T}) where {T}
    gauge = if T === SQA.StaticTime || iszero(U.gauge)
        U.gauge
    else
        -SQA._reduce_params(
            SQA._apply_rules(U.gauge, U.data.inverse_rules), U.data.relations, true,
        )
    end
    data = RuleTransformData(
        copy(U.data.inverse_rules), copy(U.data.rules), U.data.generators,
        U.data.sites, copy(U.data.relations),
    )
    return NestedTransform{T}(data, gauge, U.time)
end

function nested_compose(
        first::NestedTransform{T}, second::NestedTransform{T},
    ) where {T}
    relations = SQA._merge_relations(first.data.relations, second.data.relations)
    rules = SQA._compose_rules(first.data.rules, second.data.rules)
    inverse_rules = SQA._compose_rules(
        second.data.inverse_rules, first.data.inverse_rules,
    )
    gauge = if iszero(first.gauge)
        second.gauge
    else
        transported = SQA._reduce_params(
            SQA._apply_rules(first.gauge, second.data.rules), relations, true,
        )
        iszero(second.gauge) ? transported : SQA._add_gauges(transported, second.gauge)
    end
    generators = sort!(collect(keys(rules)))
    data = RuleTransformData(
        rules, inverse_rules, generators, SQA._site_infos(generators), relations,
    )
    return NestedTransform{T}(data, gauge, first.time)
end

function nested_render(U::NestedTransform)
    return sprint() do io
        write(io, "NestedTransform(")
        for (i, generator) in enumerate(U.data.generators)
            i > 1 && write(io, ", ")
            show(io, generator)
            write(io, " -> ")
            show(io, U.data.rules[generator])
        end
        write(io, ")")
    end
end

"""Convert an existing exact transform to closure-coordinate coefficient maps."""
function coefficient_pair(U::SQA.UnitaryTransform)::AffineRulePair
    coordinates = QTerm[]
    lookup = Dict{QTerm, Int}()
    for rules in (U.rules, U.inverse_rules), generator in U.generators
        for (term, _) in rules[generator]
            isempty(term.ops) && continue
            haskey(lookup, term) && continue
            push!(coordinates, SQA._copy_key(term))
            lookup[term] = length(coordinates)
        end
    end
    function scan(rules)
        coefficients = fill(SQA._CNUM_ZERO, length(U.generators), length(coordinates))
        offsets = fill(SQA._CNUM_ZERO, length(U.generators))
        for (i, generator) in enumerate(U.generators), (term, coefficient) in rules[generator]
            if isempty(term.ops)
                offsets[i] = coefficient
            else
                coefficients[i, lookup[term]] = coefficient
            end
        end
        return CoefficientRuleMap(U.generators, coordinates, coefficients, offsets)
    end
    return AffineRulePair(scan(U.rules), scan(U.inverse_rules))
end

function same_transform(left::SQA.UnitaryTransform, right::SQA.UnitaryTransform)
    left.generators == right.generators || return false
    left.rules == right.rules || return false
    left.inverse_rules == right.inverse_rules || return false
    left.gauge == right.gauge || return false
    return true
end

end

