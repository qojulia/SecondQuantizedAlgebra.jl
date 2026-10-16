@inline provably_real_replacement(x::Real) = true
@inline provably_real_replacement(x::Num) = SymbolicUtils.symtype(x) <: Real
@inline provably_real_replacement(x::SymbolicUtils.BasicSymbolic) =
    SymbolicUtils.symtype(x) <: Real
@inline provably_real_replacement(x) = false


const ScalarRules = Dict{RawExpression, RawExpression}
const PhaseRules = Vector{Pair{RawExpression, RawExpression}}

const EMPTY_SCALAR_RULES = ScalarRules()
const EMPTY_PHASE_SUBSTITUTIONS = PhaseRules()

raw_rule(x::RawExpression) = x
raw_rule(x::Num) = raw_rule(SymbolicUtils.unwrap(x))
raw_rule(x::Coeff) = lower_complex_slots(raw_tail(x))
raw_rule(x::Number) = SymbolicUtils.Const{SymbolicUtils.SymReal}(x)
raw_rule(x) = throw(ArgumentError("unsupported scalar substitution `$x` of type $(typeof(x))"))

scalar_substitution_rules(d::ScalarRules) = d
function scalar_substitution_rules(d::AbstractDict)
    out = ScalarRules()
    for (from, to) in d
        add_scalar_rule!(out, from, to)
    end
    return out
end

function add_scalar_rule!(rules::ScalarRules, from::Union{Number, RawExpression}, to)
    rules[raw_rule(from)] = raw_rule(to)
    return rules
end
add_scalar_rule!(rules::ScalarRules, from, to) = rules

function nonreal_phase_substitutions(d::ScalarRules)
    isempty(d) && return EMPTY_PHASE_SUBSTITUTIONS
    rules = PhaseRules()
    for (from, to) in d
        provably_real_replacement(to) || push!(rules, from => to)
    end
    return rules
end

function collect_phase_arguments!(out::Vector{RawExpression}, x::RawExpression)
    if is_phase(x)
        push!(out, only(SymbolicUtils.arguments(x)))
        return out
    end
    SymbolicUtils.iscall(x) || return out
    for argument in SymbolicUtils.arguments(x)
        collect_phase_arguments!(out, argument)
    end
    return out
end

@noinline function validate_phase_substitutions(c::Coeff, phase_rules::PhaseRules)
    isempty(phase_rules) && return c
    tail = c.tail
    phase_arguments = RawExpression[]
    if tail isa Poly
        for monomial in tail.terms, symbol in monomial.syms
            is_phase(symbol) && push!(phase_arguments, only(SymbolicUtils.arguments(symbol)))
        end
    elseif tail isa RawSymbolicCoeff
        collect_phase_arguments!(phase_arguments, tail.expr)
    end
    for argument in phase_arguments
        for (from, to) in phase_rules
            raw_depends_on(argument, from) && nonreal_phase_argument(to)
        end
    end
    return c
end

function substitute_cnum(c::CNum, d::ScalarRules, phase_rules::PhaseRules = nonreal_phase_substitutions(d))
    isempty(d) && return c
    is_native(c) && return c
    validate_phase_substitutions(c, phase_rules)
    raw = raw_rule(c)
    new_raw = SymbolicUtils.substitute(raw, d)
    isequal(new_raw, raw) && return c
    return to_cnum(new_raw)
end
substitute_cnum(c::CNum, d::AbstractDict) = substitute_cnum(c, scalar_substitution_rules(d))

SymbolicUtils.substitute(c::Coeff, rules::AbstractDict) = substitute_cnum(c, rules)

"""
    substitute(expr::QField, rules::AbstractDict; replace_adjoint=true)

Replace parameters and operators in `expr` according to `rules`.

Keys may be parameters or operators, mixed freely. Operator rules act on the
original expression all at once, so a rule whose replacement contains its own
key, such as the shift `a => a + α`, is applied exactly once. A rule for `a` also
covers `a'`, unless `a'` has a rule of its own or `replace_adjoint = false`.

# Examples

```jldoctest
julia> a = Destroy(FockSpace(:cavity), :a);

julia> @variables Δ α::Number;

julia> substitute(Δ * a' * a, Dict(Δ => 2))
2 * a' * a

julia> substitute(a' * a, Dict(a => a + α))
conj(α)*α + conj(α) * a + α * a' + a' * a
```

See also [`change_index`](@ref), [`Displace`](@ref).
"""
function SymbolicUtils.substitute(op::QSym, rules::AbstractDict; replace_adjoint = true)
    return SymbolicUtils.substitute(
        single_qadd(CNUM_ONE, Op[op]), rules; replace_adjoint = replace_adjoint
    )
end

function SymbolicUtils.substitute(q::QAdd, rules::AbstractDict; replace_adjoint::Bool = true)
    iszero(q) && return q
    op_rules, scalar_rules = split_substitution_rules(rules, replace_adjoint)
    return substitute_split(q, op_rules, scalar_rules)
end

const OpRules = Dict{Op, QAdd}

const EMPTY_OP_RULES = OpRules()

substitute_op_rules(q::QAdd, op_rules::OpRules) =
    iszero(q) ? q : substitute_split(q, op_rules, EMPTY_SCALAR_RULES)

function substitute_split(q::QAdd, op_rules::OpRules, scalar_rules::ScalarRules)
    out = QTermDict()
    indices = copy(q.indices)
    phase_rules = nonreal_phase_substitutions(scalar_rules)
    for (t, c) in q
        temp = QTermDict()
        replacement_indices = stream_substitution_once!(
            temp, t.ops, c, t.ne, op_rules, scalar_rules, phase_rules,
        )
        settle_scope!(out, temp, bound_indices(t, c, q.indices))
        indices = merge_unique(indices, replacement_indices)
    end
    return QAdd(out, drop_unused_indices(out, indices))
end

function split_substitution_rules(rules::AbstractDict, replace_adjoint::Bool)
    op_rules = OpRules()
    scalar_rules = ScalarRules()
    for (k, v) in rules
        add_substitution_rule!(op_rules, scalar_rules, k, v)
    end
    replace_adjoint && add_adjoint_rules!(op_rules)
    return op_rules, scalar_rules
end

function add_substitution_rule!(op_rules::OpRules, scalar_rules::ScalarRules, k::Op, v)
    op_rules[k] = replacement_qadd(k, v)
    return op_rules
end
function add_substitution_rule!(op_rules::OpRules, scalar_rules::ScalarRules, k, v)
    add_scalar_rule!(scalar_rules, k, v)
    return op_rules
end

function add_adjoint_rules!(op_rules::OpRules)
    for (k, v) in collect(op_rules)
        k_adj = Base.adjoint(k)
        haskey(op_rules, k_adj) || (op_rules[k_adj] = qadjoint(v))
    end
    return op_rules
end

function has_operator_rule(ops::Vector{Op}, dict::OpRules)
    for op in ops
        haskey(dict, op) && return true
    end
    return false
end

function stream_substitution_once!(
        out::QTermDict, ops::Vector{Op}, c::CNum, ne::Vector{NonEqualPair},
        op_rules::OpRules, scalar_rules::ScalarRules, phase_rules::PhaseRules,
    )
    if !has_operator_rule(ops, op_rules)
        new_c = substitute_cnum(c, scalar_rules, phase_rules)
        iszero(new_c) || stream!(out, copy(ops), new_c, ne)
        return EMPTY_INDICES
    end
    indices = EMPTY_INDICES
    for op in ops
        haskey(op_rules, op) && (indices = merge_unique(indices, op_rules[op].indices))
    end
    substitute_ops(Op[], c, ne, ops, 1, op_rules) do new_ops, new_c, new_ne
        new_c = substitute_cnum(new_c, scalar_rules, phase_rules)
        iszero(new_c) || canonicalize!(out, new_ops, new_c, new_ne)
    end
    return indices
end

function substitute_ops(
        sink::F, prefix::Vector{Op}, c::CNum, ne::Vector{NonEqualPair},
        ops::Vector{Op}, k::Int, op_rules::OpRules,
    ) where {F}
    k > length(ops) && (sink(prefix, c, ne); return)
    if !haskey(op_rules, ops[k])
        push!(prefix, ops[k])
        substitute_ops(sink, prefix, c, ne, ops, k + 1, op_rules)
        return
    end
    for (rt, rc) in op_rules[ops[k]]
        branch = append!(copy(prefix), rt.ops)
        substitute_ops(sink, branch, c * rc, merge_ne(ne, rt.ne), ops, k + 1, op_rules)
    end
    return
end

replacement_qadd(::Op, v::QAdd) = v
replacement_qadd(::Op, v::Op) = single_qadd(CNUM_ONE, Op[v])
replacement_qadd(::Op, v::Union{Number, SymbolicUtils.BasicSymbolic, Coeff}) = single_qadd(to_cnum(v), Op[])
replacement_qadd(op::Op, v) =
    throw(ArgumentError("operator replacement for `$op` has unsupported value `$v`"))
