@inline provably_real_replacement(x::Real) = true
@inline provably_real_replacement(x::Union{Num, SymbolicUtils.BasicSymbolic}) =
    SymbolicUtils.symtype(x) <: Real
@inline provably_real_replacement(x) = false

const ScalarRules = Dict{RawExpression, RawExpression}
const PhaseRules = Vector{Pair{RawExpression, RawExpression}}

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

@noinline function validate_phase_substitutions(c::Coeff, phase_rules::PhaseRules)
    isempty(phase_rules) && return c
    function check(x)
        is_phase(x) || return false
        argument = only(SymbolicUtils.arguments(x))
        for (from, to) in phase_rules
            raw_depends_on(argument, from) && nonreal_phase_argument(to)
        end
        return false
    end
    tail = c.tail
    if tail isa Poly
        foreach(check, (symbol for monomial in tail.terms for symbol in monomial.syms))
    elseif tail isa RawSymbolicCoeff
        has_subterm(check, tail.expr)
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
function SymbolicUtils.substitute(op::QSym, rules::AbstractDict; replace_adjoint::Bool = true)
    return substitute_with(to_qadd(op), SubstitutionRules(rules, replace_adjoint))
end

function SymbolicUtils.substitute(q::QAdd, rules::AbstractDict; replace_adjoint::Bool = true)
    return substitute_with(q, SubstitutionRules(rules, replace_adjoint))
end

const OpRules = Dict{Op, QAdd}

"""
    SubstitutionRules

The rules of one substitution: operator rules, scalar rules, and the scalar rules whose
replacement is not provably real, which must not reach the argument of a phase.
"""
struct SubstitutionRules
    ops::OpRules
    scalars::ScalarRules
    phases::PhaseRules
end
SubstitutionRules(ops::OpRules = OpRules(), scalars::ScalarRules = ScalarRules()) =
    SubstitutionRules(ops, scalars, nonreal_phase_substitutions(scalars))

function SubstitutionRules(rules::AbstractDict, replace_adjoint::Bool)
    op_rules = OpRules()
    scalar_rules = ScalarRules()
    for (k, v) in rules
        add_substitution_rule!(op_rules, scalar_rules, k, v)
    end
    replace_adjoint && add_adjoint_rules!(op_rules)
    return SubstitutionRules(op_rules, scalar_rules)
end

Base.isempty(r::SubstitutionRules) = isempty(r.ops) && isempty(r.scalars)

substitute_cnum(c::CNum, r::SubstitutionRules) = substitute_cnum(c, r.scalars, r.phases)

function substitute_with(q::QAdd, r::SubstitutionRules)
    (iszero(q) || isempty(r)) && return q
    out = QTermDict()
    indices = copy(q.indices)
    for (t, c) in q
        temp = QTermDict()
        replacement_indices = substitute_term!(temp, t.ops, c, t.ne, r)
        settle_scope!(out, temp, bound_indices(t, c, q.indices))
        indices = merge_unique(indices, replacement_indices)
    end
    return QAdd(out, drop_unused_indices(out, indices))
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

function substitute_term!(
        out::QTermDict, ops::Vector{Op}, c::CNum, ne::Vector{NonEqualPair}, r::SubstitutionRules,
    )
    if !any(op -> haskey(r.ops, op), ops)
        new_c = substitute_cnum(c, r)
        iszero(new_c) || stream!(out, copy(ops), new_c, ne)
        return EMPTY_INDICES
    end
    indices = EMPTY_INDICES
    for op in ops
        indices = merge_unique(indices, get(r.ops, op, ZERO_QADD).indices)
    end
    substitute_ops(Op[], c, ne, ops, 1, r.ops) do new_ops, new_c, new_ne
        new_c = substitute_cnum(new_c, r)
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
replacement_qadd(::Op, v::Op) = to_qadd(v)
replacement_qadd(::Op, v::Union{Number, SymbolicUtils.BasicSymbolic, Coeff}) = single_qadd(to_cnum(v), Op[])
replacement_qadd(op::Op, v) =
    throw(ArgumentError("operator replacement for `$op` has unsupported value `$v`"))
