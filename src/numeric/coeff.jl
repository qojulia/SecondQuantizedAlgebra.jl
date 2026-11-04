const NO_SCALAR_SUBS = ScalarRules()
const NO_TIME_PARAMETER = ()
const NO_PARAMETER = ScalarRules()

function expand_parameter(parameter::AbstractDict)
    op_rules = OpRules()
    scalar_rules = ScalarRules()
    for (k, v) in parameter
        add_parameter_rule!(op_rules, scalar_rules, k, v)
    end
    add_adjoint_rules!(op_rules)
    return SubstitutionRules(op_rules, scalar_rules)
end

function add_parameter_rule!(op_rules::OpRules, scalar_rules::ScalarRules, k::Complex, v)
    scalar_rules[raw_rule(real(k))] = raw_rule(real(v))
    ik = raw_rule(imag(k))
    SymbolicUtils.isconst(ik) || (scalar_rules[ik] = raw_rule(imag(v)))
    return op_rules
end
add_parameter_rule!(op_rules::OpRules, scalar_rules::ScalarRules, k, v) =
    add_substitution_rule!(op_rules, scalar_rules, k, v)

normalize_time_parameter(::Tuple{}) = NO_TIME_PARAMETER
function normalize_time_parameter(time_parameter::AbstractDict)
    out = Dict{Num, Any}()
    for (k, v) in time_parameter
        out[time_key(k)] = time_value(v)
    end
    return out
end

time_key(k::Num) = k
time_key(k::RawExpression) = Num(k)
function time_key(k::Complex{Num})
    u = SymbolicUtils.unwrap(k)
    SymbolicUtils.issym(u) && return Num(u)
    v = SymbolicUtils.unwrap(conj(k))
    SymbolicUtils.issym(v) && return Num(conj(v))
    throw(ArgumentError("time_parameter key `$k` must be a symbolic variable or its conjugate"))
end
time_key(k) = throw(ArgumentError("time_parameter key `$k` must be a symbolic variable"))

time_value(v::Number) = t -> v
time_value(v) = v

function tp_reads_p(@nospecialize(f))
    arities = Set{Int}()
    for m in methods(f)
        m.isva && throw(
            ArgumentError(
                "time_parameter value `$f` is variadic; write it as an explicit " *
                    "`t -> value` or `(p, t) -> value` closure.",
            ),
        )
        push!(arities, Int(m.nargs) - 1)
    end
    length(arities) == 1 || throw(
        ArgumentError(
            "time_parameter value `$f` has methods of conflicting arity $(sort!(collect(arities))); " *
                "write it as an explicit `t -> value` or `(p, t) -> value` closure.",
        ),
    )
    n = only(arities)
    n == 1 && return false
    n == 2 && return true
    throw(
        ArgumentError(
            "time_parameter value `$f` must take one argument (`t -> value`) or two arguments " *
                "(`(p, t) -> value`), but takes $n.",
        ),
    )
end

function lift_pt(@nospecialize(f))
    tp_reads_p(f) && return f
    return let f = f
        (p, t) -> f(t)
    end
end

function conj_vf(@nospecialize(f), p_aware::Bool)
    return let f = f
        p_aware ? ((p, t) -> conj(f(p, t))) : (t -> conj(f(t)))
    end
end

function time_basis(time_parameter::Dict{Num, Any})
    p_aware = any(tp_reads_p, values(time_parameter))
    basevars = Num[]
    valuefuncs = Any[]
    for (k, f) in time_parameter
        g = p_aware ? lift_pt(f) : f
        uk = SymbolicUtils.unwrap(k)
        if SymbolicUtils.issym(uk)
            push!(basevars, k)
            push!(valuefuncs, g)
            continue
        end
        vs = raw_variables(uk)
        length(vs) == 1 || throw(
            ArgumentError(
                "time_parameter key `$k` must depend on exactly one variable, got $(length(vs)).",
            ),
        )
        if isequal(uk, vs[1])
            push!(basevars, k)
            push!(valuefuncs, g)
        elseif isequal(uk, conj(vs[1]))
            push!(basevars, Num(vs[1]))
            push!(valuefuncs, conj_vf(g, p_aware))
        else
            throw(
                ArgumentError(
                    "unsupported time_parameter key `$k`; only a bare variable `v` or `conj(v)` are supported.",
                ),
            )
        end
    end
    return basevars, Tuple(valuefuncs), p_aware
end

coeff_is_const(c::Complex{Num}) =
    isempty(Symbolics.get_variables(real(c))) && isempty(Symbolics.get_variables(imag(c)))

as_cnum(x::Complex) = Complex{Num}(Num(real(x)), Num(imag(x)))
as_cnum(x) = Complex{Num}(Num(x), Num(false))

# IndexedVariable calls lack Symbolics' dependent-variable metadata.
function coefficient_atom(x::SymbolicUtils.BasicSymbolic)
    SymbolicUtils.default_is_atomic(x) && return true
    SymbolicUtils.iscall(x) || return false
    op = SymbolicUtils.operation(x)
    return op isa SymbolicUtils.BasicSymbolic && SymbolicUtils.issym(op)
end

function raw_variables(x::RawExpression)
    vars = RawExpression[]
    for v in Symbolics.get_variables(x; is_atomic = coefficient_atom)
        push!(vars, v)
    end
    return vars
end

function coefficient_variables(c::Complex{Num})
    vars = raw_variables(SymbolicUtils.unwrap(real(c)))
    append!(vars, raw_variables(SymbolicUtils.unwrap(imag(c))))
    unique!(vars)
    return vars
end

function check_time_variables(c::Complex{Num}, basevars::Vector{Num})
    vars = coefficient_variables(c)
    missing = RawExpression[]
    for v in vars
        any(b -> isequal(v, SymbolicUtils.unwrap(b)), basevars) || push!(missing, v)
    end
    isempty(missing) && return c
    throw(
        ArgumentError(
            "time-dependent coefficient `$c` depends on variables without time values: " *
                join(map(string, missing), ", "),
        ),
    )
end

function compile_coeff(c::Complex{Num}, vars...)
    f_re = build_function(real(c), vars...; expression = Val(false))
    f_im = build_function(imag(c), vars...; expression = Val(false))
    g_re = f_re isa Tuple ? first(f_re) : f_re
    g_im = f_im isa Tuple ? first(f_im) : f_im
    return (vals...) -> g_re(vals...) + im * g_im(vals...)
end

function at_time(g::G, @nospecialize(valuefuncs), p_aware::Bool) where {G}
    return let valuefuncs = valuefuncs
        p_aware ? (p, t) -> g(t, map(f -> f(p, t), valuefuncs)) : t -> g(t, map(f -> f(t), valuefuncs))
    end
end

const CoeffTerm = Tuple{CNum, Vector{Op}}

coeff_terms(q::QAdd) = CoeffTerm[(c, t.ops) for (t, c) in q.arguments]

struct SplitScalarRules
    re::ScalarRules
    im::ScalarRules
    has_imag::Bool
end

function split_scalar_subs(scalar_subs::AbstractDict)
    sub_re = ScalarRules()
    sub_im = ScalarRules()
    has_imag = false
    for (k, v) in scalar_subs
        has_imag |= add_split_scalar_sub!(sub_re, sub_im, raw_rule(k), v)
    end
    return SplitScalarRules(sub_re, sub_im, has_imag)
end

function add_split_scalar_sub!(sub_re::ScalarRules, sub_im::ScalarRules, k::RawExpression, v::Complex)
    sub_re[k] = raw_rule(real(v))
    sub_im[k] = raw_rule(imag(v))
    return !iszero(imag(v))
end
function add_split_scalar_sub!(sub_re::ScalarRules, sub_im::ScalarRules, k::RawExpression, v)
    sub_re[k] = raw_rule(v)
    sub_im[k] = RAW_ZERO
    return false
end

function apply_scalar_subs(c::CNum, s::SplitScalarRules)
    (isempty(s.re) || is_native(c)) && return c
    re, im = realimag(c)
    re_part = SymbolicUtils.unwrap(re)
    im_part = SymbolicUtils.unwrap(im)
    if !s.has_imag
        return cnum(Num(Symbolics.substitute(re_part, s.re)), Num(Symbolics.substitute(im_part, s.re)))
    end
    new_re = Symbolics.substitute(re_part, s.re) - Symbolics.substitute(im_part, s.im)
    new_im = Symbolics.substitute(re_part, s.im) + Symbolics.substitute(im_part, s.re)
    return cnum(Num(new_re), Num(new_im))
end

function reduce_const(n::Num)
    x = SymbolicUtils.unwrap(n)
    try
        return fold_const(x)
    catch err
        err isa ArgumentError || rethrow()
        isempty(Symbolics.get_variables(n)) || rethrow()
        v = SymbolicUtils.evaluate(x)
        is_float_const(v) && return const_complex(v)
        return not_const_error(v)
    end
end

function fold_const(x::RawExpression)
    SymbolicUtils.iscall(x) || return is_float_const(x) ? const_complex(x) : not_const_error(x)
    op = SymbolicUtils.operation(x)
    args = SymbolicUtils.arguments(x)
    if op === (+)
        acc = zero(ComplexF64)
        for a in args
            acc += fold_const(a)
        end
        return acc
    end
    if op === (*)
        acc = one(ComplexF64)
        for a in args
            acc *= fold_const(a)
        end
        return acc
    end
    op === conj && return conj(fold_const(first(args)))
    op === real && return real(fold_const(first(args)))
    op === imag && return imag(fold_const(first(args)))
    op === (/) && return fold_const(first(args)) / fold_const(last(args))
    op === (^) && return fold_const(first(args))^fold_const(last(args))
    if op === (-)
        length(args) == 1 && return -fold_const(first(args))
        return fold_const(first(args)) - fold_const(last(args))
    end
    op === expim && return exp(im * fold_const(only(args)))
    op === sqrt && return sqrt(fold_const(only(args)))
    if op === cbrt
        z = fold_const(only(args))
        iszero(imag(z)) && return ComplexF64(cbrt(real(z)))
    end
    return not_const_error(x)
end

@noinline not_const_error(x) = throw(ArgumentError("cannot reduce symbolic expression $x to a concrete number"))

to_complex(c::Coeff) = is_native(c) ? to_float_scalar(native_scalar(c)) : to_complex(to_num(c))
function to_complex(x)
    x isa ComplexF64 && return x
    x isa Complex{Num} && return reduce_const(real(x)) + im * reduce_const(imag(x))
    x isa Complex && return convert(ComplexF64, x)
    x isa Num && return reduce_const(x)
    x isa RawExpression && is_float_const(x) && return const_complex(x)
    x isa SymbolicUtils.BasicSymbolic && return not_const_error(x)
    return convert(ComplexF64, Complex(x, false))
end
