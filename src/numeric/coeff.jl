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
    return op_rules, scalar_rules
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
        if isequal(uk, conj(vs[1]))
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

const_coeff(c::Complex{Num}) = to_complex(c)

as_cnum(x::Complex) = Complex{Num}(Num(real(x)), Num(imag(x)))
as_cnum(x) = Complex{Num}(Num(x), Num(false))

function raw_variables(x::RawExpression)
    vars = RawExpression[]
    for v in Symbolics.get_variables(x)
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

function compile_td_coeffs(q::QAdd, time_parameter)
    basevars, valuefuncs, p_aware = time_basis(time_parameter)
    exprs = Complex{Num}[to_num(c) for (_, c) in q.arguments]
    isempty(exprs) && push!(exprs, Complex{Num}(0))
    for e in exprs
        check_time_variables(e, basevars)
    end
    return compile_coeffs(exprs, basevars, valuefuncs, p_aware)
end

function compile_coeffs(exprs::Vector{Complex{Num}}, basevars::Vector{Num}, @nospecialize(valuefuncs), p_aware::Bool)
    f_re = first(build_function(Num[real(e) for e in exprs], basevars...; expression = Val(false)))
    f_im = first(build_function(Num[imag(e) for e in exprs], basevars...; expression = Val(false)))
    return let f_re = f_re, f_im = f_im, valuefuncs = valuefuncs
        if p_aware
            (u, p, t) -> complex_values(f_re, f_im, map(f -> f(p, t), valuefuncs))
        else
            (u, p, t) -> complex_values(f_re, f_im, map(f -> f(t), valuefuncs))
        end
    end
end

function complex_values(f_re::F, f_im::G, vals::Tuple) where {F, G}
    re = f_re(vals...)
    imag_part = f_im(vals...)
    return ComplexF64[re[k] + im * imag_part[k] for k in eachindex(re)]
end

function split_scalar_subs(scalar_subs::AbstractDict)
    sub_re = ScalarRules()
    sub_im = ScalarRules()
    has_imag = false
    for (k, v) in scalar_subs
        has_imag |= add_split_scalar_sub!(sub_re, sub_im, raw_rule(k), v)
    end
    return sub_re, sub_im, has_imag
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

function apply_scalar_subs(c::CNum, sub_re::ScalarRules, sub_im::ScalarRules, has_imag::Bool)
    (isempty(sub_re) || is_native(c)) && return c
    re, im = realimag(c)
    re_part = SymbolicUtils.unwrap(re)
    im_part = SymbolicUtils.unwrap(im)
    if !has_imag
        return cnum(
            Num(Symbolics.substitute(re_part, sub_re)),
            Num(Symbolics.substitute(im_part, sub_re)),
        )
    end
    new_re = Symbolics.substitute(re_part, sub_re) -
        Symbolics.substitute(im_part, sub_im)
    new_im = Symbolics.substitute(re_part, sub_im) +
        Symbolics.substitute(im_part, sub_re)
    return cnum(Num(new_re), Num(new_im))
end

function reduce_const(n::Num)
    x = SymbolicUtils.unwrap(n)
    try
        return fold_const(x)
    catch err
        err isa ArgumentError || rethrow()
        isempty(Symbolics.get_variables(n)) || rethrow()
        return compile_const(n)
    end
end

function compile_const(n::Num)
    x = SymbolicUtils.evaluate(SymbolicUtils.unwrap(n))
    is_float_const(x) && return const_complex(x)
    throw(ArgumentError("cannot reduce symbolic expression $x to a concrete number"))
end

function fold_const(x::RawExpression)
    if SymbolicUtils.iscall(x)
        op = SymbolicUtils.operation(x)
        args = SymbolicUtils.arguments(x)
        if op === (+)
            acc = zero(ComplexF64)
            for a in args
                acc += fold_const(a)
            end
            return acc
        elseif op === (*)
            acc = one(ComplexF64)
            for a in args
                acc *= fold_const(a)
            end
            return acc
        elseif op === conj
            return conj(fold_const(first(args)))
        elseif op === real
            return real(fold_const(first(args)))
        elseif op === imag
            return imag(fold_const(first(args)))
        elseif op === (/)
            return fold_const(first(args)) / fold_const(last(args))
        elseif op === (^)
            return fold_const(first(args))^fold_const(last(args))
        elseif op === (-)
            return length(args) == 1 ? -fold_const(first(args)) :
                fold_const(first(args)) - fold_const(last(args))
        elseif op === expim
            return exp(im * fold_const(only(args)))
        elseif op === sqrt
            return sqrt(fold_const(only(args)))
        elseif op === cbrt
            z = fold_const(only(args))
            iszero(imag(z)) && return ComplexF64(cbrt(real(z)))
        end
    elseif is_float_const(x)
        return const_complex(x)
    end
    throw(ArgumentError("cannot reduce symbolic expression $x to a concrete number"))
end

to_complex(c::Coeff) = is_native(c) ? native_float(native_scalar(c)) : to_complex(to_num(c))

function to_complex(x)
    x isa ComplexF64   && return x
    x isa Complex{Num} && return reduce_const(real(x)) + im * reduce_const(imag(x))
    x isa Complex      && return convert(ComplexF64, x)
    x isa Num          && return reduce_const(x)
    if x isa SymbolicUtils.BasicSymbolic
        (x isa RawExpression && is_float_const(x)) ||
            throw(ArgumentError("cannot reduce symbolic expression $x to a concrete number"))
        return const_complex(x)
    end
    return convert(ComplexF64, Complex(x, false))
end
