# Coefficient lowering: reduce a (possibly symbolic) `CNum`/`Num` prefactor to a concrete
# `ComplexF64`, split user parameter/time substitutions, and compile time-dependent
# coefficient functions. Entirely backend-independent.

const _NO_SCALAR_SUBS = Dict{Num, Any}()

# --- parameter / time-parameter normalisation -----------------------------------------

function _expand_parameter(parameter)
    isempty(parameter) && return parameter
    out = Dict{Any, Any}()
    for (k, v) in parameter
        if k isa Complex
            out[real(k)] = real(v)
            ik = imag(k)
            SymbolicUtils.unwrap(ik) isa Number || (out[ik] = imag(v))
        else
            out[k] = v
        end
    end
    return out
end

function _normalize_time_parameter(time_parameter)
    isempty(time_parameter) && return time_parameter
    out = Dict{Any, Any}()
    for (k, v) in time_parameter
        out[k] = v isa Number ? (t -> v) : v
    end
    return out
end

function _time_basis(time_parameter)
    basevars = Any[]
    valuefuncs = Any[]
    for (k, f) in time_parameter
        uk = SymbolicUtils.unwrap(k)
        if SymbolicUtils.issym(uk)
            push!(basevars, k)
            push!(valuefuncs, f)
            continue
        end
        vs = collect(Symbolics.get_variables(k))
        length(vs) == 1 || throw(
            ArgumentError(
                "time_parameter key `$k` must depend on exactly one variable, got $(length(vs)).",
            ),
        )
        wv = Symbolics.wrap(vs[1])
        if isequal(uk, SymbolicUtils.unwrap(conj(wv)))
            push!(basevars, wv)
            push!(valuefuncs, t -> conj(f(t)))
        else
            throw(
                ArgumentError(
                    "unsupported time_parameter key `$k`; only a bare variable `v` or `conj(v)` are supported.",
                ),
            )
        end
    end
    return basevars, Tuple(valuefuncs)
end

_coeff_is_const(c::Complex{Num}) =
    isempty(Symbolics.get_variables(real(c))) && isempty(Symbolics.get_variables(imag(c)))

_const_coeff(c::Complex{Num}) = _to_complex(c)

_as_cnum(x::Complex) = Complex{Num}(Num(real(x)), Num(imag(x)))
_as_cnum(x) = Complex{Num}(Num(x), Num(false))

function _coefficient_variables(c::Complex{Num})
    vars = Any[]
    append!(vars, Symbolics.get_variables(real(c)))
    append!(vars, Symbolics.get_variables(imag(c)))
    unique!(vars)
    return vars
end

function _check_time_variables(c::Complex{Num}, basevars)
    vars = _coefficient_variables(c)
    base_unwrapped = SymbolicUtils.unwrap.(basevars)
    missing = Any[]
    for v in vars
        any(b -> isequal(v, b), base_unwrapped) || push!(missing, v)
    end
    isempty(missing) && return nothing
    throw(
        ArgumentError(
            "time-dependent coefficient `$c` depends on variables without time values: " *
                join(string.(missing), ", "),
        ),
    )
end

function _compile_coeff(c::Complex{Num}, vars...)
    f_re = build_function(real(c), vars...; expression = Val(false))
    f_im = build_function(imag(c), vars...; expression = Val(false))
    g_re = f_re isa Tuple ? first(f_re) : f_re
    g_im = f_im isa Tuple ? first(f_im) : f_im
    return (vals...) -> g_re(vals...) + im * g_im(vals...)
end

# --- indexed scalar substitutions (used by the indexed unroll) ------------------------

# Split user-supplied scalar substitutions into real and imag legs once per
# `to_numeric` call. A complex RHS `re + im*ii` propagates as separate real
# and imaginary substitutions, preserving the `Complex{Num}` invariant.
function _split_scalar_subs(scalar_subs::AbstractDict)
    sub_re = Dict{Any, Any}()
    sub_im = Dict{Any, Any}()
    has_imag = false
    for (k, v) in scalar_subs
        kraw = SymbolicUtils.unwrap(k)
        if v isa Complex
            sub_re[kraw] = real(v)
            sub_im[kraw] = imag(v)
            iszero(imag(v)) || (has_imag = true)
        else
            sub_re[kraw] = v
            sub_im[kraw] = 0
        end
    end
    return sub_re, sub_im, has_imag
end

function _apply_scalar_subs(c::CNum, sub_re::Dict, sub_im::Dict, has_imag::Bool)
    (isempty(sub_re) || _is_native(c)) && return c
    re_part = SymbolicUtils.unwrap(real(c))
    im_part = SymbolicUtils.unwrap(imag(c))
    if !has_imag
        return _cnum(
            Num(Symbolics.substitute(re_part, sub_re)),
            Num(Symbolics.substitute(im_part, sub_re)),
        )
    end
    # (re + i*im) -> (re|sub_re - im|sub_im) + i*(re|sub_im + im|sub_re)
    new_re = Symbolics.substitute(re_part, sub_re) -
        Symbolics.substitute(im_part, sub_im)
    new_im = Symbolics.substitute(re_part, sub_im) +
        Symbolics.substitute(im_part, sub_re)
    return _cnum(Num(new_re), Num(new_im))
end

# --- concrete reduction ----------------------------------------------------------------

function _reduce_const(n::Num)::ComplexF64
    x = Symbolics.value(n)
    try
        return _fold_const(x)
    catch err
        err isa ArgumentError || rethrow()
        isempty(Symbolics.get_variables(n)) || rethrow()
        return _compile_const(n)
    end
end

_compile_const(n::Num)::ComplexF64 = ComplexF64(symbolic_to_float(n))

function _fold_const(x)::ComplexF64
    x isa Number && return x
    if x isa SymbolicUtils.BasicSymbolic
        if SymbolicUtils.iscall(x)
            op = SymbolicUtils.operation(x)
            args = SymbolicUtils.arguments(x)
            if op === (+)
                acc = zero(ComplexF64)
                for a in args
                    acc += _fold_const(a)
                end
                return acc
            elseif op === (*)
                acc = one(ComplexF64)
                for a in args
                    acc *= _fold_const(a)
                end
                return acc
            elseif op === conj
                return conj(_fold_const(first(args)))
            elseif op === real
                return real(_fold_const(first(args)))
            elseif op === imag
                return imag(_fold_const(first(args)))
            elseif op === (/)
                return _fold_const(first(args)) / _fold_const(last(args))
            elseif op === (^)
                return _fold_const(first(args))^_fold_const(last(args))
            elseif op === (-)
                return length(args) == 1 ? -_fold_const(first(args)) :
                    _fold_const(first(args)) - _fold_const(last(args))
            end
        elseif SymbolicUtils.isconst(x)
            return x.val
        end
    end
    throw(ArgumentError("cannot reduce symbolic expression $x to a concrete number"))
end

_to_complex(c::Coeff) = _is_native(c) ? c.z : _to_complex(to_num(c))

# One method (union-split budget) routing every input through `convert ∘ Complex`,
# the only pattern that infers to ComplexF64 from `Any` after `Symbolics.value`.
function _to_complex(x)
    x isa ComplexF64   && return x
    x isa Complex{Num} && return _reduce_const(real(x)) + im * _reduce_const(imag(x))
    x isa Complex      && return convert(ComplexF64, x)
    x isa Num          && return _reduce_const(x)
    if x isa SymbolicUtils.BasicSymbolic
        SymbolicUtils.isconst(x) || throw(ArgumentError("cannot reduce symbolic expression $x to a concrete number"))
        return convert(ComplexF64, Complex(x.val, false))
    end
    return convert(ComplexF64, Complex(x, false))
end
