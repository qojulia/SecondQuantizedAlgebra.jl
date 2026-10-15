@inline is_one_at_zero(op) = op === exp || op === cos || op === cosh
@inline is_zero_at_zero(op) = op === sin || op === tan || op === sinh || op === tanh

const RAW_ZERO = SymbolicUtils.Const{SymbolicUtils.SymReal}(0)
const RAW_ONE = SymbolicUtils.Const{SymbolicUtils.SymReal}(1)

const ConstInt = Union{Int, Bool}
const ConstRational = Union{Int, BigInt, Bool, Rational{Int}, Rational{BigInt}}
const ConstReal = Union{ConstRational, Float64, BigFloat}

@inline is_const_zero(x::RawExpression) = SymbolicUtils.isconst(x) && SymbolicUtils._iszero(x)
@inline const_payload_isa(x::RawExpression, ::Type{T}) where {T} =
    SymbolicUtils.isconst(x) && x.val isa T

const ConstNumber = Union{
    ConstReal, Complex{Int}, ComplexF64, Complex{Bool}, Complex{Rational{Int}},
}
const ConstIrrational = Union{
    Irrational{:π}, Irrational{:ℯ}, Irrational{:γ}, Irrational{:φ}, Irrational{:catalan},
}

"""
    with_const(f, x)

`f` applied to the payload of the constant `x`, with the payload's concrete type, so that
every reader splits the untyped `x.val` the same way.
"""
function with_const(f::F, x::SymbolicUtils.BasicSymbolic) where {F}
    v = x.val
    v isa Int && return f(v)
    v isa Float64 && return f(v)
    v isa Rational{Int} && return f(v)
    v isa BigInt && return f(v)
    v isa Rational{BigInt} && return f(v)
    v isa BigFloat && return f(v)
    v isa Bool && return f(v)
    v isa Complex{Int} && return f(v)
    v isa ComplexF64 && return f(v)
    v isa Complex{Bool} && return f(v)
    v isa Complex{Rational{Int}} && return f(v)
    v isa Irrational{:π} && return f(v)
    v isa Irrational{:ℯ} && return f(v)
    v isa Irrational{:γ} && return f(v)
    v isa Irrational{:φ} && return f(v)
    v isa Irrational{:catalan} && return f(v)
    throw(ArgumentError("not a literal constant: $x"))
end

@inline is_const_int(x::RawExpression) = const_payload_isa(x, ConstInt)
const_int(x::RawExpression) = Int(x.val::ConstInt)

@inline is_rational_exponent(x::RawExpression) = const_payload_isa(x, Rational{Int})
rational_exponent(x::RawExpression) = x.val::Rational{Int}

@inline is_exact_number(x::RawExpression) = const_payload_isa(x, ConstRational)
exact_number(x::RawExpression) = with_const(Rational{BigInt}, x)

@inline is_real_const(x::RawExpression) = const_payload_isa(x, ConstReal)
real_const_float(x::RawExpression) = with_const(Float64, x)

function const_sign(x::RawExpression)
    is_real_const(x) || return 0
    v = real_const_float(x)
    return isnan(v) ? 0 : Int(sign(v))
end

is_const_minus_one(x::SymbolicUtils.BasicSymbolic) =
    SymbolicUtils.isconst(x) && x.val isa ConstReal && with_const(==(-1), x)

@inline is_numeric_const(x::RawExpression) = const_payload_isa(x, ConstNumber)
const_coeff(x::RawExpression) = with_const(to_cnum, x)

@inline is_float_const(x::RawExpression) = const_payload_isa(x, Union{ConstNumber, ConstIrrational})
const_complex(x::RawExpression) = with_const(ComplexF64, x)

"""
    raw_coeff(x, known_real)

The coefficient of the raw expression `x`: a number when `x` is a numeric constant, the raw
expression otherwise.
"""
@inline function raw_coeff(x::RawExpression, known_real::Bool)
    is_numeric_const(x) && return const_coeff(x)
    return symbolic_coeff(x, known_real)
end

"""
    raw_leaf(x, known_real = false)

The coefficient of a raw expression that has no polynomial form. Unlike `raw_coeff`, it
also recognizes a real expression and a lone unit phase.
"""
function raw_leaf(x::RawExpression, known_real::Bool = false)
    is_phase(x) && return phase_coeff(only(SymbolicUtils.arguments(x)))
    real_part, imaginary_part = raw_realimag(x)
    is_const_zero(imaginary_part) && return raw_coeff(real_part, true)
    return raw_coeff(x, known_real)
end

@inline function is_atom(b)
    b isa SymbolicUtils.BasicSymbolic || return false
    SymbolicUtils.issym(b) && return true
    SymbolicUtils.iscall(b) || return false
    op = SymbolicUtils.operation(b)
    op === getindex && return true
    op === expim && return true
    (op === (+) || op === (*) || op === (^) || op === (/) || op === complex) &&
        return false
    args = SymbolicUtils.arguments(b)
    length(args) == 1 || return false
    argument = only(args)
    is_atom(argument) && return true
    is_call_of(argument, *) || return false
    SymbolicUtils.symtype(argument) <: Real || return false
    return all(is_atom, SymbolicUtils.arguments(argument))
end

@inline atom_coeff(x::SymbolicUtils.BasicSymbolic) = poly_coeff(
    Poly(
        Monomial{ExactComplex}[
            Monomial{ExactComplex}(EXACT_ONE, RawExpression[x], Rational{Int}[1]),
        ],
    ),
)

function rational_power(basearg, r::Rational{Int}, x, real_root::Bool)
    if is_exact_number(basearg)
        value = exact_number(basearg)
        value >= 0 && return numeric_radical(value, r, x)
        (real_root || denominator(r) == 2) && return negative_radical(value, r, x, real_root)
    end
    real_root && return raw_leaf(x)
    base = to_cnum(basearg)
    is_native(base) && return native(to_float_scalar(native_scalar(base))^r)
    if base.tail isa Poly && length(base.tail.terms) == 1
        return monomial_power(only(base.tail.terms), r, x)
    end
    return raw_leaf(x)
end

function monomial_power(m::Monomial, r::Rational{Int}, x)
    scalar = monomial_scalar(m)
    if any(is_phase, m.syms)
        length(m.syms) == 1 && isone(scalar) &&
            throw(ArgumentError("a unit phase cannot have a fractional power"))
        return raw_leaf(x)
    end
    length(m.syms) == 1 && only(m.exps) == 1 || return raw_leaf(x)
    iszero(imag(scalar)) || return raw_leaf(x)
    value = real(scalar)
    value > 0 || return raw_leaf(x)
    magnitude = if value isa Float64
        native(ComplexF64(value^Float64(r)))
    else
        numeric_radical(Rational{BigInt}(value), r, x)
    end
    magnitude.tail isa RawSymbolicCoeff && return raw_leaf(x)
    return magnitude * tiered(monomial_terms, EXACT_ONE, m.syms, Rational{Int}[r])
end

function to_cnum(x::RawExpression)
    is_numeric_const(x) && return const_coeff(x)
    is_imaginary_unit(x) && return CNUM_IM
    SymbolicUtils.issym(x) && return atom_coeff(x)
    SymbolicUtils.iscall(x) || return raw_leaf(x)
    op = SymbolicUtils.operation(x)
    args = SymbolicUtils.arguments(x)
    if op === (+)
        return recognize_sum(args)
    elseif op === (*)
        return recognize_prod(args)
    elseif op === (^) || op === sqrt || op === cbrt
        op === (^) && length(args) == 2 && is_const_int(args[2]) &&
            return pow_cnum_integer(to_cnum(args[1]), const_int(args[2]))
        base, r = radical_form(x)
        iszero(r) && return raw_leaf(x)
        return rational_power(base, r, x, op === cbrt)
    elseif op === getindex
        return atom_coeff(x)
    elseif op === expim
        length(args) == 1 && return phase_coeff(only(args))
        return raw_leaf(x)
    elseif op === conj
        inner = to_cnum(args[1])
        (is_native(inner) || is_phase(args[1])) && return conj(inner)
        return is_atom(x) ? atom_coeff(x) : raw_leaf(x)
    elseif op === (/)
        length(args) == 2 || return raw_leaf(x)
        numerator = to_cnum(args[1])
        is_const_int(args[2]) && return numerator * to_cnum(1 // const_int(args[2]))
        return numerator / to_cnum(args[2])
    elseif op === complex
        length(args) == 2 && return cnum(Num(args[1]), Num(args[2]))
        return raw_leaf(x)
    end
    if length(args) == 1 && (is_one_at_zero(op) || is_zero_at_zero(op))
        a1 = to_cnum(only(args))
        (is_native(a1) && iszero(native_scalar(a1))) && return is_one_at_zero(op) ? CNUM_ONE : CNUM_ZERO
    end
    return is_atom(x) ? atom_coeff(x) : raw_leaf(x)
end

function fold_tree(leaf::F, x::RawExpression) where {F}
    SymbolicUtils.iscall(x) || return leaf(x)
    op = SymbolicUtils.operation(x)
    args = SymbolicUtils.arguments(x)
    if op === (+) || op === (*)
        result = fold_tree(leaf, args[1])
        for k in 2:lastindex(args)
            result = op(result, fold_tree(leaf, args[k]))
        end
        return result
    elseif op === (/) && length(args) == 2
        return fold_tree(leaf, args[1]) / fold_tree(leaf, args[2])
    elseif op === (^) && length(args) == 2 && is_const_int(args[2])
        return fold_tree(leaf, args[1])^const_int(args[2])
    end
    return leaf(x)
end

function recognize_sum(args)
    parts = Coeff[]
    znative = CNUM_ZERO
    sym = CNUM_ZERO
    have_sym = false
    for a in args
        ca = to_cnum(a)
        t = ca.tail
        if t isa Native
            znative = znative + ca
        elseif t isa Poly
            push!(parts, ca)
        else
            sym = sym + ca
            have_sym = true
        end
    end
    iszero(znative) || push!(parts, znative)
    poly = tiered(summed_terms, parts)
    return have_sym ? poly + sym : poly
end

function summed_terms(::Type{E}, parts::Vector{Coeff}) where {E}
    terms = Monomial{E}[]
    for c in parts
        t = c.tail
        if t isa Native
            push!(terms, constant_term(E, native_scalar(c)))
        else
            append!(terms, as_tier(E, t::Poly))
        end
    end
    return canonical_terms!(terms)
end

function recognize_prod(args)
    factors = Coeff[]
    other = CNUM_ONE
    have_other = false
    for a in args
        ca = to_cnum(a)
        t = ca.tail
        if t isa Native || (t isa Poly && length(t.terms) == 1)
            push!(factors, ca)
        else
            other = mul_coeff_slow(other, ca)
            have_other = true
        end
    end
    mono = tiered(product_terms, factors)
    return have_other ? mul_coeff_slow(mono, other) : mono
end

function product_terms(::Type{E}, factors::Vector{Coeff}) where {E}
    scalar::Union{ComplexF64, E} = one(E)
    syms = RawExpression[]
    exps = Rational{Int}[]
    for c in factors
        t = c.tail
        if t isa Native
            scalar = scalar_mul(scalar, as_tier(E, native_scalar(c)))
        else
            m = only(as_tier(E, t::Poly))
            scalar = scalar_mul(scalar, monomial_scalar(m))
            append!(syms, m.syms)
            append!(exps, m.exps)
        end
    end
    iszero(scalar) && return Monomial{E}[]
    return Monomial{E}[monomial(E, scalar, syms, exps)]
end

function raw_product(x::RawExpression, y::RawExpression)
    rx, restx = split_radicals(x)
    ry, resty = split_radicals(y)
    (is_native(rx) && is_native(ry)) && return x * y
    result = raw_tail(rx * ry)
    for f in restx
        result = result * f
    end
    for f in resty
        result = result * f
    end
    return result
end
