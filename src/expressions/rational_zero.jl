# Exact zero test for a sum of rational functions of the coefficient atoms. It lets a
# cancellation across quotients with different denominators, such as `D*conj(D)` against its
# expanded product, drop the term eagerly instead of waiting for `simplify`.
struct RationalForm
    num::Coeff
    den::Coeff
end

# Failure sentinel: a valid form never has a zero denominator.
const RATIONAL_FAIL = RationalForm(CNUM_ZERO, CNUM_ZERO)
# Bounds on the term pairs of a polynomial product, so the test stays cheap; exceeding one
# reports "not shown zero", never "zero". Bringing two quotients inside one operand to a
# common denominator gets the tight bound: without it, a long accumulated sum of distinct
# quotients pays a growing failed attempt on every addition.
const RATIONAL_TERM_BUDGET = 256
const RATIONAL_INNER_SUM_BUDGET = 16

@inline is_failed_form(r::RationalForm) = iszero_cnum(r.den)
@inline polynomial_term_count(c::Coeff) = c.tail isa Poly ? length(c.tail.terms) : 1
@inline within_term_budget(a::Coeff, b::Coeff, budget::Int = RATIONAL_TERM_BUDGET) =
    polynomial_term_count(a) * polynomial_term_count(b) <= budget

function checked_rational_form(num::Coeff, den::Coeff)
    (is_polynomial_coeff(num) && is_polynomial_coeff(den)) || return RATIONAL_FAIL
    polynomial_term_count(num) <= RATIONAL_TERM_BUDGET || return RATIONAL_FAIL
    polynomial_term_count(den) <= RATIONAL_TERM_BUDGET || return RATIONAL_FAIL
    return RationalForm(num, den)
end

function rational_form_mul(a::RationalForm, b::RationalForm)
    (is_failed_form(a) || is_failed_form(b)) && return RATIONAL_FAIL
    (within_term_budget(a.num, b.num) && within_term_budget(a.den, b.den)) ||
        return RATIONAL_FAIL
    return checked_rational_form(mul_cnum(a.num, b.num), mul_cnum(a.den, b.den))
end

function rational_form_add(a::RationalForm, b::RationalForm, budget::Int)
    (is_failed_form(a) || is_failed_form(b)) && return RATIONAL_FAIL
    isequal(a.den, b.den) && return checked_rational_form(add_cnum(a.num, b.num), a.den)
    (
        within_term_budget(a.num, b.den, budget) && within_term_budget(b.num, a.den, budget) &&
            within_term_budget(a.den, b.den, budget)
    ) || return RATIONAL_FAIL
    num = add_cnum(mul_cnum(a.num, b.den), mul_cnum(b.num, a.den))
    return checked_rational_form(num, mul_cnum(a.den, b.den))
end

@inline rational_form_inv(a::RationalForm) =
    (is_failed_form(a) || iszero_cnum(a.num)) ? RATIONAL_FAIL : RationalForm(a.den, a.num)

function rational_form_pow(base::RationalForm, n::Int)
    n < 0 && return rational_form_pow(rational_form_inv(base), -n)
    result = RationalForm(CNUM_ONE, CNUM_ONE)
    for _ in 1:n
        result = rational_form_mul(result, base)
        is_failed_form(result) && return RATIONAL_FAIL
    end
    return result
end

function rational_form(c::Coeff)
    tail = c.tail
    tail isa RawSymbolicCoeff && return rational_form(tail.expr)
    return RationalForm(c, CNUM_ONE)
end

function rational_form(x::RawExpression)
    is_numeric_const(x) && return checked_rational_form(const_coeff(x), CNUM_ONE)
    if SymbolicUtils.iscall(x)
        op = SymbolicUtils.operation(x)
        args = SymbolicUtils.arguments(x)
        if op === (+)
            result = rational_form(args[1])
            for k in 2:lastindex(args)
                is_failed_form(result) && return RATIONAL_FAIL
                result = rational_form_add(result, rational_form(args[k]), RATIONAL_INNER_SUM_BUDGET)
            end
            return result
        elseif op === (*)
            result = rational_form(args[1])
            for k in 2:lastindex(args)
                is_failed_form(result) && return RATIONAL_FAIL
                result = rational_form_mul(result, rational_form(args[k]))
            end
            return result
        elseif op === (/) && length(args) == 2
            num = rational_form(args[1])
            is_failed_form(num) && return RATIONAL_FAIL
            return rational_form_mul(num, rational_form_inv(rational_form(args[2])))
        elseif op === (^) && length(args) == 2 && is_const_int(args[2])
            n = const_int(args[2])
            n == typemin(Int) && return RATIONAL_FAIL
            return rational_form_pow(rational_form(args[1]), n)
        end
    end
    return checked_rational_form(recognize(x), CNUM_ONE)
end
rational_form(x) = RATIONAL_FAIL

# `true` only when `a + b` is identically zero as a rational function of its atoms, that
# is wherever it is defined. A sum whose form cannot be built is reported as not zero.
function rational_sum_is_zero(a::Coeff, b::Coeff)
    fa = rational_form(a)
    is_failed_form(fa) && return false
    sum = rational_form_add(fa, rational_form(b), RATIONAL_TERM_BUDGET)
    return !is_failed_form(sum) && iszero_cnum(sum.num)
end
