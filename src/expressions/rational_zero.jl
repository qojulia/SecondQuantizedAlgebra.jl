struct RationalForm
    num::Coeff
    den::Coeff
end

const RATIONAL_FAIL = RationalForm(CNUM_ZERO, CNUM_ZERO)
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

function Base.:*(a::RationalForm, b::RationalForm)
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

@inline Base.inv(a::RationalForm) =
    (is_failed_form(a) || iszero_cnum(a.num)) ? RATIONAL_FAIL : RationalForm(a.den, a.num)

function Base.:^(base::RationalForm, n::Int)
    n == typemin(Int) && return RATIONAL_FAIL
    n < 0 && return inv(base)^(-n)
    result = RationalForm(CNUM_ONE, CNUM_ONE)
    for _ in 1:n
        result *= base
        is_failed_form(result) && return RATIONAL_FAIL
    end
    return result
end

Base.:+(a::RationalForm, b::RationalForm) = rational_form_add(a, b, RATIONAL_INNER_SUM_BUDGET)
Base.:/(a::RationalForm, b::RationalForm) = a * inv(b)

function rational_form(c::Coeff)
    tail = c.tail
    tail isa RawSymbolicCoeff && return rational_form(tail.expr)
    return RationalForm(c, CNUM_ONE)
end

rational_form(x::RawExpression) = fold_tree(rational_leaf, x)
rational_leaf(x::RawExpression) = checked_rational_form(to_cnum(x), CNUM_ONE)
rational_form(x) = RATIONAL_FAIL

function rational_sum_is_zero(a::Coeff, b::Coeff)
    fa = rational_form(a)
    is_failed_form(fa) && return false
    sum = rational_form_add(fa, rational_form(b), RATIONAL_TERM_BUDGET)
    return !is_failed_form(sum) && iszero_cnum(sum.num)
end
