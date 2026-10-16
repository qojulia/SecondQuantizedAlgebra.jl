struct RationalForm
    num::Coeff
    den::Coeff
end

const RATIONAL_FAIL = RationalForm(CNUM_ZERO, CNUM_ZERO)
const RATIONAL_TERM_BUDGET = 256
const RATIONAL_INNER_SUM_BUDGET = 16

@inline is_failed_form(r::RationalForm) = iszero(r.den)
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
    return checked_rational_form(a.num * b.num, a.den * b.den)
end

function rational_form_add(a::RationalForm, b::RationalForm, budget::Int)
    (is_failed_form(a) || is_failed_form(b)) && return RATIONAL_FAIL
    isequal(a.den, b.den) && return checked_rational_form(a.num + b.num, a.den)
    (
        within_term_budget(a.num, b.den, budget) && within_term_budget(b.num, a.den, budget) &&
            within_term_budget(a.den, b.den, budget)
    ) || return RATIONAL_FAIL
    num = a.num * b.den + b.num * a.den
    return checked_rational_form(num, a.den * b.den)
end

@inline Base.inv(a::RationalForm) =
    (is_failed_form(a) || iszero(a.num)) ? RATIONAL_FAIL : RationalForm(a.den, a.num)

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

"""
    sum_is_zero(a, b)

Whether `a + b` is identically zero: two raw expressions that are structural negatives of each
other, or two coefficients whose rational forms cancel within the term budget.
"""
function sum_is_zero(a::Coeff, b::Coeff)
    ta, tb = a.tail, b.tail
    if ta isa RawSymbolicCoeff && tb isa RawSymbolicCoeff && ta.known_real == tb.known_real
        (raw_is_negative_of(ta.expr, tb.expr) || raw_is_negative_of(tb.expr, ta.expr)) &&
            return true
    end
    fa = rational_form(a)
    is_failed_form(fa) && return false
    sum = rational_form_add(fa, rational_form(b), RATIONAL_TERM_BUDGET)
    return !is_failed_form(sum) && iszero(sum.num)
end

@inline function raw_is_negative_of(
        positive::SymbolicUtils.BasicSymbolic, negative::SymbolicUtils.BasicSymbolic,
    )
    SymbolicUtils.isconst(positive) && SymbolicUtils.isconst(negative) &&
        return is_const_zero(positive + negative)
    SymbolicUtils.iscall(negative) || return false
    op = SymbolicUtils.operation(negative)
    args = SymbolicUtils.arguments(negative)
    if op === (*) && length(args) >= 2 && is_const_minus_one(args[1])
        if is_call_of(positive, *)
            positive_args = SymbolicUtils.arguments(positive)
            length(args) - 1 == length(positive_args) || return false
            @inbounds for i in eachindex(positive_args)
                isequal(args[i + 1], positive_args[i]) || return false
            end
            return true
        end
        return length(args) == 2 && isequal(args[2], positive)
    end
    if op === (+) && is_call_of(positive, +)
        positive_args = SymbolicUtils.arguments(positive)
        length(args) == length(positive_args) || return false
        @inbounds for i in eachindex(positive_args)
            raw_is_negative_of(positive_args[i], args[i]) || return false
        end
        return true
    end
    if op === (/) && length(args) == 2 &&
            is_call_of(positive, /)
        positive_args = SymbolicUtils.arguments(positive)
        length(positive_args) == 2 || return false
        return isequal(args[2], positive_args[2]) &&
            (
            raw_is_negative_of(positive_args[1], args[1]) ||
                raw_is_negative_of(args[1], positive_args[1])
        )
    end
    return false
end
