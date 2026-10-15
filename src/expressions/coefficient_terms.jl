"""
    CoefficientTerm

One term `scalar * ∏ xᵢ^eᵢ` of a coefficient, as returned by [`coefficient_terms`](@ref). The
`scalar` is a `Complex{Rational{BigInt}}` when exact and a `ComplexF64` when inexact. `powers`
holds the `atom => exponent` pairs with nonzero rational exponents, possibly negative. An atom is
a parameter, a prime under a radical (`√2` is `2 => 1//2`), or a symbolic subexpression the
coefficient keeps whole, such as `sin(g)`. The order of `powers` is not specified.

# Examples

```jldoctest
julia> using SecondQuantizedAlgebra: coefficient_terms, CoefficientTerm

julia> using Symbolics: @variables

julia> @variables p;

julia> only(coefficient_terms(3p^2))
CoefficientTerm(3//1 + 0//1*im, Pair{Symbolics.Num, Rational{Int64}}[p => 2])
```

See also [`coefficient_terms`](@ref).
"""
struct CoefficientTerm
    scalar::Union{Complex{Rational{BigInt}}, ComplexF64}
    powers::Vector{Pair{Num, Rational{Int}}}
end

public_scalar(z::GaussianRational) = Complex{Rational{BigInt}}(real(z), imag(z))
public_scalar(z::ComplexF64) = z

function public_powers(m::Monomial)
    return Pair{Num, Rational{Int}}[Num(s) => e for (s, e) in zip(m.syms, m.exps)]
end

"""
    has_coefficient_terms(c) -> Bool

Whether the coefficient `c` is stored as a sum of [`CoefficientTerm`](@ref)s, so that
[`coefficient_terms`](@ref) can read it. Numbers and Laurent polynomials in parameters, radicals
and function calls such as `sin(q)` are. A quotient by a sum or a root of a sum, such as
`p / (1 + q)` or `sqrt(p + q)`, is kept as one symbolic expression and is not.

# Examples

```jldoctest
julia> using SecondQuantizedAlgebra: has_coefficient_terms

julia> using Symbolics: @variables

julia> @variables p q;

julia> has_coefficient_terms(2p^2 + sin(q)), has_coefficient_terms(p / (1 + q))
(true, false)
```

See also [`coefficient_terms`](@ref).
"""
has_coefficient_terms(c::Coeff) = is_polynomial_coeff(c)
has_coefficient_terms(x::Coefficient) = has_coefficient_terms(to_cnum(x))

"""
    coefficient_terms(c) -> Vector{CoefficientTerm}

The terms of the coefficient `c`, without a zero term. A number gives one term with no powers,
and zero gives no term. Throws an `ArgumentError` when `c` has no
[`has_coefficient_terms`](@ref) form.

# Examples

```jldoctest
julia> using SecondQuantizedAlgebra: coefficient_terms

julia> using Symbolics: @variables

julia> @variables p;

julia> [(t.scalar, t.powers) for t in coefficient_terms(2p^2 + 1 // 3)]
2-element Vector{Tuple{Complex{Rational{BigInt}}, Vector{Pair{Symbolics.Num, Rational{Int64}}}}}:
 (1//3 + 0//1*im, [])
 (2//1 + 0//1*im, [p => 2])
```

See also [`CoefficientTerm`](@ref), [`has_coefficient_terms`](@ref).
"""
function coefficient_terms(c::Coeff)
    has_coefficient_terms(c) || throw(
        ArgumentError("coefficient `$(to_num(c))` is not stored as a sum of terms")
    )
    if is_native(c)
        scalar = native_scalar(c)
        iszero(scalar) && return CoefficientTerm[]
        return [CoefficientTerm(public_scalar(scalar), Pair{Num, Rational{Int}}[])]
    end
    return CoefficientTerm[
        CoefficientTerm(public_scalar(monomial_scalar(m)), public_powers(m))
            for m in (c.tail::Poly).terms if !iszero(monomial_scalar(m))
    ]
end
coefficient_terms(x::Coefficient) = coefficient_terms(to_cnum(x))

"""
    coefficient(expr, term) -> Coeff

The coefficient of the operator product `term` in the polynomial `expr`, or zero when `expr`
has no such term. `term` is an operator, a product such as `a' * a`, or `1` for the scalar
part. It must canonicalize to a single product with unit coefficient, so `2a`, `a + a'` and
`a * a'` (which is `a'a + 1`) are refused. Index inequalities are part of a term: a product
that `expr` stores under other constraints, or under several, throws instead of reading as
zero; build `term` with [`assume_distinct_index`](@ref) to match a constrained product.
Summed expressions are refused.

Unlike [`coefficient_terms`](@ref), which splits a scalar into its parameter monomials,
`coefficient` reads the scalar that multiplies an operator product.

```jldoctest
julia> using SecondQuantizedAlgebra: coefficient

julia> @qnumbers a::Destroy(FockSpace(:c));

julia> H = 3 * a' * a + 2 * a' + 2 * a + 5;

julia> coefficient(H, a' * a), coefficient(H, a'), coefficient(H, 1), coefficient(H, a' * a')
(3, 2, 5, 0)
```
"""
function coefficient(expr::Union{QField, Coefficient}, term::Union{QField, Coefficient})
    q = coefficient_operand(expr, "expression")
    isempty(q.indices) || summed_coefficient_error()
    key = monomial_key(coefficient_operand(term, "term"))
    found = CNUM_ZERO
    matched = false
    for (t, c) in q.arguments
        isequal(t.ops, key.ops) || continue
        (matched || !isequal(t.ne, key.ne)) && ambiguous_coefficient_error()
        found = c
        matched = true
    end
    return found
end

coefficient_operand(x::QAdd, ::String) = x
coefficient_operand(x::QSym, ::String) = single_qadd(CNUM_ONE, Op[x])
coefficient_operand(x::Coefficient, ::String) = single_qadd(to_cnum(x), Op[])
coefficient_operand(::QField, what::String) = throw(
    ArgumentError(
        "`coefficient` needs a polynomial $what; lower a formal expression with `taylor` first"
    )
)

function monomial_key(q::QAdd)
    isempty(q.indices) || summed_coefficient_error()
    length(q.arguments) == 1 || throw(
        ArgumentError(
            "`coefficient` needs a single operator product with unit coefficient, but the " *
                "term has $(length(q.arguments)) canonical terms"
        )
    )
    t, c = first(q.arguments)
    isone(c) || throw(
        ArgumentError(
            "`coefficient` needs a single operator product with unit coefficient, but the " *
                "term has coefficient $(to_num(c))"
        )
    )
    return t
end

@noinline summed_coefficient_error() = throw(
    ArgumentError("`coefficient` does not read summed expressions; iterate the terms of the sum instead")
)

@noinline ambiguous_coefficient_error() = throw(
    ArgumentError(
        "the operator product is stored under other or several index constraints; build the " *
            "term with `assume_distinct_index` to match, or iterate the expression's terms"
    )
)
