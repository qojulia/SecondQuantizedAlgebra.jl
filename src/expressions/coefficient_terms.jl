"""
    CoefficientTerm

One term `scalar * ∏ xᵢ^eᵢ` of a coefficient, as returned by [`coefficient_terms`](@ref). The
`scalar` is a `Complex{Rational{BigInt}}` when exact and a `ComplexF64` when inexact. `powers`
holds the `atom => exponent` pairs with nonzero rational exponents, possibly negative. An atom is
a parameter, a prime under a radical (`√2` is `2 => 1//2`), or a symbolic subexpression the
coefficient keeps whole, such as `sin(g)`. The order of `powers` is not specified.
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
[`coefficient_terms`](@ref) can read it. Numbers and Laurent polynomials in atoms are; a
symbolic expression that SecondQuantizedAlgebra keeps as an unexpanded tree is not.
"""
has_coefficient_terms(c::Coeff) = is_polynomial_coeff(c)
has_coefficient_terms(x::Coefficient) = has_coefficient_terms(to_cnum(x))

"""
    coefficient_terms(c) -> Vector{CoefficientTerm}

The terms of the coefficient `c`, without a zero term. A number gives one term with no powers,
and zero gives no term. Throws an `ArgumentError` when `c` has no
[`has_coefficient_terms`](@ref) form.

```jldoctest
julia> using SecondQuantizedAlgebra: coefficient_terms

julia> using Symbolics: @variables

julia> @variables p;

julia> [(t.scalar, t.powers) for t in coefficient_terms(2p^2 + 1 // 3)]
2-element Vector{Tuple{Complex{Rational{BigInt}}, Vector{Pair{Symbolics.Num, Rational{Int64}}}}}:
 (1//3 + 0//1*im, [])
 (2//1 + 0//1*im, [p => 2])
```
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
