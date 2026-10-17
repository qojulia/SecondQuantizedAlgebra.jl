```@meta
CollapsedDocStrings = true
CurrentModule = SecondQuantizedAlgebra
```

# Internals

Docstrings of internal types and functions, for reading next to the code. They describe the current implementation and may change without notice. The rules they implement are on the [Architecture](@ref) page.

```@contents
Pages = ["internals.md"]
Depth = 2:3
```

## Terms and canonicalization

A `QAdd` stores its terms in a `QTermDict`. Every operation that changes an operator product ends in `canonicalize!`, which sorts with `partial_sort!` and then streams the product through the passes of `stream!`.

```@docs
QTermDict
SiteCmp
partial_sort!
canonicalize!
stream!
```

## Exact coefficients

A `Coeff` holds an exact number, a polynomial in the parameters, or a raw SymbolicUtils expression. The exact numbers and polynomials are built from the following types.

```@docs
GaussianRational
ExactComplex
Monomial
Poly
```

## Orderings

The total, identity-faithful ordering behind [`canonical_isless`](@ref).

```@docs
order_key
term_order_key
qadd_order_key
superterm_order_key
superoperator_order_key
channel_order_key
```

## Averages and sums

An average is a SymbolicUtils term with the operation `AvgFunc`, and the average of an indexed sum is a term whose scope is a `SumScope` argument. The accessors read that scope.

```@docs
AvgFunc
SumFunc
SumScope
is_indexed_sum
has_sum_metadata
get_sum_indices
get_sum_non_equal
get_sum_body
indexed_sum
```

## Accumulation

```@docs
SuperBuilder
LindbladianBuilder
accumulate!
build
```

## Numeric conversion

```@docs
NumericContext
```

## Superoperator storage

```@docs
superterm
poly
formal_terms
```
