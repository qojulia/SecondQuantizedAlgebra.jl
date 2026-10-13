```@meta
CollapsedDocStrings = true
CurrentModule = SecondQuantizedAlgebra
```

# [Open Systems API](@id API: Open Systems)

```@docs
Lindbladian
hamiltonian_action
dissipator
liouvillian
SuperOperator
SuperExpr
spre
spost
sprepost
compose
```

## Introspection and builders

```@docs
SuperTerm
superterm
left_operator
right_operator
Channel
hamiltonian
channels
jump
partner
scope
formal_terms
poly
SuperBuilder
LindbladianBuilder
accumulate!
build
superterm_order_key
superoperator_order_key
channel_order_key
```

## GKSL coordinates

```@docs
JumpBasis
independent_basis
support_basis
kossakowski
Lindbladian(::SuperOperator, ::JumpBasis)
GKSLFormError
```

## Numerics

```@docs
to_numeric_master
numeric_superoperator
numeric_matrix
numeric_scaled_td
```
