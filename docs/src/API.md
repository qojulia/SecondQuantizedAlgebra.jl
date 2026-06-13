```@meta
CollapsedDocStrings = true
CurrentModule = SecondQuantizedAlgebra
```

# API

```@contents
Pages = ["API.md"]
Depth = 2:3
```

## [Hilbert Spaces](@id API: Hilbert Spaces)

```@docs
HilbertSpace
```

```@docs
ProductSpace
```

```@docs
FockSpace
```

```@docs
NLevelSpace
```

```@docs
PauliSpace
```

```@docs
SpinSpace
```

```@docs
PhaseSpace
```

```@docs
⊗
```

```@docs
tensor
```


## [Operators](@id API: Operators)

```@docs
QField
```

```@docs
QSym
```

```@docs
QTerm
```

```@docs
QTermDict
```

```@docs
QAdd
```

```@docs
@qnumbers
```

```@docs
Destroy
```

```@docs
Create
```

```@docs
Transition
```

```@docs
Pauli
```

```@docs
Spin
```

```@docs
Position
```

```@docs
Momentum
```


## [Algebra](@id API: Algebra)

```@docs
simplify
```

```@docs
normal_order
```

```@docs
normal_to_symmetric
```

```@docs
symmetric_to_normal
```

```@docs
commutator
```

```@docs
anticommutator
```

```@docs
expand
```

```@docs
expand_completeness
```

```@docs
assume_distinct_index
```

```@docs
qadjoint
```

```@docs
inner_adjoint
```


## [Average](@id API: Average)

```@docs
average
```

```@docs
undo_average
```

```@docs
is_average
```


## [Utility Functions](@id API: Utils)

```@docs
acts_on
```

```@docs
find_operators
```

```@docs
unique_ops
```

```@docs
unique_ops!
```

```@docs
fundamental_operators
```

```@docs
prefactor
```

```@docs
operators
```

```@docs
constraint_pairs
```

```@docs
sorted_arguments
```

```@docs
order_key
term_order_key
qadd_order_key
```

```@docs
substitute
```

```@docs
transition_superscript
```

```@docs
to_numeric
```

```@docs
numeric_average
```


## [Symbolic Summations](@id API: Sums)

```@docs
Index
```

```@docs
IndexedOperator
```

```@docs
IndexedVariable
```

```@docs
DoubleIndexedVariable
```

```@docs
Σ
```

```@docs
∑
```

```@docs
change_index
```

```@docs
get_indices
```

```@docs
has_index
```

```@docs
has_sum_metadata
```

```@docs
get_sum_indices
```

```@docs
get_sum_non_equal
```
