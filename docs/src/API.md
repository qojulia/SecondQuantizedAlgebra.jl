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
CollectiveNLevelSpace
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
CollectiveTransition
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

```@docs
Op
OpKind
```

```@docs
operator_name
```

```@docs
operator_index
```

```@docs
set_acts_on
```

```@docs
optype
is_destroy
is_create
is_transition
is_collective_transition
is_pauli
is_spin
is_position
is_momentum
```

### Internal representation

```@docs
Coeff
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
make_time_dependent
```

```@docs
is_average
```

```@docs
is_indexed_sum
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
to_num
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
index_slot
```

```@docs
index_name
```

```@docs
index_range
```

```@docs
index_sym
```

```@docs
rename
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

```@docs
get_sum_body
```

```@docs
indexed_sum
```
