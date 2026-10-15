```@meta
CurrentModule = SecondQuantizedAlgebra
```

# Extending the package

The package has two extension points. A new operator kind is added inside the package, since every operator is one concrete struct. A numeric backend is added from outside, by implementing a small interface.

```@contents
Pages = ["extending.md"]
Depth = 2:3
```

## Adding an operator kind

Every operator is an `Op`, and its role is the `kind` field of type [`OpKind`](@ref). The algebra never inspects an operator directly. It talks to operators through six functions in `src/operators/operators.jl`, each a single method on `(::Op, ::Op)` that branches on `kind`:

| Function | Returns | Purpose |
| :-- | :-- | :-- |
| `site_compare(a, b, ne)` | `Less`, `Equal`, `Undetermined` or `Greater` | site relation of two operators, which drives the sort |
| `can_commute(a, b)` | `Bool` | whether a same-site pair may swap without a residual |
| `commute_pair(a, b)` | `(b, a, c1, ops1, c2, ops2)` | the swap and up to two residual terms ``c_1\,\text{ops}_1 + c_2\,\text{ops}_2`` |
| `may_reduce(a, b)` | `Bool` | whether `reduce_pair` can do anything, a cheap pre-check |
| `reduce_pair(a, b)` | `(kind, op, c)` | a local identity: the pair becomes one operator, a number, or stays |
| `ground_state_expand(op)` | `(is_gs, g, n, site)` | the completeness data of a ground-state projector |

Two residual slots exist because a collective commutator ``[S^{ij}, S^{kl}] = \delta_{jk} S^{il} - \delta_{li} S^{kj}`` has two terms; the other kinds fill only the first.

A new kind touches these places:

1. **The kind.** An `OP_*` member of `OpKind` in `src/operators/op.jl`, an `is_*` predicate, and an export.
2. **The constructor.** A function in the file of its Hilbert space that validates the arguments and returns an `Op`. Levels or axes go in the shared fields `l1`, `l2`, `g`, `nlev`; a kind that needs none leaves them zero. Write the constructor for the single space only; its `ProductSpace` methods are one-liners through `on_factor`, which checks the factor and moves the operator onto it.
3. **The algebra.** A branch in each function of the table that the kind needs. Kinds on one space that interact, like annihilation and creation, form a family in `site_family`.
4. **Adjoint.** A branch of `adjoint` in `src/operators/operators.jl`. `order_key` reads the shared fields and needs no change, as long as distinct operators of the kind differ in those fields.
5. **Numerics.** A branch of `numeric_operator` in every backend under `ext/`, or a deliberate `ArgumentError` where a backend has no representation.
6. **Printing.** The Unicode printer in `src/printing/printing.jl` and the LaTeX recipe in `src/printing/latexify_recipes.jl`.
7. **Tests.** Construction in `test/fundamentals/`, the algebra rules in `test/algebra/`, the conversion in `test/numeric/`, and the output in `test/presentation/`. The `order_key` contract test in `test/fundamentals/ordering_test.jl` loops over a list of operators, which needs one of the new kind.

An extra payload field on `Op`, for instance a `Union{Nothing, …}`, is not an option: `Op` must stay a concrete bits type, which the quality tests enforce.

## Adding a numeric backend

The numeric core in `src/numeric/` converts an expression into a list of terms, each a `ComplexF64` coefficient with a list of operator factors, and hands them to a backend. A backend is a singleton subtype of [`NumericBackend`](@ref) with methods for the following functions:

| Function | Responsibility |
| :-- | :-- |
| [`numeric_basis`](@ref)`(be, h, dims)` | build the full basis from a Hilbert space and truncations |
| [`numeric_num_subsystems`](@ref)`(be, basis)` | number of tensor factors |
| [`numeric_subbasis`](@ref)`(be, basis, slot)` | one tensor factor |
| [`numeric_operator`](@ref)`(be, op, subbasis)` | the matrix of one operator |
| [`numeric_embed`](@ref)`(be, basis, slot, leaf)` | a factor placed in the full space |
| [`numeric_identity`](@ref)`(be, basis)` | the identity on the full space |
| [`numeric_assemble`](@ref)`(be, basis, terms)` | the sum of the terms, possibly lazy |
| [`numeric_materialize`](@ref)`(be, assembled, op_type)` | the representation returned to the user |

`numeric_materialize` gives `op_type` a fixed meaning: [`EagerOperator`](@ref)`()` returns the backend's ordinary operator, `identity` returns the assembled object unchanged, and any other function is applied to it. Materialization happens once, after assembly, so a single operator, a product and a sum come out in the same representation.

A minimal backend:

```julia
import SecondQuantizedAlgebra as SQA

struct MyBackend <: SQA.NumericBackend end

SQA.numeric_basis(::MyBackend, h::SQA.HilbertSpace, dims) = error("build the basis")
SQA.numeric_num_subsystems(::MyBackend, basis) = error("count subsystems")
SQA.numeric_subbasis(::MyBackend, basis, slot::Int) = error("select a subsystem")
SQA.numeric_operator(::MyBackend, op::SQA.Op, subbasis) = error("construct a leaf")
SQA.numeric_embed(::MyBackend, basis, slot::Int, leaf) = error("embed a leaf")
SQA.numeric_identity(::MyBackend, basis) = error("construct the identity")
SQA.numeric_assemble(::MyBackend, basis, terms) = error("assemble terms")
SQA.numeric_materialize(::MyBackend, assembled, ::SQA.EagerOperator) =
    error("make an eager operator")
SQA.numeric_materialize(::MyBackend, assembled, ::typeof(identity)) = assembled
```

It is used as `to_numeric(H, h, dims; backend = MyBackend())`; only the two bundled backends are found automatically. Further features need further methods:

| Feature | Additional methods |
| :-- | :-- |
| time-dependent parameters | [`numeric_assemble_td`](@ref) |
| conversion with a state in place of a basis | [`numeric_backend`](@ref)`(state)`, `numeric_basis(state)` |
| `numeric_average` and `expect` | the state methods and [`numeric_expect`](@ref) |
| superoperators and `to_numeric_master` | [`numeric_superoperator`](@ref), [`numeric_matrix`](@ref), and [`numeric_scaled_td`](@ref) for time-dependent collapse operators |

Two rules keep a backend fast and inferable. The type of the assembled object must not depend on the number of terms, so a vector-backed lazy sum is used rather than a tuple. And a time-dependent coefficient that reads the solver parameters is a two-argument function `(p, t) -> value`, recognized by its arity rather than by a wrapper type. The two bundled extensions under `ext/` are complete examples.

```@docs
NumericBackend
numeric_basis
numeric_num_subsystems
numeric_subbasis
numeric_operator
numeric_embed
numeric_identity
numeric_assemble
numeric_assemble_td
numeric_materialize
EagerOperator
numeric_expect
numeric_backend
numeric_superoperator
numeric_matrix
numeric_scaled_td
NumericContext
```
