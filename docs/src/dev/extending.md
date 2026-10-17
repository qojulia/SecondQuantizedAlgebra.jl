```@meta
CurrentModule = SecondQuantizedAlgebra
```

# Building on the package

Packages such as QuantumCumulants.jl derive equations from the expressions of this package, and others add operator kinds or numeric backends. All of them work through the public API: an expression is read term by term, a kind is described by data, and a backend implements a fixed set of hooks.

```@contents
Pages = ["extending.md"]
Depth = 2:3
```

## Reading expressions

An operator expression is a [`QAdd`](@ref), and iterating it gives its terms as [`QTerm`](@ref) `=> coefficient` pairs. A term holds its operators in canonical order in `ops`, and the pairs of indices it requires to differ in `ne`. A coefficient becomes a Symbolics expression with [`to_num`](@ref).

```@example reading
using SecondQuantizedAlgebra
using SecondQuantizedAlgebra: to_num
h = FockSpace(:cavity) ⊗ NLevelSpace(:atom, 2)
a = Destroy(h, :a, 1)
σ(i, j) = Transition(h, :σ, i, j, 2)
@variables Δ g κ γ
H = Δ * a' * a + g * (a' * σ(1, 2) + a * σ(2, 1))
[(term.ops, to_num(c)) for (term, c) in H]
```

Every operator is an [`Op`](@ref), whatever its kind, so its name, its subsystems and its kind are read with functions.

```@example reading
using SecondQuantizedAlgebra: optype
operator_name(σ(1, 2)), acts_on(σ(1, 2)), optype(σ(1, 2))
```

Equations of motion often contain an expression and its adjoint as one unknown, ``\langle X \rangle`` and ``\langle X^\dagger \rangle = \langle X \rangle^*``. Keeping one of the two needs a choice that is the same in every session, which [`canonical_isless`](@ref) makes.

```@example reading
using SecondQuantizedAlgebra: canonical_isless
x = a' * σ(1, 2)
canonical_isless(x, x') ? x : x'
```

Indexed sums carry their summation indices, and the average of a sum its summand and excluded sites, read with the accessors shown under [Many emitters](@ref).

## Reading maps

A [`Lindbladian`](@ref) splits into its [`hamiltonian`](@ref) and its [`channels`](@ref) with their rates. A channel ``\sum \mathcal{D}[J, K]`` is read with [`jump`](@ref), [`partner`](@ref) and [`scope`](@ref).

```@example reading
using SecondQuantizedAlgebra: hamiltonian, channels, jump, partner
L = hamiltonian_action(H) + κ * dissipator(a) + γ * dissipator(σ(1, 2))
[(jump(ch), partner(ch), rate) for (ch, rate) in channels(L)]
```

Any [`SuperOperator`](@ref) iterates as sandwiches with their coefficients, and the two sides of a sandwich are read with [`left_operator`](@ref) and [`right_operator`](@ref).

```@example reading
using SecondQuantizedAlgebra: SuperOperator, left_operator, right_operator
[(left_operator(t), right_operator(t), c) for (t, c) in SuperOperator(dissipator(a))]
```

## Adding an operator kind

An algebra that no built-in operator covers is added as an [`OperatorKind`](@ref): a finite set of generators on one type of space, their adjoints, the relations that rewrite their products, and a realization in built-in operators for numerics. Its operators are made with [`CustomOperator`](@ref) and are ordinary [`Op`](@ref)s, so they need no further methods.

The generators of ``\mathfrak{su}(1,1)`` obey ``[K_0, K_\pm] = \pm K_\pm`` and ``[K_-, K_+] = 2K_0``. Ordered as ``K_+``, ``K_0``, ``K_-``, an unordered pair becomes the ordered pair plus its commutator. On a Fock space they are realized by ``K_+ = a^{\dagger 2}/2``, ``K_- = a^2/2`` and ``K_0 = (a^\dagger a + 1/2)/2``.

```@example kind
using SecondQuantizedAlgebra
using SecondQuantizedAlgebra: OperatorKind, CustomOperator

SU11 = OperatorKind(
    :su11, FockSpace, 3;
    adjoint = (3, 2, 1),
    relations = function (K)
        Kp, K0, Km = K
        return ((K0, Kp) => Kp * K0 + Kp, (Km, K0) => K0 * Km + Km, (Km, Kp) => Kp * Km + 2 * K0)
    end,
    realization = let a = Destroy(FockSpace(:mode), :a)
        (a' * a' / 2, (a' * a + 1 // 2) / 2, a * a / 2)
    end,
    suffixes = ("₊", "₀", "₋"), latex = ("_{+}", "_{0}", "_{-}"),
)

h = FockSpace(:mode)
Kp, K0, Km = (CustomOperator(SU11, h, :K, g) for g in 1:3)
Km * Kp
```

The Casimir ``C = K_0^2 - (K_+ K_- + K_- K_+)/2`` commutes with every generator.

```@example kind
C = K0 * K0 - (Kp * Km + Km * Kp) / 2
commutator(C, Kp), commutator(C, K0), commutator(C, Km)
```

The relations apply to generators of one kind with the same name and site, and all other operators act independently of them, as two modes with different names do. Rules and realization are stored as data, so operators of the new kind run through the same compiled code as built-in ones. A kind lives in the current session like an operator name, so a package creates it in `__init__`.

```@docs
OperatorKind
CustomOperator
operator_kind
operator_label
```

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

The type of the assembled object must not depend on the number of terms, so a vector-backed lazy sum is used rather than a tuple. And a time-dependent coefficient that reads the solver parameters is a two-argument function `(p, t) -> value`, recognized by its arity rather than by a wrapper type. The two bundled extensions under `ext/` are complete examples.

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
```
