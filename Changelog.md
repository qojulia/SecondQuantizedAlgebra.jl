# Changelog

All notable changes to [`SecondQuantizedAlgebra.jl`](https://github.com/qojulia/SecondQuantizedAlgebra.jl) will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## [Unreleased]

### Changed (breaking)

- The operator-kind predicates (`optype`, `is_destroy`, `is_create`, `is_transition`, `is_collective_transition`, `is_pauli`, `is_spin`, `is_position`, `is_momentum`), `is_indexed_sum` and the numeric backend interface (`NumericBackend`, `EagerOperator`, `numeric_backend` and the `numeric_*` methods) are public but no longer exported. Import them explicitly, e.g. `using SecondQuantizedAlgebra: is_destroy`.
- `canonical_isless(a, b)` replaces the public order keys (`order_key`, `term_order_key`, `qadd_order_key`, `superterm_order_key`, `superoperator_order_key`, `channel_order_key`), which are now internal. `Base.isless` is no longer defined on `Op`, `QAdd` or `QExpr`, so `sort`, `max` and `<` on operators throw; use `sort(v; lt = canonical_isless)`.
- `QTermDict`, `NumericContext` and the builders `SuperBuilder`, `LindbladianBuilder`, `accumulate!` and `build` are internal. Sum many maps with `sum` over an array.
- Custom backends implement `numeric_superoperator(be, basis, terms)`, where each term is `(c, L, R)` with `L` and `R` already backend matrices, instead of receiving the conversion context.

## [v0.14.0]

Adds symbolic dissipative dynamics and exact GKSL extraction, with breaking API changes for type stability.

### Added

- Symbolic `Lindbladian`, `SuperOperator` and `SuperExpr` maps, including collective dissipators, composition, Heisenberg generators, indexed sums and frame transformations. Use `to_numeric` for backend superoperators and `to_numeric_master` for Hamiltonian/collapse operators.
- Exact GKSL extraction through `Lindbladian(S[, basis])` and `kossakowski(S, basis)`, with `JumpBasis`, `support_basis`, `independent_basis` and `GKSLFormError`. No Hilbert-space cutoff is needed; accepting floating-point roundoff requires an explicit `rtol`.
- `coefficient` reads an operator-product or superoperator-sandwich coefficient; ambiguous constraints and summed expressions throw. `map_coefficients` transforms coefficients, while `coefficient_terms`, `CoefficientTerm` and `has_coefficient_terms` expose their scalar factors and atom powers.
- `harmonics`, `harmonic` and `secular` extract exact Fourier components from coefficients, operators and dissipative maps. Noninteger harmonics throw.
- `displacement(U, a)` reads a frame's displacement amplitude; `UnitaryTransform(G, θ)` supports symmetric hopping and squeezing generators.
- Introspection through `get_sum_indices(::QAdd)`, `n_levels`, `ground_state`, `transition_levels` and `operator_axis`. Indexed `to_numeric` accepts `parameter` and `operators` keywords.

### Changed (breaking)

- `index_slot` returns `0` instead of `nothing` for an abstract index; use `iszero(index_slot(x))`.
- `to_numeric`, `numeric_average` and `expect` default to `op_type = EagerOperator()`. Custom backends must dispatch `numeric_materialize` on `::EagerOperator` instead of `::Nothing`.
- `fundamental_operators` and `find_operators` use `names = Symbol[]` instead of `nothing` for automatic names. `find_operators` returns `Vector{Union{Op, QAdd}}` instead of `Vector{QField}`.
- `channels(L)` returns a sorted vector of `Channel => rate` pairs in `channel_order_key` order, replacing unordered dictionary pairs.
- For custom numeric backends, the `numeric_operator(be, ::Val{kind}, op, subbasis)` fallback is removed; unsupported roles throw. `numeric_assemble_td` receives plain `(p, t)` closures instead of `PControlCoeff` wrappers.
- `numeric_basis` returns a tuple of dimensions for composite QuantumToolbox states.

### Changed

- Faster symbolic coefficient arithmetic, sparse QuantumOptics conversion on composite bases, `LazyKet` expectation values and QuantumToolbox time-dependent coefficient evaluation. Runtime-dispatch checks now cover both backend extensions.
- `substitute` rejects parameter values that are neither numeric nor symbolic scalars with `ArgumentError`.

### Fixed

- Indexed sums preserve site counts after substitution or coefficient simplification, respect excluded sites and count overlapping diagonal cases once, including products of sums and dependent binders.
- Coefficient arithmetic cancels positive numeric factors in roots without changing symbolic root branches, recognizes zero sums across different denominators and accepts bare symbolic scalars alongside phases. More large-prime radicals are canonical ([#291]).
- Time-dependent squeezing requires `Squeeze(a, r, ϕ, t)`; ambiguous `Squeeze(a, r * t, t)` calls now throw.
- Indexed numeric conversion applies parameters even without sums; `numeric_average` handles constant imaginary terms. `qadjoint` and `inner_adjoint` accept averages and phase coefficients.
- Single-operator expressions compare correctly (`2a - a == a`), `find_operators` removes duplicates, and LaTeX renders zero expressions and brackets complex prefactors correctly.

## [v0.13.0]

### Changed (breaking)

- Exactness follows the numeric type: `2.0` stays floating point, including when combined with rationals. `ComplexF64` coefficients arise only from floating-point input or division by an exact zero.
- Internal coefficient storage replaces `Coeff.z` with `Coeff.slot`, accessed through `native_scalar`, and makes `Monomial` parametric in its exact tier.

### Added

- Scalar broadcasting, generic reductions and linear algebra for coefficient arrays, including `det`, adjoints, transposes, inner products and norms. LinearAlgebra is now a hard dependency.

### Fixed

- Exact arithmetic promotes beyond `Int` without overflow or conversion to floats. Raw symbolic rational coefficients and moving-frame conjugation preserve exact complex values.
- Exact numeric radicals canonicalize when every prime factor is below `2^32`; larger factors remain symbolic ([#291]).
- Coefficient and operator arrays render correctly through `latexify` and `text/latex`.

## [v0.12.0]

### Added

- Exact unitary transformations for bosonic and phase-space displacements, Bogoliubov maps, mode mixing, squeezing, spin/Pauli rotations and N-level basis changes, with composition, inversion, rotating/displacement frames and moving-frame gauges.
- Formal operator functions `sin(A)`, `cos(A)` and Hermitian-only `expim(A)`, with exact structural rewrites and unitary conjugation. Use `taylor(expr, 0:n)` to lower them to polynomials before numeric conversion.
- `get_variables` discovers scalar parameters, including filtered and buffer forms of the Symbolics interface.

### Changed (breaking)

- Rename `prefactor` to `get_prefactor`, and `unique_ops` / `unique_ops!` to `unique_up_to_adjoint` / `unique_up_to_adjoint!`.
- Replace single-term `operators(::QAdd)` with expression-wide `get_operators`, which returns unique operators in canonical order. For repeated factors and their order, iterate the expression and read each `QTerm.ops`.

### Fixed

- Type-stable QuantumToolbox conversion and standalone coefficient LaTeX rendering.
- Import `tensor` and `⊗` from TensorCore for QuantumInterface 0.4.4 compatibility; zero-argument Hilbert-space tensor products throw instead of recursing indefinitely.

## [v0.11.0]

### Fixed

- Preserve compact `expim` phases through simplification, substitution, differentiation, unitary transformations, printing and averaging.

## [v0.10.1]

### Added

- Exact named unitary transformations for Fock, phase-space, spin/Pauli and N-level operators, including composition, inversion and time-dependent gauges.
- Canonical `expim` phase arithmetic, with substitution, differentiation and exponential/trigonometric conversion; simplification of exact trigonometric and hyperbolic identities.
- A unitary-transformation guide and moving-frame resonator example.

### Fixed

- `dagger(::QField)` extends `QuantumInterface.dagger` and forwards to `qadjoint`. Use `qadjoint` or `qconj` for scalar and symbolic adjoints.

## [v0.10.0]

Numeric conversion gains extensible QuantumOpticsBase and QuantumToolbox backends.

### Added

- `QuantumOpticsBackend()` and `QuantumToolboxBackend()`, selected through `to_numeric(op, h::HilbertSpace, dims; backend, parameter, time_parameter, operators, adjoint_ops, op_type)`.
- Public `numeric_*` hooks for basis construction, operator assembly, materialization and expectation values, allowing downstream numeric backends.
- QuantumToolbox coefficients of the form `(p, t) -> value` for differentiable solver controls; QuantumOptics rejects this form.
- Expanded precompilation workloads to reduce first-call latency.

### Changed (breaking)

- Numeric APIs require explicitly loading `QuantumOpticsBase` or `QuantumToolbox`; QuantumOpticsBase is now a weak dependency and QuantumInterface a hard dependency.
- Time-dependent conversion returns native `TimeDependentSum` or `QobjEvo` objects instead of `t -> op(t)` closures. It accepts only `op_type = nothing` or `identity`.
- Product-space `dims` requires one entry per symbolic subspace. Indexed layouts reject duplicate or out-of-range physical slots.

### Changed

- Averages of provably Hermitian operators carry real symbolic type, allowing their conjugates to simplify directly.

### Fixed

- Coefficients fold elementary functions at literal zero to exact values.
- Leaf commutators distinguish different operator names on the same space, agreeing with `a*b - b*a`.

## [v0.9.4]

### Added

- `CollectiveNLevelSpace` and `CollectiveTransition` describe identical multilevel systems exactly in the permutation-symmetric subspace, including QuantumOpticsBase `ManyBodyBasis` conversion. Collective transitions cannot be site-indexed; their diagonal sum is `N I`, not `I`.
- `OP_COLLECTIVE_TRANSITION = 7`, preserving existing `OpKind` values.
- Public indexed-average construction and access through `indexed_sum(body, indices; non_equal)` and `get_sum_body`.
- Operator/index accessors and rebuilders: `operator_index`, `acts_on(::Index)`, `set_acts_on` and `rename`.

### Changed

- Commutation supports both residual terms of collective-transition algebra.
- Constructors reject string names with a clear `ArgumentError`; names must be `Symbol`s.

### Fixed

- Index suffixes such as `i_2_1` render as `i_{2,1}`, avoiding invalid LaTeX double subscripts.

### Documentation

- Indexed versus collective descriptions and a finite-N Tavis-Cummings example.

## [v0.9.3]

### Fixed

- Products crossing indexed sums preserve diagonal exclusions when the summed index occurs only in the coefficient, preventing spurious diagonal contributions.

## [v0.9.2]

### Changed (breaking)

- `to_numeric` defaults to concrete sparse operators consistently across expression forms. Use `op_type = dense` for dense output or `op_type = identity` to retain the previous lazy representation.

## [v0.9.1]

### Fixed

- Sums of averages containing different operator types print without comparison errors.

## [v0.9.0]

### Added

- Simultaneous operator and scalar substitution. Replacements are applied once; adjoint rules are inferred unless `replace_adjoint = false`.
- `to_numeric(op, basis; parameter, time_parameter, operators, adjoint_ops, op_type)` and vector conversion. Time-dependent parameters may be constants or functions, producing `t -> op(t)` closures.
- Public `to_num(::Coeff)`, `Coeff` and `CNum`, with API declarations through SciMLPublic.

### Changed

- Expressions display in Hilbert-space order, then operator-kind order.

## [v0.8.3]

### Changed

- `average` of a many-term operator is now O(n) instead of O(n²): a 500-term sum drops from about 31 ms to 9.6 ms.
- Wide symbolic prefactors build in O(n log n) instead of O(n²), with measured speedups of about 12.5× for sums and 17× for products at 800 terms.

## [v0.8.2]

### Changed

- Non-indexed commutators allocate less and run about twice as fast on measured many-mode workloads, using about five times less memory.

### Fixed

- Double conjugation returns the original coefficient, and conjugated expressions that become constants fold numerically.

## [v0.8.1]

### Added

- In-place `sum` and `reduce(+, …)` for `Vector{QAdd}`, integrated with MutableArithmetics. Use `sum([terms...])`; bare generators retain the generic fold. Measured many-mode construction is about 4–12× faster with 5–11× less memory.

### Changed

- Faster hashing of non-indexed operators reduces product, power, commutator and substitution costs without changing canonical results.

## [v0.8.0]

### Changed (breaking)

- Operator roles share one concrete `Op` type. `Destroy`, `Create`, `Transition`, `Pauli`, `Spin`, `Position` and `Momentum` remain constructors, but type dispatch such as `isa Destroy` must use role predicates (`is_destroy`, `is_transition`, etc.) or `optype(op)::OpKind`.
- Transition fields `i`, `j`, `ground_state`, `n_levels` become `l1`, `l2`, `g`, `nlev`; the Pauli/Spin `axis` field becomes `l1`.
- `Op` and `Index` use interned IDs. Replace direct `.name`, `.range` and `.sym` reads with `operator_name`, `index_name`, `index_range` and `index_sym`. Subspace indices are `Int32`; interned IDs must not be serialized across sessions.

### Changed

- Concrete operator storage and faster coefficient arithmetic reduce dispatch and allocations in products, powers and normal ordering.

## [v0.7.2]

### Changed

- Complex parameters and scalar functions stay on the native polynomial coefficient path, including conjugation. A measured three-cavity SLH construction is about 18× faster.
- Single-parameter radicals canonicalize, so `sqrt(p)^2` folds to `p`. Roots of composite expressions remain undistributed to preserve branches.

## [v0.7.1]

### Changed

- Cached term hashes and faster operator-reduction checks improve measured product and power workloads by about 1.5–1.8×, with roughly half the allocations.

## [v0.7.0]

### Changed

- `make_time_dependent` preserves per-site moments inside indexed sums instead of creating one collective unknown.
- Concrete coefficients keep numeric and parameter-polynomial arithmetic native, expanding polynomials canonically. Measured algebra workloads improve by about 2–3.5×.

### Fixed

- Averaged indexed sums retain their full body and summation scope through symbolic arithmetic; differently scoped sums no longer cancel incorrectly. Added `is_indexed_sum`; existing sum metadata accessors remain supported.

## [v0.6.5]

### Added

- `@variables η::Number` keeps a complex parameter atomic, with symbolic `conj(η)` and correct adjoint/numeric conversion, avoiding the real/imaginary expansion of `::Complex` parameters.

## [v0.6.4]

### Fixed

- `change_index` zeros `DoubleIndexedVariable` terms marked `identical = false` when their indices coincide, including occurrences nested in sums and products.

## [v0.6.3]

### Changed

- Time-dependent averages display as `⟨op⟩(t)` in text and LaTeX instead of internal variable names.

## [v0.6.2]

### Added

- `index_slot` recovers the position of a per-slot index created with `(i::Index)(k)` without parsing its name.

## [v0.6.1]

### Added

- Public `order_key`, `term_order_key` and `qadd_order_key` provide reproducible structural ordering, with equal keys exactly when expressions are `isequal`.

## [v0.6.0]

### Changed

- Averages use symbolic type `Number`; `make_time_dependent` lifts them to time-dependent unknowns while preserving average metadata.

## [v0.5.2]

### Fixed

- Commuting expressions with symbolic rational coefficients cancel exactly instead of leaving spurious terms in derived equations. Disjoint-subspace commutators avoid forming cancelling products.

## [v0.5.1]

### Changed

- Averaged expressions use angle brackets in LaTeX; averaged indexed sums retain their summation prefix.

### Fixed

- Diagonal splitting preserves index-inequality constraints, preventing double-counting in nested and collective sums.

### Documentation

- Explain explicit bound-index naming and why alpha-equivalent sums are not automatically collected ([#156]).

## [v0.5.0]

Breaking redesign of the algebraic core. The following changes apply when migrating from SecondQuantizedAlgebra v0.4 or the corresponding QuantumCumulants API.

### Renamed

| Previous API | Replacement |
|---|---|
| `reorder(expr, pairs)` | `assume_distinct_index(expr, pairs)` |
| `simplify(expr, h)` | `expand_completeness(simplify(expr))` |
| `normal_order(expr, h)` | `expand_completeness(normal_order(expr))` |
| `@cnumbers a b` | `@variables a b` |
| `@rnumbers a b` | `@variables a::Real b::Real` |
| `Parameter(:x)` | `@variables x` |

### Removed

| Removed construct | Migration |
|---|---|
| `ClusterSpace` | Use `Index` and `Σ`; pass truncation order to downstream moment expansion, e.g. `meanfield(...; order = k)`. |
| `SingleSum`, `DoubleSum`, `SpecialIndexedTerm` | Use `Σ`, returning a `QAdd`; inspect bound indices and each term's `ne` constraints. |
| `NumberedOperator`, `insert_index` | Keep indices symbolic; instantiate sites in the numeric layer. |
| `IndexedAverageSum`, `IndexedAverageDoubleSum` | Averaging handles indexed sums through symbolic metadata. |
| `QMul` | Arithmetic returns `QAdd`; iteration yields `QTerm => CNum` pairs with operator strings in `term.ops` and constraints in `term.ne`. |
| `order_by_index` | Declare unresolved index inequalities with `assume_distinct_index`. |

### Changed

- Ground-state projectors remain atomic through arithmetic, `normal_order` and `simplify`. Use `expand_completeness` to expand them explicitly.
- `simplify` also simplifies scalar coefficients and removes unused sum indices.
- `Transition` stores its ground state and level count, inferred from its `NLevelSpace`.
- Free indices with unknown relationships retain physical order; use `assume_distinct_index(expr, pairs)` to enable sorting.

### Migration example

Use `@variables` instead of `@cnumbers`, and specify the atom's subsystem when constructing transitions on a product space:

```julia
hc = FockSpace(:cavity)
ha = NLevelSpace(:atom, 2)
h = hc ⊗ ha
@variables N Δ
@qnumbers a::Destroy(h, 1)
i = Index(h, :i, N, ha)
g(i) = IndexedVariable(:g, i)
σ(α, β, k) = IndexedOperator(Transition(h, :σ, α, β, 2), k)
H = -Δ * a' * a + Σ(g(i) * (a' * σ(1, 2, i) + a * σ(2, 1, i)), i)
```

The final `2` in `Transition` selects the atom's subsystem in `h`, not its level count. Other constructor names and the commutator, averaging and numeric-conversion entry points remain available.


<!-- Links generated by Changelog.jl -->

[v0.5.0]: https://github.com/qojulia/SecondQuantizedAlgebra.jl/releases/tag/v0.5.0
[v0.5.1]: https://github.com/qojulia/SecondQuantizedAlgebra.jl/releases/tag/v0.5.1
[v0.5.2]: https://github.com/qojulia/SecondQuantizedAlgebra.jl/releases/tag/v0.5.2
[v0.6.0]: https://github.com/qojulia/SecondQuantizedAlgebra.jl/releases/tag/v0.6.0
[v0.6.1]: https://github.com/qojulia/SecondQuantizedAlgebra.jl/releases/tag/v0.6.1
[v0.6.2]: https://github.com/qojulia/SecondQuantizedAlgebra.jl/releases/tag/v0.6.2
[v0.6.3]: https://github.com/qojulia/SecondQuantizedAlgebra.jl/releases/tag/v0.6.3
[v0.6.4]: https://github.com/qojulia/SecondQuantizedAlgebra.jl/releases/tag/v0.6.4
[v0.6.5]: https://github.com/qojulia/SecondQuantizedAlgebra.jl/releases/tag/v0.6.5
[v0.7.0]: https://github.com/qojulia/SecondQuantizedAlgebra.jl/releases/tag/v0.7.0
[v0.7.1]: https://github.com/qojulia/SecondQuantizedAlgebra.jl/releases/tag/v0.7.1
[v0.7.2]: https://github.com/qojulia/SecondQuantizedAlgebra.jl/releases/tag/v0.7.2
[v0.8.0]: https://github.com/qojulia/SecondQuantizedAlgebra.jl/releases/tag/v0.8.0
[v0.8.1]: https://github.com/qojulia/SecondQuantizedAlgebra.jl/releases/tag/v0.8.1
[v0.8.2]: https://github.com/qojulia/SecondQuantizedAlgebra.jl/releases/tag/v0.8.2
[v0.8.3]: https://github.com/qojulia/SecondQuantizedAlgebra.jl/releases/tag/v0.8.3
[v0.9.0]: https://github.com/qojulia/SecondQuantizedAlgebra.jl/releases/tag/v0.9.0
[v0.9.1]: https://github.com/qojulia/SecondQuantizedAlgebra.jl/releases/tag/v0.9.1
[v0.9.2]: https://github.com/qojulia/SecondQuantizedAlgebra.jl/releases/tag/v0.9.2
[v0.9.3]: https://github.com/qojulia/SecondQuantizedAlgebra.jl/releases/tag/v0.9.3
[v0.9.4]: https://github.com/qojulia/SecondQuantizedAlgebra.jl/releases/tag/v0.9.4
[v0.10.0]: https://github.com/qojulia/SecondQuantizedAlgebra.jl/releases/tag/v0.10.0
[v0.10.1]: https://github.com/qojulia/SecondQuantizedAlgebra.jl/releases/tag/v0.10.1
[v0.11.0]: https://github.com/qojulia/SecondQuantizedAlgebra.jl/releases/tag/v0.11.0
[v0.12.0]: https://github.com/qojulia/SecondQuantizedAlgebra.jl/releases/tag/v0.12.0
[v0.13.0]: https://github.com/qojulia/SecondQuantizedAlgebra.jl/releases/tag/v0.13.0
[v0.14.0]: https://github.com/qojulia/SecondQuantizedAlgebra.jl/releases/tag/v0.14.0
[#156]: https://github.com/qojulia/SecondQuantizedAlgebra.jl/issues/156
[#291]: https://github.com/qojulia/SecondQuantizedAlgebra.jl/issues/291
