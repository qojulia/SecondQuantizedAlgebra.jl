# Changelog

All notable changes to [`SecondQuantizedAlgebra.jl`](https://github.com/qojulia/SecondQuantizedAlgebra.jl) will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


##  [v0.11.1]

### Fixed

- Follow QuantumInterface 0.4.4's move of `tensor` and `⊗` to TensorCore, declaring TensorCore directly so the shared tensor-product API continues to load and passes explicit-import ownership checks.
- Make numeric conversion for QuantumToolbox.jl type-stable

##  [v0.11.0]

### Fixed

- Improve compact expim phases through simplification, substitution, differentiation, unitary transformations, printing, and averaging.

## [v0.10.1]

### Fixed

- `dagger(::QField)` now extends `QuantumInterface.dagger`, providing the shared QuantumOptics-ecosystem adjoint verb for operator expressions. It forwards to `qadjoint`; use `qadjoint` or `qconj` for scalar and symbolic adjoints.

## [v0.10.1]

### Added

- Exact named unitary transformations for Fock, phase-space, spin/Pauli, and ordinary N-level operators, with `conjugate`, time-aware `transform`, inversion, composition, and analytically derived gauges.
- `expim(x)` provides a canonical exact unit-phase algebra: products and integer powers add their real arguments, inverse phases cancel, and substitution, differentiation, projections, and explicit exponential/trigonometric conversion preserve the phase representation.
- Coefficient simplification reduces exact trigonometric and hyperbolic identities, including supported composite arguments and higher powers.
- A unitary-transformation guide and resonator-rotation example show how to prepare an exact moving-frame Hamiltonian for a later Floquet expansion.


## [v0.10.0]

Numeric conversion (`to_numeric`/`numeric_average`/`expect`) was redesigned to be extensible, type-stable, and multi-backend. This is a breaking release.

### Added

- A second numeric backend, [QuantumToolbox.jl](https://github.com/qutip/QuantumToolbox.jl), alongside QuantumOpticsBase.jl. Both are wired in through Julia package extensions and selected with the backend singletons `QuantumOpticsBackend()` / `QuantumToolboxBackend()`.
- A uniform Hilbert-space entry point `to_numeric(op, h::HilbertSpace, dims; backend, parameter, time_parameter, operators, adjoint_ops, op_type)`. It builds the backend basis from `dims` (Fock cutoff / spin number) and is the only form that works for both backends. The backend defaults to the single loaded one.
- Open backend hooks `numeric_operator`, `numeric_basis`, `numeric_subbasis`, `numeric_embed`, `numeric_identity`, `numeric_num_subsystems`, `numeric_assemble`, `numeric_assemble_td`, `numeric_materialize`, `numeric_expect`, and `numeric_backend` are exported. Downstream packages can implement another backend and add numeric support for existing operator roles; `OpKind` remains a closed symbolic-role set.
- Differentiable control: a `time_parameter` value may be a `(p, t) -> value` function of the solver parameter vector. On QuantumToolbox this yields a `QobjEvo` differentiable with respect to `p` (SciMLSensitivity with Enzyme/Mooncake), enabling gradient-based optimal control via `sesolve`/`mesolve`; QuantumOptics rejects the form. See the Kerr-resonator control example.

- Some TTFX motivated changes and more complicated precompile workload. [#223](https://github.com/qojulia/SecondQuantizedAlgebra.jl/issues/223)

### Changed (breaking)

- `QuantumOpticsBase` moved from a hard dependency to a weak dependency. Using the numeric API now requires loading a backend: add `using QuantumOpticsBase` (or `using QuantumToolbox`) next to `using SecondQuantizedAlgebra`. The lightweight `QuantumInterface.jl` is a new hard dependency (it owns the `⊗`/`tensor`/`expect`/`basis` generics that the algebra extends).
- The time-dependent form (`time_parameter` non-empty) returns the backend's **native** time-dependent operator: a `TimeDependentSum` (QuantumOptics) or a `QobjEvo` (QuantumToolbox), both directly consumable by `mesolve`/`master_dynamic`/`sesolve`. It is no longer a `t -> op(t)` closure.
- Time-dependent conversion accepts only `op_type=nothing` or `identity`; eager `op_type` materializers are rejected because a time-varying operator cannot be materialized once during conversion.
- Product-space `dims` must have exactly one entry per symbolic subspace. Indexed numeric layouts validate that every physical slot is unique and in range instead of silently collapsing multiple sites onto a simple basis.

### Changed

- Averages of provably Hermitian operators (`adjoint(A) == A`) now carry `symtype === Real` instead of `Number`. This gives a faster `simplify` path and makes `conj(⟨a'a⟩)` fold to `⟨a'a⟩` rather than an inert `conj(...)` wrapper; indexed sums and lifted time-dependent variables inherit the typing, which survives `substitute`. Resolves [#171](https://github.com/qojulia/SecondQuantizedAlgebra.jl/issues/171).

### Fixed

- An elementary function of a literal zero left unevaluated by Symbolics (e.g. the `exp(0)` factor produced when `exp(im*ω*t)` is Euler-expanded to `cos(ω t)*exp(0) + exp(0)*im*sin(ω t)`) is now folded to its exact value in the coefficient algebra, so it no longer leaks into printed equations. Only exact identities at argument `0` fold (`exp/cos/cosh → 1`, `sin/tan/sinh/tanh → 0`); non-zero or non-constant arguments (`exp(2)`, `sin(π)`) stay symbolic.

- `commutator` on two operator leaves no longer disagrees with `*`. Its fast path gated on a same-site test that ignored the operator name, so two distinct modes sharing one `FockSpace` (`a` and `b` from `Destroy(h, :a)`, `Destroy(h, :b)`) returned the `[a,a†] = 1` residual while `a*b' - b'*a` correctly gave `0`. The gate is now the `Equal` site comparison the per-operator hooks are documented to require.


## [v0.9.4]

### Fixed

- LaTeX rendering of an index whose name carries numeric per-slot suffixes no longer produces an invalid double subscript. As the pipeline transforms an index it accumulates suffixes through `(i::Index)(k)` (for example an index reaches `i_2_1` when a second distinct atom `i(2)` is later collapsed to its position-1 representative), and `_latex_index_suffix` emitted the name verbatim as `_{i_2_1}`. The unbraced `_2_1` is a double subscript that MathJax rejects with "Double subscripts: use braces to clarify", so the whole equation was left unrendered in Documenter pages and notebooks. The slot numbers are now joined into a single comma subscript (`i_2_1` renders as `i_{2,1}`); bare index names are unchanged.


### Added

- `CollectiveNLevelSpace` and `CollectiveTransition` provide an exact collective description of ``N`` identical multilevel systems in the permutation-symmetric subspace. A collective transition represents ``S^{ij} = \sum_k |i\rangle_k\langle j|`` and obeys the closed ``\mathfrak{su}(n)`` algebra ``[S^{ij}, S^{kl}] = \delta_{jk}S^{il} - \delta_{li}S^{kj}``, including the two-term case ``[S^{12}, S^{21}] = S^{11} - S^{22}``. Integer and symbolic levels, simple and product-space constructors, adjoints, fundamental-operator generation, Unicode/LaTeX printing, and the `is_collective_transition` predicate are supported. Collective transitions deliberately cannot be site-indexed and do not use the single-site completeness expansion: ``\sum_i S^{ii} = N I``, rather than ``I``.
- `to_numeric` converts a `CollectiveTransition` on a `QuantumOpticsBase.ManyBodyBasis` with an `NLevelBasis` one-body basis through `manybodyoperator`. This gives exact finite-``N`` dynamics in a symmetric space of dimension ``\binom{N+n-1}{n-1}`` and works unchanged inside composite bases. For two levels it reproduces the spin-``N/2`` representation, with ``S^{12}=S_+``, ``S^{21}=S_-``, and ``(S^{11}-S^{22})/2=S_z`` in the QuantumOpticsBase convention.
- `OP_COLLECTIVE_TRANSITION` extends the public `OpKind` tags with value `7`. It is appended after `OP_MOMENTUM`, preserving the values `0:6` and structural ordering keys of every existing operator role.
- A public constructor and body accessor for the moment-layer indexed-sum node, completing the read/write API alongside the existing `is_indexed_sum`, `has_sum_metadata`, `get_sum_indices`, and `get_sum_non_equal`. `SecondQuantizedAlgebra.indexed_sum(body, indices; non_equal = NonEqualPair[])` wraps an averaged `body` in a sum node over `indices` (mirroring what `average` builds for a summed `QAdd`), and `SecondQuantizedAlgebra.get_sum_body(x)` reads the summed body back out. This lets downstream code (for example cumulant expansion factorizing a summed moment into a product form) re-wrap a body while preserving its summation scope. Resolves [#209](https://github.com/qojulia/SecondQuantizedAlgebra.jl/issues/209).
