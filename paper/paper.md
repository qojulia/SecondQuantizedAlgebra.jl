---
title: 'SecondQuantizedAlgebra.jl: Reusable symbolic operator algebra for quantum optics and many-body physics'
tags:
  - Julia
  - quantum optics
  - symbolic computation
  - many-body physics
  - open quantum systems
authors:
  - name: Orjan Ameye
    affiliation: 1
  - name: Christoph Hotter
    affiliation: 2
  - name: David Plankensteiner
    affiliation: 3
affiliations:
  - name: Department of Physics, University of Konstanz, Germany
    index: 1
  - name: Department of Physics, Massachusetts Institute of Technology, Cambridge, Massachusetts, USA
    index: 2
  - name: QuEra Computing Inc., USA
    index: 3
date: 17 October 2026
bibliography: paper.bib
---

# Summary

Theoretical descriptions of quantum mechanical systems rely on noncommuting operators and their algebraic relations. Performing calculations within these models requires repeated application of operator identities, with expressions growing in complexity as the number of modes, atomic levels, and interactions increases. SecondQuantizedAlgebra.jl is a Julia package [@Bezanson2017] that automates this symbolic operator algebra for quantum optics and many-body physics. It supports bosonic, finite-level transition, spin, Pauli, and position-momentum operators in composite systems and symbolically indexed ensembles. Researchers can construct conservative and disipative systems and evaluate commutators, manipulate sums, transform operators, and convert symbolic expressions to numerical operators. These capabilities provide a shared foundation for research software, connecting analytic derivations and numerical simulations through a common model definition.

# Statement of need

Theoretical research in quantum optics and many-body physics often proceeds from an analytical derivation to numerical evaluation. Perturbative methods such as Feshbach projection, Schrieffer–Wolff transformations, and high-frequency expansions yield effective descriptions whose predictions are then studied numerically. At higher orders, these derivations involve increasingly large operator expressions, nested commutators, and cancellations. Efficient symbolic operator algebra is therefore essential for extending calculations beyond low-order approximations. Translating the resulting expressions into separate numerical implementations adds work and creates opportunities for inconsistencies between the derived and simulated models.

SecondQuantizedAlgebra.jl addresses these needs through two complementary roles: a high-performance engine for symbolic operator calculations and a common representation of quantum models for use across research packages. Researchers can manipulate Hamiltonians and dissipative generators while retaining symbolic parameters and operator structure, then convert the resulting expressions to numerical operators for QuantumOptics.jl [@Kraemer2018] or QuantumToolbox.jl [@QuantumToolbox]. This connects analytical calculations with numerical evaluation and comparison against the original model.

The package also provides a foundation for specialized methods and toolboxes. Its algebra originated in QuantumCumulants.jl [@Plankensteiner2022] and now supports applications including cumulant equations, quantum input-output modeling [@QuantumInputOutput], and Floquet expansions [@FloquetExansions]. A shared operator representation allows these packages to build on the same algebraic rules and creates a basis for combining their methods within broader workflows. The intended users are researchers performing demanding symbolic derivations and developers building specialized tools that connect analytical models to numerical simulation.

# State of the field


Several tools support symbolic calculations with quantum operators. SNEG provides second-quantized operator algebra in Mathematica, while SymPy and OpenFermion provide related functionality in Python. These tools address overlapping needs in analytical quantum physics. Our target workflow combines interactive model development with repeated operator calculations in higher-order perturbative expansions and subsequent numerical simulation. Julia provides an interactive programming environment with compiled execution and an established ecosystem for scientific computing, making it a suitable setting for connecting these stages.

Within Julia, QuantumAlgebra.jl provides closely related capabilities, including normal ordering, indexed operator sums, and expectation and cumulant equations. It also supports Symbolics.jl coefficients and numerical code generation. SecondQuantizedAlgebra.jl places particular emphasis on interoperability through the operator and expectation representations used by QuantumCumulants.jl, its ModelingToolkit-based equation workflow, and direct conversion to QuantumOptics and QuantumToolbox operators. QuantumSymbolics.jl likewise emphasizes interoperability, providing symbolic states, operators, and channels with translation to several numerical representations. Its general expression and rewriting framework differs from SecondQuantizedAlgebra’s eagerly canonicalized operator arithmetic: products are expanded, local operator identities are applied, and equal terms are collected during construction. This ensures that subsequent commutators and truncations act on already reduced expressions.

SecondQuantizedAlgebra was therefore extracted from QuantumCumulants to make this specialized algebra available independently of cumulant expansions while preserving its established symbolic interfaces. Developing the existing infrastructure retained these conventions and downstream connections, while enabling dedicated improvements to algebraic performance and numerical conversion. Its contribution is the combination of eager operator reduction, constrained indexed sums, and integration with analytical and numerical research workflows. The accompanying reproducible benchmarks characterize its performance relative to alternative tools for specific workloads.
<!-- Several established tools address related tasks. QuantumAlgebra.jl supports noncommuting operators, symbolic indices, normal forms, and equations of motion in Julia [@QuantumAlgebra]. SymPy provides a general symbolic-computing environment with quantum-operator functionality [@Meurer2017]. Numerical frameworks such as QuantumOptics.jl [@Kraemer2018] and QuantumToolbox.jl [@QuantumToolbox] provide representations and solvers for quantum dynamics. These tools overlap with different parts of the workflow, and the appropriate choice depends on the operator families, symbolic conventions, and numerical methods needed by a project.

SecondQuantizedAlgebra's contribution is the combination of explicit composite Hilbert spaces, eagerly canonicalized arithmetic, constrained indexed sums, symbolic transformations, and conversion to existing numerical backends. Its standalone development preserves the operator conventions and symbolic interfaces needed by QuantumCumulants while making them available to other methods. Building on the existing algebra also provides a shared place to maintain identities and index-scope rules across those methods. The package therefore complements numerical solvers and supplies a reusable component within the Julia quantum-software ecosystem. Its documentation includes a reproducible comparison suite. -->

# Software design

The central design choice is to maintain a canonical polynomial expression after each arithmetic operation. Products apply the relevant local operator identities and commutation relations, then collect equal terms. Subsequent commutators and transformations therefore operate on expressions that have already been reduced. This choice performs simplification during expression construction and commits users to a particular normal form. It is useful for derivations with many intermediate expressions, where a delayed reduction would leave redundant terms for later steps to process.

Indexed sums retain their ranges and inequality constraints as part of the expression. When algebraic operations require separating coincident indices, the canonicalization pipeline splits diagonal and off-diagonal contributions. This supports calculations with a symbolic ensemble size rather than an explicit list of atoms. For maps acting on a density operator from both sides, one shared summation scope preserves the meaning of a collective or correlated decay term. Ambiguous index reuse is rejected so that a transformation cannot silently change the user's intended contractions.

Exact arithmetic and approximation have separate roles. Exact coefficients retain rational and radical structure, while floating-point inputs remain inexact. Formal functions of operators occupy a separate expression layer; they are expanded into polynomial algebra only through an explicit Taylor expansion. Similarly, a Lindblad representation preserves Hamiltonian and decay-channel information, while a general superoperator represents the linear map itself. These distinctions let downstream algorithms specify where they truncate an expansion or change representation.

Numerical conversion is organized around a backend-independent core with extensions for QuantumOpticsBase and QuantumToolbox. Researchers can derive an expression before choosing numerical bases and cutoffs, then use the result with existing simulation tools. Symbolic averages connect the operator algebra to scalar symbolic expressions. Public interfaces are documented through a manual and executable examples, including Schrieffer-Wolff transformations, Bogoliubov transformations, collective atom-cavity models, and optomechanics. Automated tests cover algebraic identities, indexing, numeric conversion, and workflows that cross these boundaries; additional quality checks examine inference, imports, and concrete field types.

# Research impact statement

SecondQuantizedAlgebra is a declared dependency of three research packages that use operator algebra for different purposes. QuantumCumulants uses it as the symbolic foundation for generalized mean-field calculations. The published QuantumCumulants framework establishes the scientific motivation and provenance of this algebra [@Plankensteiner2022]; its current standalone dependency documents continued reuse [@QuantumCumulantsSoftware]. Research performed with earlier versions of QuantumCumulants predates the extraction and should be understood as the predecessor workflow's impact.

QuantumInputOutput.jl uses SecondQuantizedAlgebra for symbolic models of systems interacting with quantum pulses, combining scattering, coupling, and Hamiltonian descriptions with numerical dynamics . Its use in the preprint *Input-Output Analysis of Quantum Dot SUPER Excitation* provides a concrete application to quantum-emitter dynamics driven by traveling light pulses [@Kerber2026]. FloquetExpansions.jl uses the algebra in an ongoing research project on high-frequency expansions for periodically driven quantum systems, deriving effective generators and micromotion from symbolic inputs [@FloquetExpansions]. These integrations demonstrate reuse across mean-field equations, quantum input-output modeling, and periodically driven dynamics. The repository also defines downstream CI jobs for all three packages, making compatibility with these research workflows an explicit maintenance target. The accompanying examples give researchers reproducible starting points for applying the algebra in their own models.

# AI usage disclosure

During the development various LLM models where consulted including Antropic's claude and OpenAI GPT models.

# Acknowledgements

**AUTHOR ACTION BEFORE SUBMISSION:** Supply funding acknowledgements, grant identifiers, other acknowledgements, and any required statement about funder involvement, or confirm that no funding acknowledgement is applicable.

# References
