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

Quantum models describe light, atoms, and their interactions using operators:
mathematical objects whose multiplication depends on order. Deriving equations
of motion or changing a model's representation requires repeatedly applying
relations between these operators. Doing this by hand becomes tedious as
the number of modes, atomic levels, or terms grows.
SecondQuantizedAlgebra.jl automates these manipulations in Julia
[@Bezanson2017]. It supports bosonic modes, finite-level transitions, spin,
Pauli, and position-momentum operators, including composite systems and
symbolically indexed ensembles. Researchers can construct Hamiltonians,
evaluate commutators, manipulate sums, transform operators and dissipative
generators, and convert expressions to numerical operators. The package
provides a common algebra layer for quantum-optics and many-body research
software, allowing analytic derivations and numerical simulations to share
the same model.

# Statement of need

Many research workflows begin with operator identities and end with equations
or finite-dimensional simulations. A cavity coupled to an atomic ensemble,
for example, combines bosonic commutation relations with finite-level
multiplication rules and sums over atoms. Deriving an observable's equation
of motion must distinguish contributions in which two atom indices coincide
from those in which they differ. Transforming the same model into a displaced
or rotating frame must preserve both its coherent dynamics and its decay
channels. Independent implementations of these steps duplicate effort and
make conventions harder to check across projects.

SecondQuantizedAlgebra addresses this shared need by separating operator
algebra from the choice of approximation or dynamics solver. Its algebra
originated in QuantumCumulants.jl, a framework for deriving generalized
mean-field equations [@Plankensteiner2022], and was extracted into a standalone
package. This separation lets research projects reuse and extend the algebra
without depending on the cumulant-expansion workflow. The intended users are
researchers deriving quantum-optical and many-body models, and developers
building methods that require repeated symbolic operator manipulation.

# State of the field

Several established tools address related tasks. QuantumAlgebra.jl supports
noncommuting operators, symbolic indices, normal forms, and equations of
motion in Julia [@QuantumAlgebra]. SymPy provides a general symbolic-computing
environment with quantum-operator functionality [@Meurer2017]. Numerical
frameworks such as QuantumOptics.jl [@Kraemer2018] and QuantumToolbox.jl
[@QuantumToolbox] provide representations and solvers for quantum dynamics.
These tools overlap with different parts of the workflow, and the appropriate
choice depends on the operator families, symbolic conventions, and numerical
methods needed by a project.

SecondQuantizedAlgebra's contribution is the combination of explicit composite
Hilbert spaces, eagerly canonicalized arithmetic, constrained indexed sums,
symbolic transformations, and conversion to existing numerical backends.
Its standalone development preserves the operator conventions and symbolic
interfaces needed by QuantumCumulants while making them available to other
methods. Building on the existing algebra also provides a shared place to
maintain identities and index-scope rules across those methods. The package
therefore complements numerical solvers and supplies a reusable component
within the Julia quantum-software ecosystem. Its documentation includes a
reproducible comparison suite.

# Software design

The central design choice is to maintain a canonical polynomial expression
after each arithmetic operation. Products apply the relevant local operator
identities and commutation relations, then collect equal terms. Subsequent
commutators and transformations therefore operate on expressions that have
already been reduced. This choice performs simplification during expression
construction and commits users to a particular normal form. It is useful for
derivations with many intermediate expressions, where a delayed reduction
would leave redundant terms for later steps to process.

Indexed sums retain their ranges and inequality constraints as part of the
expression. When algebraic operations require separating coincident indices,
the canonicalization pipeline splits diagonal and off-diagonal contributions.
This supports calculations with a symbolic ensemble size rather than an
explicit list of atoms. For maps acting on a density operator from both sides,
one shared summation scope preserves the meaning of a collective or correlated
decay term. Ambiguous index reuse is rejected so that a transformation cannot
silently change the user's intended contractions.

Exact arithmetic and approximation have separate roles. Exact coefficients
retain rational and radical structure, while floating-point inputs remain
inexact. Formal functions of operators occupy a separate expression layer;
they are expanded into polynomial algebra only through an explicit Taylor
expansion. Similarly, a Lindblad representation preserves Hamiltonian and
decay-channel information, while a general superoperator represents the
linear map itself. These distinctions let downstream algorithms specify
where they truncate an expansion or change representation.

Numerical conversion is organized around a backend-independent core with
extensions for QuantumOpticsBase and QuantumToolbox. Researchers can derive
an expression before choosing numerical bases and cutoffs, then use the
result with existing simulation tools. Symbolic averages connect the operator
algebra to scalar symbolic expressions. Public interfaces are documented
through a manual and executable examples, including Schrieffer-Wolff
transformations, Bogoliubov transformations, collective atom-cavity models,
and optomechanics. Automated tests cover algebraic identities, indexing,
numeric conversion, and workflows that cross these boundaries; additional
quality checks examine inference, imports, and concrete field types.

# Research impact statement

SecondQuantizedAlgebra is a declared dependency of three research packages
that use operator algebra for different purposes. QuantumCumulants uses it
as the symbolic foundation for generalized mean-field calculations. The
published QuantumCumulants framework establishes the scientific motivation
and provenance of this algebra [@Plankensteiner2022]; its current standalone
dependency documents continued reuse [@QuantumCumulantsSoftware]. Research
performed with earlier versions of QuantumCumulants predates the extraction
and should be understood as the predecessor workflow's impact.

QuantumInputOutput.jl uses SecondQuantizedAlgebra for symbolic models of
systems interacting with quantum pulses, combining scattering, coupling, and
Hamiltonian descriptions with numerical dynamics [@QuantumInputOutput].
Its use in the preprint *Input-Output Analysis of Quantum Dot SUPER Excitation*
provides a concrete application to quantum-emitter dynamics driven by
traveling light pulses [@Kerber2026].
FloquetExpansions.jl uses the algebra in an ongoing research project on
high-frequency expansions for periodically driven quantum systems, deriving
effective generators and micromotion from symbolic inputs
[@FloquetExpansions]. These integrations demonstrate reuse across mean-field
equations, quantum input-output modeling, and periodically driven dynamics.
The repository also defines downstream CI jobs for all three packages,
making compatibility with these research workflows an explicit maintenance
target. The accompanying examples give researchers reproducible starting
points for applying the algebra in their own models.

# AI usage disclosure

During the development various LLM models where consulted including Antropic's claude and OpenAI GPT models.

# Acknowledgements

**AUTHOR ACTION BEFORE SUBMISSION:** Supply funding acknowledgements, grant
identifiers, other acknowledgements, and any required statement about funder
involvement, or confirm that no funding acknowledgement is applicable.

# References
