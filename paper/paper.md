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
date: 16 November 2026
bibliography: paper.bib
---

# Summary

Theoretical descriptions of quantum mechanical systems rely on noncommuting operators and their algebraic relations. Performing calculations within these models requires repeated application of operator identities, with expressions growing in complexity as the number of modes, atomic levels, and interactions increases. SecondQuantizedAlgebra.jl is a Julia package [@Bezanson2017] that automates this symbolic operator algebra for quantum optics and many-body physics. It supports bosonic, finite-level transition, spin, Pauli, and position-momentum operators in composite systems and symbolically indexed ensembles. Researchers can construct conservative and dissipative systems and evaluate commutators, manipulate sums, transform operators, and convert symbolic expressions to numerical operators. These capabilities provide a shared foundation for research software, connecting analytic derivations and numerical simulations through a common model definition.

# Statement of need

Theoretical research in quantum optics and many-body physics often proceeds from an analytical derivation to numerical evaluation. Perturbative methods such as Feshbach projection, Schrieffer–Wolff transformations, and high-frequency expansions yield effective descriptions whose predictions are then studied numerically. At higher orders, these derivations involve increasingly large operator expressions, nested commutators, and cancellations. Efficient symbolic operator algebra is therefore essential for extending calculations beyond low-order approximations. Translating the resulting expressions into separate numerical implementations adds work and creates opportunities for inconsistencies between the derived and simulated models.

SecondQuantizedAlgebra.jl addresses these needs through two complementary roles: a high-performance engine for symbolic operator calculations and a common representation of quantum models for use across research packages. Researchers can manipulate Hamiltonians and dissipative generators while retaining symbolic parameters and operator structure, then convert the resulting expressions to numerical operators for QuantumOptics.jl [@Kraemer2018] or QuantumToolbox.jl [@QuantumToolbox]. This connects analytical calculations with numerical evaluation and comparison against the original model.

The package also provides a foundation for specialized methods and toolboxes. Its algebra originated in QuantumCumulants.jl [@Plankensteiner2022] and now supports applications including cumulant equations, quantum input-output modeling [@Hotter2026QIO], and Floquet expansions [@FloquetExpansions]. A shared operator representation allows these packages to build on the same algebraic rules and creates a basis for combining their methods within broader workflows. The intended users are researchers performing demanding symbolic derivations and developers building specialized tools that connect analytical models to numerical simulation.

# State of the field

Several tools support symbolic calculations with quantum operators. SNEG provides second-quantized operator algebra in Mathematica [@Zitko2011], while SymPy [@Meurer2017] and OpenFermion [@McClean2020] provide related functionality in Python. These tools address overlapping needs in analytical quantum physics. Our target workflow combines interactive model development with repeated operator calculations in higher-order perturbative expansions and subsequent numerical simulation. Julia provides an interactive programming environment with compiled execution and an established ecosystem for scientific computing, making it a suitable setting for connecting these stages.

Within Julia, QuantumAlgebra.jl [@QuantumAlgebra] provides closely related capabilities, including normal ordering, indexed operator sums, and expectation and cumulant equations. It also supports Symbolics.jl [@Gowda2021] coefficients and numerical code generation. SecondQuantizedAlgebra.jl places particular emphasis on interoperability through the operator and expectation representations used by QuantumCumulants.jl, its ModelingToolkit-based equation workflow [@Ma2021], and direct conversion to QuantumOptics and QuantumToolbox operators. QuantumSymbolics.jl [@QuantumSymbolics] likewise emphasizes interoperability, providing symbolic states, operators, and channels with translation to several numerical representations. Its general expression and rewriting framework differs from SecondQuantizedAlgebra’s eagerly canonicalized operator arithmetic: products are expanded, local operator identities are applied, and equal terms are collected during construction. This ensures that subsequent commutators and truncations act on already reduced expressions.

SecondQuantizedAlgebra was therefore extracted from QuantumCumulants to make this specialized algebra available independently of cumulant expansions. It retains the operator types and modeling conventions of its predecessor, so existing models carry over with few changes, while its algebraic core and numerical conversion were redesigned for performance. Its contribution is the combination of eager operator reduction, constrained indexed sums, and integration with analytical and numerical research workflows. The accompanying reproducible benchmarks characterize its performance relative to alternative tools for specific workloads.

# Software design

Operator algebra is mechanical but outgrows hand calculation quickly. Every product of noncommuting operators must be reordered with commutation relations, and the number of terms grows with each order of perturbation theory, each atomic level, and each mode. Numerical simulation avoids the algebra, but it fixes the parameters and truncates the Hilbert space, so it returns numbers rather than the dependence of a result on couplings and on the number of emitters. A symbolic calculator reaches high orders, keeps parameters and system size symbolic, and turns a derivation into a program that is rerun when the model changes. SecondQuantizedAlgebra.jl makes this reliable through one rule: it applies every exact rewrite while an expression is built, and changes basis, approximates, or renames a symbol, e.g., an index, only on request.

High orders require expressions to stay small and comparable. A product orders its factors, composes transitions on the same atom, e.g., $\sigma^{ij}\sigma^{kl}=\delta_{jk}\sigma^{il}$, applies commutation relations, and collects equal terms. The terms of nested commutators in a perturbative expansion therefore cancel as they appear, and expressions equal under these rules compare equal, so a moment hierarchy recognizes a moment it has already derived. We implemented the alternative, deferring all reduction to an explicit call, and removed it before release. Its intermediate results carried terms that cancel only at the end, so products grew before they shrank, and an unordered $aa^\dagger$ differed from the equal $a^\dagger a+1$, so equal terms went uncollected. In practice it was slower on every workload we compared, increasingly so as expressions grew. Eager reduction has its own cost. It expands every product even when little cancels, does work that a later projection discards, and fixes one normal order, from which the symmetric order of phase-space methods is reached only by an explicit conversion.

Eager reduction also requires canonical coefficients, so that equal values compare equal and cancelling terms vanish at once. We first represented every coefficient with Symbolics.jl [@Gowda2021], which offers no simplification suited to this purpose. Its automatic normalization is not canonical, and its explicit simplification is too costly to repeat after every operation and still incomplete. Symbolics is, in addition, designed primarily for real-valued expressions, whereas coefficients in quantum mechanics are generically complex. We therefore make coefficients canonical by construction. Exact numbers are stored as complex rationals of any size, with radicals reduced to products of prime roots, and the dependence on real or complex parameters as sums of monomials in a fixed order. Sums and products formed in this representation are already canonical, so no rewriting is needed after any operation. All other coefficients remain Symbolics expressions. Any scalar expression is thus admissible, although exact cancellation is guaranteed only for the classes above. The price is a scalar arithmetic of our own and a translation layer between the two representations, which must follow changes in Symbolics. We maintain this layer deliberately, because compatibility with Symbolics allows users to define parameters as Symbolics variables and equation-based tools of the Julia ecosystem to consume our results directly.

# Research impact statement

SecondQuantizedAlgebra.jl has been used in research since its extraction from QuantumCumulants.jl [@Plankensteiner2022] in June 2025, both directly and through packages built on it. QuantumInputOutput.jl constructs its symbolic models of quantum pulses in the SLH formalism with our algebra, and its accompanying article runs on version 0.15.0 [@Hotter2026QIO]. A study of SUPER excitation of quantum dots combines both packages with version 0.10.0 [@Kerber2026]. QuantumCumulants computes its operator algebra with our package. Through it, an independent group derived fourth-order cumulant equations for superradiant photon bunching in nanolasers [@BundgaardNielsen2026], and a non-Gaussianity criterion based on photon correlations was established [@Hotter2026].

Before the extraction, the same algebra served research within and beyond the developers' group as part of QuantumCumulants. Most of these studies derive mean-field and cumulant equations for collective light-matter systems, such as superradiant lasers, solid-state masers, and multilevel atoms in a cavity [@Hernandez2024]. Several carry the expansion to fourth and fifth order, either to test its validity [@FowlerWright2023] or to resolve nonreciprocal synchronization of quantum spins [@Nadolny2025a]. Others use the derived equations as a reference for different approximations, e.g., the discrete truncated Wigner approximation [@DTWA2023], or build new methods on the algebra itself, as in the systematic time coarse-graining of driven systems with QuantumGraining.jl [@Bello2025]. 

# AI usage disclosure

Generative AI tools assisted the development of the software, its documentation, and this manuscript. We used Anthropic's Claude models through Claude Code and OpenAI's GPT models through Codex for code suggestions, refactoring, test and benchmark scripts, documentation drafts, and language editing. All assisted output was reviewed, tested, and revised by the authors. The design decisions and scientific claims are our own, and we take full responsibility for the content.

# Acknowledgements

**AUTHOR ACTION BEFORE SUBMISSION:** Supply funding acknowledgements, grant identifiers, other acknowledgements, and any required statement about funder involvement, or confirm that no funding acknowledgement is applicable.

# References
