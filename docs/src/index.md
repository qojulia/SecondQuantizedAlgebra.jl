```@raw html
---
layout: home

hero:
  name: SecondQuantizedAlgebra.jl
  text: Symbolic algebra for quantum operators
  tagline: Build and manipulate second-quantized operator expressions in Julia with canonical arithmetic, indexed sums, and exact transformations.
  actions:
    - theme: brand
      text: Get started
      link: implementation/
    - theme: alt
      text: Examples
      link: examples/schrieffer_wolff/
    - theme: alt
      text: API
      link: API/
    - theme: alt
      text: View on GitHub
      link: https://github.com/qojulia/SecondQuantizedAlgebra.jl
  image:
    src: assets/logo.svg
    alt: SecondQuantizedAlgebra.jl logo

features:
  - icon: "[,]"
    title: Canonical operator algebra
    details: Apply commutation relations, local identities, normal ordering, and simplification directly to symbolic operator expressions.
  - icon: "⊗"
    title: Multiple quantum algebras
    details: Combine bosonic, N-level, Pauli, spin, and phase-space operators in composite Hilbert spaces.
  - icon: "Σ"
    title: Indexed many-body systems
    details: Work with symbolic sums and indexed operator families, including automatic diagonal splitting and free-index constraints.
  - icon: "↻"
    title: Exact transformations
    details: Construct displacement, rotation, squeezing, Bogoliubov, and generator-derived unitary transformations symbolically.
  - icon: "→"
    title: Numerical bridges
    details: Convert symbolic operators to QuantumOpticsBase or QuantumToolbox representations when numerical evaluation is needed.
---
```

```@meta
CurrentModule = SecondQuantizedAlgebra
```

`SecondQuantizedAlgebra.jl` provides the noncommutative symbolic layer used to build and transform quantum-operator expressions before numerical simulation. The algebra originated in [`QuantumCumulants.jl`](https://github.com/qojulia/QuantumCumulants.jl) and was separated into a reusable package as its scope expanded [Plankensteiner2022](@cite).

## Quick start

Install the package with Julia's package manager:

```julia-repl
pkg> add SecondQuantizedAlgebra
```

Construct a composite cavity–atom space and manipulate its operators directly:

```julia
using SecondQuantizedAlgebra

hc = FockSpace(:cavity)
ha = NLevelSpace(:atoms, 2)
h = hc ⊗ ha

@qnumbers b::Destroy(h, 1)
σ(i, j) = Transition(h, :σ, i, j, 2)

@variables g Δ
H = Δ * b' * b + g * (b * σ(2, 1) + b' * σ(1, 2))

simplify(commutator(H, b))
```

The postfix `'` follows Julia's [`Base.adjoint`](@extref Julia) convention.

The [Implementation](implementation.md) guide introduces the algebraic model and canonicalization rules. The [examples](examples/schrieffer_wolff.md) show complete workflows, while the [API](API.md) collects the exported interface.
