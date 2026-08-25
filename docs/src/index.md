# SecondQuantizedAlgebra

[![docs](https://img.shields.io/badge/docs-online-blue.svg)](https://qojulia.github.io/SecondQuantizedAlgebra.jl/)
[![codecov](https://codecov.io/gh/qojulia/SecondQuantizedAlgebra.jl/branch/main/graph/badge.svg)](https://app.codecov.io/gh/qojulia/SecondQuantizedAlgebra.jl)
[![Benchmarks](https://github.com/qojulia/SecondQuantizedAlgebra.jl/actions/workflows/Benchmarks.yaml/badge.svg?branch=main)](https://qojulia.github.io/SecondQuantizedAlgebra.jl/benchmark/)

[![Code Style: Blue](https://img.shields.io/badge/blue%20style%20-%20blue-4495d1.svg)](https://github.com/JuliaDiff/BlueStyle)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![jet](https://img.shields.io/badge/%F0%9F%9B%A9%EF%B8%8F_tested_with-JET.jl-233f9a)](https://github.com/aviatesk/JET.jl)

A Julia package for symbolic manipulation and algebraic computation with second-quantized operators. SecondQuantizedAlgebra.jl provides a flexible framework for working with quantum operators, their commutation relations, and algebraic expressions common in quantum many-body theory and quantum optics.

The package provides:
- Symbolic operators across Fock (`Destroy`, `Create`), N-level (`Transition`), Pauli (`Pauli`), spin (`Spin`), and phase-space (`Position`, `Momentum`) Hilbert spaces, composed via `⊗` into `ProductSpace`.
- Eager canonical-form arithmetic: every `*` applies (anti-)commutation, local algebraic identities, and `NLevelSpace` completeness in one pass, so the return type is always a canonical `QAdd`.
- Explicit pipeline functions when you want piecewise control: `normal_order`, `simplify`, `commutator`, `anticommutator`, `expand`, `expand_completeness`.
- Symbolic summations via `Index` and `Σ` for indexed families, with automatic diagonal splitting and `assume_distinct_index` for free-index constraints.
- Averaging to symbolic scalars via `average` / `undo_average`, and extensible numeric conversion via QuantumOpticsBase or QuantumToolbox with `to_numeric` / `numeric_average`.
- Hermitian conjugation across mixed operator + symbolic expressions via `qadjoint` (aliased as `qconj`; `dagger` extends `QuantumInterface.dagger` for operators) and the average-aware `inner_adjoint`.
- Extensible for custom operator types via five small hooks — see the developer docs.

The code was refactored out of [QuantumCumulants.jl](https://github.com/qojulia/QuantumCumulants.jl).

### Installation

Install with Julia's package manager:
```julia
pkg> add SecondQuantizedAlgebra
```

### Usage

```julia
using SecondQuantizedAlgebra

hc = FockSpace(:cavity)
ha = NLevelSpace(:atoms, 2)
h = hc ⊗ ha

@qnumbers b::Destroy(h, 1)
σ(i, j) = Transition(h, :σ, i, j, 2)

@variables g Δ

H = Δ * b' * b + g * (b * σ(2, 1) + b' * σ(1, 2))

@show b * b'
@show σ(2, 1) * σ(1, 1)

simplify(commutator(H, b))
```

See the [documentation](https://qojulia.github.io/SecondQuantizedAlgebra.jl/) for more details and advanced usage.

## Contributing

Contributions and suggestions are welcome! Please open issues or pull requests on [GitHub](https://github.com/qojulia/SecondQuantizedAlgebra.jl).

## License

This project is licensed under the MIT License.
