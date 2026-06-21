# SecondQuantizedAlgebra

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://qojulia.github.io/SecondQuantizedAlgebra.jl/dev/)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://qojulia.github.io/SecondQuantizedAlgebra.jl/stable/)
[![codecov](https://codecov.io/gh/qojulia/SecondQuantizedAlgebra.jl/branch/main/graph/badge.svg)](https://app.codecov.io/gh/qojulia/SecondQuantizedAlgebra.jl)
[![Benchmarks](https://github.com/qojulia/SecondQuantizedAlgebra.jl/actions/workflows/Benchmarks.yaml/badge.svg?branch=main)](https://qojulia.github.io/SecondQuantizedAlgebra.jl/benchmark/)

[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![jet](https://img.shields.io/badge/%F0%9F%9B%A9%EF%B8%8F_tested_with-JET.jl-233f9a)](https://github.com/aviatesk/JET.jl)

A Julia package for symbolic manipulation and algebraic computation with second-quantized operators. SecondQuantizedAlgebra.jl provides a flexible framework for working with quantum operators, their commutation relations, and algebraic expressions common in quantum many-body theory and quantum optics.

The package provides:
- Bosonic, N-level, Pauli, spin, and phase-space operators in composite Hilbert spaces
- Automatic commutation relations and canonical-form arithmetic
- Normal ordering, simplification, and completeness expansion
- Symbolic summations with automatic diagonal splitting
- Averaging and numeric conversion via QuantumOpticsBase
- Extensible for custom operator types

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
