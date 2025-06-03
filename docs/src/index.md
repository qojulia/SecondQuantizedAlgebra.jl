# SecondQuantizedAlgebra

[![docs](https://img.shields.io/badge/docs-online-blue.svg)](https://qojulia.github.io/SecondQuantizedAlgebra.jl/)
[![codecov](https://codecov.io/gh/qojulia/SecondQuantizedAlgebra.jl/branch/main/graph/badge.svg)](https://app.codecov.io/gh/qojulia/SecondQuantizedAlgebra.jl)
[![Benchmarks](https://github.com/qojulia/SecondQuantizedAlgebra.jl/actions/workflows/Benchmarks.yaml/badge.svg?branch=main)](https://qojulia.github.io/SecondQuantizedAlgebra.jl/benchmarks/)

[![Code Style: Blue](https://img.shields.io/badge/blue%20style%20-%20blue-4495d1.svg)](https://github.com/JuliaDiff/BlueStyle)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![jet](https://img.shields.io/badge/%F0%9F%9B%A9%EF%B8%8F_tested_with-JET.jl-233f9a)](https://github.com/aviatesk/JET.jl)

A Julia package for symbolic manipulation and algebraic computation with second quantized operators.  
SecondQuantizedAlgebra.jl provides a flexible framework for working with creation and annihilation operators, commutation relations, and algebraic expressions common in quantum many-body theory and quantum optics.

## Features

- Symbolic representation of bosonic fock creation and annihilation operators
- Automatic (anti-)commutation relation handling
- Algebraic simplification and normal ordering
- Support for atom and spin operators
- Extensible for custom operator types

## Installation

Install with Julia's package manager:
```julia
pkg> add SecondQuantizedAlgebra
```

## Usage

```julia
using SecondQuantizedAlgebra

# Define operators
a = Destroy(:a)
ad = Create(:a)

# Commutator
comm = commutator(ad, a)  # returns 1 for bosons, or as defined

# Normal ordering
expr = ad * a * ad
normal_ordered = normalorder(expr)
```

See the [documentation](https://qojulia.github.io/SecondQuantizedAlgebra.jl/) for more details and advanced usage.

## Contributing

Contributions and suggestions are welcome! Please open issues or pull requests on [GitHub](https://github.com/qojulia/SecondQuantizedAlgebra.jl).

## License

This project is licensed under the MIT License.