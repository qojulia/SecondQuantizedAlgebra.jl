# Cross-package benchmark contract

This file is the canonical definition of the comparison suite. Every benchmark
script in this directory implements the scenarios below and must use the
**canonical keys verbatim**; the result merger (`make_table.jl`) joins the
per-package JSON files on them.

## Packages and scripts

| Package | Language | Script | Run with |
|---|---|---|---|
| SecondQuantizedAlgebra.jl (SQA) | Julia | `julia_bench.jl` | `julia --project=benchmark benchmark/comparison/julia_bench.jl` |
| QuantumAlgebra.jl | Julia | `julia_bench.jl` (same run) | see above |
| SymPy (`sympy.physics.quantum`) | Python | `sympy_bench.py` | `python sympy_bench.py` (needs `pip install sympy`) |
| OpenFermion | Python | `openfermion_bench.py` | `python openfermion_bench.py` (needs `pip install openfermion`) |
| sneg | Mathematica | `sneg_bench.wls` | `wolframscript -f sneg_bench.wls` (needs [sneg](http://auger.ijs.si/sneg/)) |

Each script writes `results/<package>.json`:

```json
{
  "meta": {"package": "...", "version": "...", "language": "...",
           "runtime_version": "...", "date": "YYYY-MM-DD"},
  "results": {"<canonical key>": {"time_ns": 1234.5, "allocs": 42}},
  "capped": {"<canonical key>": 4.2}
}
```

`allocs` is optional (Julia only). A key absent from both `results` and
`capped` means the scenario is not expressible in that package; a key in
`capped` means a single trial evaluation exceeded the cap (value in seconds).

## Fairness contract

1. **Same canonical result, same physical input.** Every package is timed
   producing the normal-ordered canonical form of the same physical
   expression, written idiomatically in that package. Lazy packages
   (QuantumAlgebra, SymPy, OpenFermion) are timed *through* their
   normalization call; eager packages (SQA, sneg) are timed on the operation
   itself.
2. **Eager fold for repeated products.** For products of n > 2 factors the
   lazy packages normalize after each multiplication
   (`foldl((x, y) -> normal_form(x * y), factors)`), matching the incremental
   canonicalization eager packages perform. Timing a lazy package on the fully
   expanded product would penalize a workflow choice, not an algorithm.
3. **Known-answer checks before timing.** Each script verifies, in its own
   package: `a a† = 1 + a†a`, `a a a† = a†a a + 2a`, and the two-level
   completeness `σ⁻σ⁺ + σ⁺σ⁻ = 1` (in the package's own spin encoding, where
   spins are supported). A failing check aborts the run.
4. **Cap.** Any scenario whose single trial evaluation exceeds **3 s** is
   reported in `capped` instead of being timed.
5. **Disclosed helpers.** SymPy's `normal_ordered_form` does not know that
   boson and Pauli operators act on different Hilbert spaces (and hence
   commute), so `sympy_bench.py` includes a small reordering helper that sorts
   boson factors before Pauli factors in each product; its cost is part of
   SymPy's timed path. OpenFermion supports only numeric coefficients, so it
   uses ω = 1.0, J = 0.5, U = 0.25 where the others use symbols.

## Timing methodology

| Language | Harness | Statistic |
|---|---|---|
| Julia | BenchmarkTools `@benchmarkable`, 2 s budget | median time, allocations |
| Python | `timeit.Timer.autorange` + `repeat(5)` | minimum per-call time |
| Mathematica | `RepeatedTiming` | trimmed mean |

Cross-language ratios are order-of-magnitude indicators, not decimal-precise:
harnesses, statistics, and runtime warm-up behavior differ.

## Scenarios

Symbols ω_c, ω_a, g, ω, J, U, g_i, N are symbolic coefficients (numeric in
OpenFermion, see above). `a†`/`a` are bosonic ladder operators; σ operators
act on a two-level system in the package's native encoding.

### 1. Jaynes–Cummings model (`jc_*`)

`H_JC = ω_c a†a + ω_a σ⁺σ⁻ + g (a†σ⁻ + a σ⁺)` with the boson and the
two-level system on distinct Hilbert spaces.

| Key | Operation |
|---|---|
| `jc_build` | canonical form of `H_JC` built from scratch |
| `jc_H2` | canonical form of `H_JC · H_JC` |
| `jc_heisenberg` | canonical form of `[H_JC, a]` |
| `jc_nested_d1` … `jc_nested_d6` | `C₀ = σ⁻`, `Cᵢ = canonical([H_JC, Cᵢ₋₁])`; time the full recursion to depth 1, 2, …, 6. The docs table shows depths 2/4/6; all six points feed the scaling figure |
| `jc_pow_n2` … `jc_pow_n6` | canonical form of `H_JCⁿ` for n = 2, …, 6, built as an eager left fold of n factors `H_JC`. Product-width analog of the nested-commutator depth sweep (n = 2 coincides with `jc_H2`); all five points feed the scaling figure, none appear in the docs table |

### 2. Bosonic normal ordering (`fock_reorder_*`)

| Key | Operation |
|---|---|
| `fock_reorder_n2`, `fock_reorder_n4`, …, `fock_reorder_n10` | canonical form of `(a·a†)ⁿ` for n = 2, 4, 6, 8, 10, built as an eager left fold of n factors `a·a†`. The docs table shows n = 4/8; all five points feed the scaling figure |

### 3. Bose–Hubbard chain (`bose_hubbard_*`, `bh_chain_*`)

Open chain of M modes:
`H_BH = Σₖ ω a†ₖaₖ + J Σₖ (a†ₖaₖ₊₁ + a†ₖ₊₁aₖ) + (U/2) Σₖ a†ₖa†ₖaₖaₖ`.

| Key | Operation |
|---|---|
| `bose_hubbard_build` | canonical form of `H_BH` built from scratch at M = 8 |
| `bose_hubbard_H2` | canonical form of `H_BH · H_BH` at M = 8 |
| `bh_chain_M2`, `bh_chain_M4`, `bh_chain_M16` | build-cost of `H_BH` at chain length M = 2, 4, 16. Together with `bose_hubbard_build` (M = 8) these form a system-size sweep M = 2, 4, 8, 16 that feeds the scaling figure; only the M = 8 build and H² appear in the docs table |

### 4. Tavis–Cummings symbolic sum (`tavis_cummings_*`)

`H_TC = Σᵢ g_i (a†σ⁻ᵢ + a σ⁺ᵢ)` with a **symbolic** upper bound N
(only expressible in SQA and QuantumAlgebra).

| Key | Operation |
|---|---|
| `tavis_cummings_build` | canonical form of `H_TC` |
| `tavis_cummings_comm` | canonical form of `[H_TC, σ⁺ⱼσ⁻ⱼ]` for a free index j (forces the i = j / i ≠ j diagonal split) |

### 5. Mean-field expectation value (`meanfield_average`)

| Key | Operation |
|---|---|
| `meanfield_average` | symbolic expectation value ⟨H_JC⟩ (SQA `average`, QuantumAlgebra `expval`). Not expressible in SymPy/OpenFermion. sneg's `vev` is a *vacuum* expectation value, a different operation, so sneg is excluded too. |
