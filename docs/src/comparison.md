# Comparison with other packages

Several packages, across several languages, do symbolic second-quantized operator algebra. This page benchmarks SecondQuantizedAlgebra (SQA) head-to-head against four of them on a compact suite of practical tasks: building Hamiltonians, deriving Heisenberg equations of motion, nested commutators, operator powers, normal ordering, multi-mode chains, symbolic indexed sums, and mean-field expectation values.

| Package | Language | Canonicalization | Bosons | Two-level / spin | Symbolic coefficients | Indexed sums (symbolic N) | Symbolic ⟨·⟩ |
|---|---|---|---|---|---|---|---|
| SecondQuantizedAlgebra.jl | Julia | eager | ✓ | ✓ (N-level, Pauli, spin-S) | ✓ | ✓ | ✓ |
| [QuantumAlgebra.jl](https://github.com/jfeist/QuantumAlgebra.jl) | Julia | lazy (`normal_form`) | ✓ | ✓ (spin-½) | ✓ | ✓ | ✓ |
| [SymPy](https://docs.sympy.org/latest/modules/physics/quantum/index.html) (`sympy.physics.quantum`) | Python | lazy (`normal_ordered_form`) | ✓ | ✓ (Pauli) | ✓ | ✗ | ✗ |
| [OpenFermion](https://quantumai.google/openfermion) (`BosonOperator`) | Python | lazy (`normal_ordered`) | ✓ | ✗ (separate qubit type) | ✗ (numeric only) | ✗ | ✗ |
| [sneg](http://auger.ijs.si/sneg/) | Mathematica | eager (`nc`) | ✓ | ✓ (fermionic mode) | ✓ | ✗ | ✗ |

The suite definition, canonical benchmark keys, and one standalone script per language live in the repository under `benchmark/comparison/`; see [Reproducing](@ref) below.

## Methodology

!!! note "Fairness contract"
    The packages have different evaluation strategies, so a naive comparison is misleading. Every row times each package **producing the same canonical (normal-ordered) result from the same physical input**, written idiomatically in that package.

    * **Eager packages** (SQA, sneg) canonicalize on every multiplication, so the operation itself is timed.
    * **Lazy packages** (QuantumAlgebra, SymPy, OpenFermion) keep products symbolic until a normalization call, so the arithmetic is timed *together with* the normalization that reaches the canonical form the eager packages hand you for free.
    * **Products of three or more factors** use the eager workflow on the lazy side too: normalize after each multiply rather than expanding the whole product first. Timing a lazy package on the fully expanded product would penalize a workflow choice, not an algorithm (see the ergonomics note below).

    Before timing, every script checks textbook identities in its own package (``a a† = 1 + a†a``, ``a a a† = a†a a + 2a``, two-level completeness ``σ⁻σ⁺ + σ⁺σ⁻ = 1``). A failing check aborts that package's run. Any benchmark whose single trial evaluation exceeds 3 s is reported as *capped* instead of timed.

Further caveats, each disclosed rather than silently absorbed:

* **Cross-language ratios are order-of-magnitude indicators.** Julia uses BenchmarkTools (median), Python uses `timeit` (best of 5 repeats), Mathematica uses `RepeatedTiming` (trimmed mean). Harness overhead, statistics, and warm-up behavior differ across runtimes, so read `86×` as "roughly two orders of magnitude", not as a precise figure. The Julia-vs-Julia column *is* precise.
* **OpenFermion uses numeric coefficients** (ω = 1.0, J = 0.5, U = 0.25) because `BosonOperator` does not accept symbols. That is *less* work than the symbolic-coefficient arithmetic every other package performs, so OpenFermion's numbers are, if anything, flattered. It also cannot mix bosons with qubit operators in one expression, so the Jaynes–Cummings rows are not expressible.
* **SymPy needs a small disclosed helper.** `normal_ordered_form` does not know that boson and Pauli operators act on different Hilbert spaces (and hence commute), so mixed products fail to merge out of the box. The benchmark script includes a ~15-line helper that sorts boson factors before Pauli factors in each product; its cost is included in SymPy's timed path.
* **sneg encodes the two-level atom as a single fermionic mode.** For one site the fermion algebra (``f f† + f†f = 1``, ``f² = 0``) is isomorphic to the σ algebra and bosons commute with fermions, so all Jaynes–Cummings results are identical to a spin encoding; this is also sneg's idiomatic usage in quantum impurity physics.
* **Different canonical bases.** QuantumAlgebra expands ``σ⁺σ⁻`` into ``σˣ/σʸ/σᶻ`` form and SymPy works in its Pauli basis, whereas SQA keeps the transition form. Results are physically equivalent but not term-by-term identical.

## Results

Single Linux workstation, 2026-07-17: SQA 0.9.0 and QuantumAlgebra 1.6.0 on Julia 1.12.6, SymPy 1.14.0 and OpenFermion 1.8.1 on Python 3.14.6, sneg 2.0.21 on Mathematica 14.2. The SQA column is absolute time; every other cell shows that package's time with its ratio to SQA in parentheses, so `(9.1×)` means 9.1× slower than SQA. `n/a` means the scenario is not expressible in that package.

| Benchmark | SQA | QuantumAlgebra.jl | SymPy | OpenFermion | sneg |
|---|---:|---:|---:|---:|---:|
| Jaynes–Cummings: build H | 2.840 μs | 22.916 μs (8.1×) | 81.913 μs (29×) | n/a | 311.997 μs (110×) |
| Jaynes–Cummings: H² | 5.600 μs | 240.141 μs (43×) | 486.348 μs (87×) | n/a | 8.487 ms (1500×) |
| Heisenberg equation: [H, a] | 1.800 μs | 22.980 μs (13×) | 275.281 μs (150×) | n/a | 2.688 ms (1500×) |
| Nested [H, [H, … σ⁻]] depth 2 | 8.790 μs | 138.911 μs (16×) | 809.507 μs (92×) | n/a | 13.389 ms (1500×) |
| Nested [H, [H, … σ⁻]] depth 4 | 65.220 μs | 1.305 ms (20×) | 6.453 ms (99×) | n/a | 162.178 ms (2500×) |
| Nested [H, [H, … σ⁻]] depth 6 | 246.102 μs | 6.340 ms (26×) | 388.800 ms (1600×) | n/a | 721.133 ms (2900×) |
| Normal-order (a·a†)⁴ | 6.110 μs | 45.045 μs (7.4×) | 1.118 ms (180×) | 197.608 μs (32×) | 8.518 ms (1400×) |
| Normal-order (a·a†)⁸ | 64.180 μs | 149.180 μs (2.3×) | 15.432 ms (240×) | 1.254 ms (20×) | 115.303 ms (1800×) |
| Bose–Hubbard (M = 8): build H | 36.110 μs | 82.466 μs (2.3×) | 885.219 μs (25×) | 221.818 μs (6.1×) | 4.644 ms (130×) |
| Bose–Hubbard (M = 8): H² | 478.542 μs | 1.570 ms (3.3×) | capped (>3 s) | 5.602 ms (12×) | 2.500 s (5200×) |
| Tavis–Cummings Σᵢ: build H | 2.280 μs | 12.390 μs (5.4×) | n/a | n/a | n/a |
| Tavis–Cummings: [H, σ⁺ⱼσ⁻ⱼ] | 35.250 μs | 376.002 μs (11×) | n/a | n/a | n/a |
| Mean-field ⟨H⟩ | 39.630 μs | 1.530 μs (0.039×) | n/a | n/a | n/a |

![Scaling comparison](assets/comparison_scaling.png)

### Reading the results

* **SQA leads every operator-algebra row, against every package.** The closest competitor throughout is QuantumAlgebra (2.3× to 43×), which shares SQA's Julia substrate and a specialized bosonic normal-ordering engine. The general-purpose CAS approaches (SymPy, Mathematica-based sneg) sit two to three orders of magnitude behind on the same physics.
* **The lead widens exactly where real derivations spend their time.** Nested commutators (the core of Heisenberg / Schrieffer–Wolff / cumulant derivations) grow from 16× to 26× vs QuantumAlgebra between depth 2 and 6, from 92× to 1600× vs SymPy, and put SymPy's depth-6 point at 389 ms per evaluation. Eager canonicalization keeps every intermediate compact, so the cost compounds slower with depth.
* **Operator powers separate width from depth.** The top-right panel of the figure sweeps ``Hⁿ`` (an eager fold, same contract as the table). This is a single wide product-and-collect rather than a deep recursion, and it is where general-CAS normal ordering degrades fastest: SymPy's ``Hⁿ`` accelerates past 100 ms by n = 4 and climbs to 3.3 s by n = 6, while SQA and QuantumAlgebra stay near log-linear. It is the product-expansion counterpart of the nested-commutator panel and stresses like-term collection instead of cancellation.
* **Chain length probes system size, not expression size.** The bottom-right panel builds the Bose–Hubbard Hamiltonian on M = 2, 4, 8, 16 modes; unlike the other three sweeps it grows the number of Hilbert-space factors rather than the depth or width of a single expression. The Hamiltonian has ``O(M)`` terms, so every package scales near-linearly and the panel becomes a constant-factor race: SQA leads with QuantumAlgebra within ~2.5×, OpenFermion and SymPy roughly a decade higher, and sneg another decade above. Building a Hamiltonian is the one axis where the field bunches, because it is the operation with the least canonicalization work per term.
* **OpenFermion is the strongest non-Julia contender on pure bosonic work** (6.1× to 32× behind SQA), while computing with plain numeric coefficients, a strictly easier task than the symbolic arithmetic every other column performs.
* **The one row SQA loses is measuring different objects.** QuantumAlgebra's `expval` tags operators inside its own term type (a near-trivial field move), whereas SQA's `average` materializes a Symbolics.jl `Num` expression (splitting coefficients into real and imaginary parts and building a CAS term tree). That object is what feeds `substitute`, `simplify`, numeric conversion, and ModelingToolkit downstream, so the 40 μs buys the bridge into the Julia symbolics ecosystem.
* **Symbolic indexed sums are a two-package race.** Only SQA and QuantumAlgebra can express ``Σᵢ gᵢ(a†σ⁻ᵢ + aσ⁺ᵢ)`` with a *symbolic* atom number N, and the ``[H, σ⁺ⱼσ⁻ⱼ]`` diagonal split (SQA 11× faster) is the operation that turns such sums into mesoscopic equations of motion.

### Ergonomics: eager by default vs the lazy-product footgun

The contract above uses the *eager* workflow for lazy packages. Their *default* lazy workflow has a sharp footgun on powers, building the full product before canonicalizing: with QuantumAlgebra, `normal_form(H⁴)` for the Jaynes–Cummings H costs about 99 ms and 2.2 M allocations lazily versus 2.7 ms and 84 k allocations eagerly, a roughly 35× time difference from one workflow choice. SymPy's cap on Bose–Hubbard H² above is the same effect in a different package. SQA canonicalizes on every `*`, so this cliff does not exist for its users.

## SQA-only capabilities

Features with no counterpart in the packages above are excluded from the table rather than shown as one-sided wins:

**General ``N``-level transitions**, e.g. a three-level Λ-system (QuantumAlgebra and SymPy are two-level only):

```julia
h = FockSpace(:c) ⊗ NLevelSpace(:Λ, 3, 1)
a = Destroy(h, :a, 1)
σ(i, j) = Transition(h, :σ, i, j, 2)
H = Δ_p * σ(2, 2) + (Δ_p - Δ_c) * σ(3, 3) +
    Ω_p * (a' * σ(1, 2) + a * σ(2, 1)) + Ω_c * (σ(2, 3) + σ(3, 2))
```

**Arbitrary spin-``S`` operators:**

```julia
h = SpinSpace(:s)
Sx = Spin(h, :S, 1)   # full spin algebra, not just Pauli
Sz = Spin(h, :S, 3)
```

**Numeric conversion to QuantumOpticsBase**, turning a symbolic average into a number on a concrete Hilbert space:

```julia
numeric_average(average(a' * s), ψ; level_map = levelmap)
```

## Related packages not benchmarked

[QuantumCumulants.jl](https://github.com/qojulia/QuantumCumulants.jl) is the package SQA was refactored from and targets the full cumulant-expansion workflow rather than the algebra layer alone. [QNET](https://github.com/mabuchilab/QNET) and its successor [QAlgebra](https://github.com/QAlgebra/qalgebra) were the closest Python analogs (SymPy-based quantum optics algebra) but are archived and do not install on current Python. [QuAlg](https://arxiv.org/abs/2008.06467) (quantum-information Fock algebra) has been dormant since 2020. The fermionic quantum-chemistry engines ([drudge](https://github.com/DrudgeCAS/drudge), [wicked](https://github.com/fevangelista/wicked), [pdaggerq](https://github.com/edeprince3/pdaggerq)) solve a different problem (Wick contractions for coupled-cluster theory) and share no expressible benchmark with SQA, which has no fermionic operators. Maple's Physics package supports creation/annihilation algebra but is commercial and was not available for this comparison. Generic noncommutative systems (NCAlgebra, FORM) would require hand-written commutation rules, so a benchmark would measure those rules rather than the package.

## [Reproducing](@id Reproducing)

Each package has a standalone script in `benchmark/comparison/` that validates the known-answer identities, runs its scenarios, and writes `results/<package>.json`; `make_table.jl` merges whatever JSONs exist into the table and figure above. The scenario definitions and canonical keys live in `benchmark/comparison/BENCHMARKS.md`.

```sh
# Julia (SQA + QuantumAlgebra)
julia --project=benchmark benchmark/comparison/julia_bench.jl

# Python (each needs: pip install sympy openfermion)
python benchmark/comparison/sympy_bench.py
python benchmark/comparison/openfermion_bench.py

# Mathematica (needs sneg installed)
wolframscript -f benchmark/comparison/sneg_bench.wls

# Merge into the docs table + figure
julia --project=benchmark benchmark/comparison/make_table.jl
```
