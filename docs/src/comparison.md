# Comparison with QuantumAlgebra.jl

[QuantumAlgebra.jl](https://github.com/jfeist/QuantumAlgebra.jl) is the other
mature Julia package for symbolic second-quantized operator algebra. It and
SecondQuantizedAlgebra (SQA) overlap on bosonic ladder operators, two-level
systems, indexed sums, and expectation values, so a direct performance
comparison is meaningful — as long as it is done fairly. This page documents a
fair benchmark and the methodology behind it.

The benchmark script is [`benchmark/quantumalgebra_comparison.jl`](https://github.com/qojulia/SecondQuantizedAlgebra.jl/blob/main/benchmark/quantumalgebra_comparison.jl);
re-run it with `julia --project=benchmark benchmark/quantumalgebra_comparison.jl`.

## Methodology

!!! note "Fairness contract"
    The two packages have different evaluation strategies, so a naive
    comparison is misleading. Every row below times each package **producing
    the same canonical result from the same physical input**, written
    idiomatically in each package.

    * **SQA canonicalizes eagerly.** Every `*` and `+` already returns the
      normal-ordered, fully-simplified expression. The work happens at
      construction time, so for SQA we time the construction/arithmetic itself.
    * **QuantumAlgebra canonicalizes on demand.** Products stay symbolic until
      `normal_form` is called. So for QA we time the same arithmetic brought to
      canonical form: `normal_form` for single operations, `normal_form ∘ comm`
      at each level for commutators.

    This charges each package for the total work needed to reach a usable,
    canonical expression.

!!! warning "Eager vs lazy matters for products — and the fair choice is eager"
    For a product of many factors (`H²`, `Hⁿ`, `(a·a†)ⁿ`) QA has two workflows
    that reach the same answer but cost wildly different amounts:

    * **Lazy** (`normal_form(H^n)`): build the entire un-normal-ordered product,
      then canonicalize once. The intermediate expression explodes before it
      collapses.
    * **Eager** (`auto_normal_form(true)`, or `normal_form` after each `*`):
      canonicalize at every multiply, so intermediates never blow up. This is
      *exactly* what SQA does on every `*`.

    Comparing SQA's eager `*` against QA's **lazy** product would penalise QA for
    a workflow choice, not an algorithmic difference — e.g. `JC H⁴` is
    `83.96 ms` (2.2M allocations) lazily but `2.42 ms` (84k) eagerly. The
    head-to-head therefore uses QA's **eager** workflow for all products and
    powers, the genuine apples-to-apples comparison. (The lazy blowup is a real
    *ergonomics* difference — SQA is eager by default so users never hit it — but
    that is a usability point, not a raw-speed claim, and is reported separately.)

Further caveats:

* **Validated against known answers.** Before timing, the script checks a set of
  textbook operator identities (e.g. ``a a† = 1 + a†a``, ``a a a† = a†a a + 2a``,
  two-level completeness ``σgg + σee = 1``) *independently in each package*. This
  confirms both sides compute the correct canonical result — the premise of a
  fair comparison — rather than asserting their internal representations are
  byte-identical (they are not; see the next point).
* **Different canonical basis.** QA expands the two-level excited projector
  ``σ⁺σ⁻`` into ``σˣ/σʸ/σᶻ`` form, whereas SQA keeps the ``σ``-transition form.
  Results are physically equivalent but not term-by-term identical, which can
  shift term counts on either side.
* **Only mutually-expressible operations are compared.** See
  [SQA-only capabilities](@ref) for physics QuantumAlgebra cannot express, which
  is excluded from the timing table rather than presented as a rigged win.
* **Large cases are time-capped.** The scaling sweeps register each benchmark
  only if a single trial evaluation stays under a few seconds, so a runaway
  term-blowup on either side degrades gracefully instead of stalling the suite.

## Results

Measured with `BenchmarkTools` (median of each side), Julia 1.12,
SecondQuantizedAlgebra v0.5.1, QuantumAlgebra v1.6.0, on a single Linux
workstation (2026-06-09). QA uses its eager workflow for products and powers, as
described in the fairness contract above. `QA / SQA` is the QuantumAlgebra time
divided by the SecondQuantizedAlgebra time, so `4.1×` means QA took 4.1× as long
and `0.3×` means QA was faster. Absolute numbers are machine-dependent; the
ratios and *scaling trends* are the point.

| Benchmark | SQA | QuantumAlgebra | QA / SQA | SQA allocs | QA allocs |
|---|---:|---:|---:|---:|---:|
| **Core algebra** | | | | | |
| Jaynes–Cummings: build H | 7.5 μs | 30.4 μs | 4.1× | 231 | 921 |
| Jaynes–Cummings: H² | 111.9 μs | 288.4 μs | 2.6× | 1354 | 8831 |
| Schrieffer–Wolff: [S, V] | 67.3 μs | 111.5 μs | 1.7× | 837 | 2958 |
| Schrieffer–Wolff: [S, [S, H₀]] | 116.7 μs | 272.0 μs | 2.3× | 1655 | 10474 |
| Two cavities: build H | 9.6 μs | 15.2 μs | 1.6× | 273 | 373 |
| Two cavities: H² | 130.8 μs | 44.2 μs | 0.3× | 1868 | 1184 |
| Multi-mode 6-op chain | 32.4 μs | 17.1 μs | 0.5× | 857 | 443 |
| **Indexed sums** | | | | | |
| Tavis–Cummings Σ construction | 6.5 μs | 16.2 μs | 2.5× | 221 | 679 |
| Dicke Σ construction | 10.9 μs | 15.6 μs | 1.4× | 379 | 479 |
| Dicke diagonal-collapse [H, σ_j] | 77.6 μs | 59.8 μs | 0.8× | 1514 | 1226 |
| Double Σ_ij spin–spin | 7.5 μs | 3.6 μs | 0.5× | 221 | 197 |
| **Expectation values** | | | | | |
| ⟨a† σ⁻⟩ | 7.2 μs | 2.4 μs | 0.3× | 125 | 121 |
| ⟨H_JC⟩ | 29.6 μs | 1.6 μs | 0.1× | 350 | 81 |

### Reading the results

The two packages have genuinely different sweet spots:

* **SQA leads on building and transforming Hamiltonians.** Assembling coupled
  light–matter Hamiltonians and running Schrieffer–Wolff transformations is
  consistently faster (JC build 4.1×, H² 2.6×, SW 1.7–2.3×, Tavis–Cummings
  construction 2.5×), with far fewer allocations — SQA's eager pipeline reuses
  its dict-based term storage instead of re-deriving normal order per call.
* **QuantumAlgebra leads on pure/structured bosonic normal-ordering.** Squaring
  multi-mode bosonic Hamiltonians and the two-cavity / 6-op chains favour QA
  (0.3–0.5×): with little atomic structure, its specialized normal-ordering does
  less bookkeeping than SQA's general pipeline.
* **Expectation values look lopsided for a structural reason.** QA's `expval`
  is a lazy wrapper that defers evaluation, whereas SQA's `average` immediately
  builds the averaged symbolic object. These measure different amounts of work;
  treat the row as "constructing an expectation-value object," not end-to-end
  equation-of-motion cost.

### Scaling and large systems

Point measurements at small sizes are dominated by constant overhead. The more
honest comparison is how each package *scales*. (QA uses its eager workflow
throughout, per the fairness contract.)

![SQA vs QuantumAlgebra scaling](assets/qa_comparison_scaling.png)

| Benchmark | SQA | QuantumAlgebra | QA / SQA | SQA allocs | QA allocs |
|---|---:|---:|---:|---:|---:|
| **Fock (a·a†)ⁿ** | | | | | |
| n=2 | 10.1 μs | 6.6 μs | 0.6× | 278 | 262 |
| n=4 | 85.9 μs | 43.1 μs | 0.5× | 1705 | 1110 |
| n=6 | 290.7 μs | 147.6 μs | 0.5× | 5649 | 2642 |
| n=8 | 683.9 μs | 141.7 μs | 0.2× | 13960 | 5050 |
| n=10 | 1.37 ms | 249.9 μs | 0.2× | 29061 | 8508 |
| **Nested [H, σ] depth** | | | | | |
| depth=1 | 20.7 μs | 49.3 μs | 2.4× | 458 | 1785 |
| depth=3 | 503.3 μs | 457.1 μs | 0.9× | 6666 | 17197 |
| depth=5 | 1.78 ms | 2.86 ms | 1.6× | 23037 | 94849 |
| depth=7 | 4.06 ms | 10.72 ms | 2.6× | 53180 | 335052 |
| depth=8 | 5.71 ms | 17.93 ms | 3.1× | 74870 | 548664 |
| **Many-mode chain: H²** | | | | | |
| M=2 | 126.7 μs | 38.9 μs | 0.3× | 1852 | 1312 |
| M=4 | 699.9 μs | 159.9 μs | 0.2× | 10129 | 6254 |
| M=8 | 3.06 ms | 657.5 μs | 0.2× | 44891 | 25200 |
| M=16 | 12.59 ms | 2.56 ms | 0.2× | 187333 | 99372 |
| **Jaynes–Cummings Hⁿ** | | | | | |
| n=2 | 102.9 μs | 222.5 μs | 2.2× | 1356 | 8832 |
| n=3 | 430.0 μs | 788.6 μs | 1.8× | 5608 | 28850 |
| n=4 | 1.08 ms | 2.48 ms | 2.3× | 14055 | 83970 |
| **Dicke (3 spins) Hⁿ** | | | | | |
| n=2 | 742.6 μs | 149.8 μs | 0.2× | 9528 | 5570 |
| n=3 | 5.83 ms | 1.11 ms | 0.2× | 75757 | 39556 |
| n=4 | 23.52 ms | 4.86 ms | 0.2× | 295065 | 158742 |

What the scaling reveals:

* **Commutator nesting is SQA's clearest scaling win.** `[H, [H, … σ]]` starts
  even, then SQA pulls steadily ahead — 1.6× at depth 5 up to **3.1× at depth
  8**, with an order-of-magnitude allocation gap (75 k vs 549 k). The eager
  pipeline keeps each intermediate canonical and compact.
* **Powers of the Jaynes–Cummings Hamiltonian also favour SQA** (~2.2–2.3×
  across `n = 2…4`), the mixed boson + two-level case where SQA's term storage
  pays off even against QA's eager workflow.
* **QA wins pure/multi-mode bosonic powers.** `(a·a†)ⁿ` grows in QA's favour
  (0.6× at `n=2` to 0.2× at `n=10`), many-mode `H²` is a flat ~5× across
  `M = 2…16`, and the 3-spin Dicke `Hⁿ` is ~5× — wherever the structure is
  predominantly bosonic, QA's specialized normal-ordering scales better.

The headline: **SQA leads on building Hamiltonians and on commutators (the
operations that dominate deriving equations of motion), and the commutator lead
grows with problem size.** QuantumAlgebra's specialized normal-ordering scales
better on predominantly-bosonic products and powers. Both reach identical
physics; pick by which operations dominate your workload.

### Ergonomics: eager-by-default vs the lazy-product footgun

The fairness contract uses QA's *eager* product workflow. Its *default* lazy
workflow has a sharp footgun on high powers — building the full product before
canonicalizing:

| `JC H⁴` | time | allocations |
|---|---:|---:|
| QA lazy (`normal_form(H^4)`) | 83.96 ms | 2,199,379 |
| QA eager (`normal_form` per `*`) | 2.42 ms | 83,970 |

A 35× time and 26× memory difference from one workflow choice. SQA canonicalizes
eagerly on every `*`, so this cliff simply does not exist for its users — a
usability advantage, distinct from the raw-speed numbers above.

The `auto_normal_form(true)` toggle makes QA eager globally; on the JC build it
gives SQA **7.0 µs** vs QA-eager **29.2 µs**, confirming the build-time gap is
intrinsic to the normal-ordering cost, not the eager/lazy choice.

## [SQA-only capabilities](@id SQA-only capabilities)

QuantumAlgebra's two-level support is **Pauli / spin-½ only**. It has no general
``N``-level transition operators, no arbitrary spin-``S`` operators, and no
numeric conversion. The following SQA features therefore have no QuantumAlgebra
counterpart and are excluded from the timing table above rather than shown as
one-sided wins.

**General ``N``-level transitions** — e.g. a three-level Λ-system:

```julia
h = FockSpace(:c) ⊗ NLevelSpace(:Λ, 3, 1)
a = Destroy(h, :a, 1)
σ(i, j) = Transition(h, :σ, i, j, 2)
H = Δ_p * σ(2, 2) + (Δ_p - Δ_c) * σ(3, 3) +
    Ω_p * (a' * σ(1, 2) + a * σ(2, 1)) + Ω_c * (σ(2, 3) + σ(3, 2))
```

**Arbitrary spin-``S`` operators** (QA only has spin-½):

```julia
h = SpinSpace(:s)
Sx = Spin(h, :S, 1)   # full spin algebra, not just Pauli
Sz = Spin(h, :S, 3)
```

**Numeric conversion to QuantumOpticsBase** — turn a symbolic average into a
number on a concrete Hilbert space:

```julia
numeric_average(average(a' * s), ψ; level_map = levelmap)
```

These reflect SQA's design goal: a single algebra spanning Fock, ``N``-level,
Pauli, spin, and phase-space operators, with a bridge to numerics — at a small,
constant per-operation overhead relative to a boson-specialized library.
