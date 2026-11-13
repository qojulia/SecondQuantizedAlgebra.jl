# Cross-package benchmark contract

Contract: **canonical-symbolic-v2**. Results from the earlier contract are
historical and must not be merged with v2 runs.

## Workloads and timed boundaries

All algebra rows use symbolic real coefficients and native canonical inputs.
Elementary operators (including adjoints, population operators, and indexed
factors), Hilbert spaces, parameters, and input Hamiltonians are prepared
outside timing. Each package returns its own canonical form of the
same physical operator; internal bases and term counts need not match.

Construction rows time assembling and canonicalizing a Hamiltonian from
predefined elementary operators. Multiplication rows start with already
canonical Hamiltonians. Every repeated product uses a canonicalizing left
fold. The normal-ordering sweep starts with the prepared factor `a a†` and
times its repeated multiplication. Canonicalization of intermediate products
is included for all packages. QuantumAlgebra's automatic normalization is
explicitly disabled, its xyz spin basis is fixed, and the required
normalization calls are timed.

| Family | Physical input and operation | Keys |
|---|---|---|
| Jaynes–Cummings | `H = ω_c a†a + ω_a σ⁺σ⁻ + g(a†σ⁻ + aσ⁺)`; build, square, `[H,a]` | `jc_build`, `jc_H2`, `jc_heisenberg` |
| Nested commutators | `C₀ = σ⁻`, `C_d = [H,C_{d-1}]`, including the full recursion | `jc_nested_d1` through `jc_nested_d7` |
| Powers | Canonicalizing left fold for `Hⁿ` | `jc_pow_n2` through `jc_pow_n7` |
| Bosonic products | Canonicalizing left fold for `(a a†)ⁿ` | `fock_reorder_n2`, `n4`, `n6`, `n8`, `n10`, `n12`, `n14` |
| Bose–Hubbard | Open chain `H = Σω a†ₖaₖ + JΣ(a†ₖaₖ₊₁ + h.c.) + (U/2)Σa†ₖ²aₖ²`; build at M=8 and square | `bose_hubbard_build`, `bose_hubbard_H2` |
| Chain size | Same Hamiltonian, construction at M=2,4,16,32; M=8 above | `bh_chain_M2`, `M4`, `M16`, `M32` |
| Indexed ensemble | `H = Σᵢ gᵢ(a†σ⁻ᵢ + aσ⁺ᵢ)`; construct and compute `[H,σ⁺ⱼσ⁻ⱼ]` | `tavis_cummings_build`, `tavis_cummings_comm` |

The indexed comparison measures the contraction needed by the commutator;
it does not measure eliminating a sum to obtain a factor of the ensemble size.
The SQA representation carries an explicit symbolic bound, while QA uses its
native formal indexed sum.

`meanfield_average` constructs `⟨H_JC⟩`. It is an **interface workload**:
SQA creates Symbolics expressions; QA creates its internal expectation
representation. Absolute times are reported separately without a speed ratio.
This is not a benchmark of deriving, closing, or solving moment equations.

The Python and Mathematica runners do not implement every family. A missing
implementation is not evidence that a package cannot express the calculation.
SNEG uses a single fermionic mode for the two-level system; the single-site
algebra is equivalent here. SymPy includes a timed helper to commute factors
on independent boson and Pauli spaces. OpenFermion uses SymPy coefficients,
including exact `U/2`, rather than the numeric coefficients used historically.

## Untimed verification

A run aborts on failing primitive identities or scenario checks. Exact
normal-order coefficients use the independent Stirling-number formula
`(aa†)ⁿ = Σₖ S(n+1,k+1) a†ᵏaᵏ`. Representative JC nested commutators are checked
against the binomial formula for repeated commutation; JC and Bose–Hubbard
Heisenberg drifts have independent analytical targets. Julia also checks
Bose–Hubbard Hermiticity and the indexed TC contraction
`gⱼ(a†σ⁻ⱼ − aσ⁺ⱼ)`. Python checks Hermiticity and drift; SNEG checks the exact
normal-order coefficients and both drifts. These are representative family
checks, not a proof of every library operation or a complete perturbative
calculation.

## Measurement

The primary statistic remains minimum warmed runtime for continuity with the
historical suite. Compilation and validation are outside timing. Julia uses
BenchmarkTools and records raw sample times, median, allocation count and
allocated bytes. Python uses pyperf with one worker and records sample times.
SNEG uses calibrated batches. Cross-language minima have different sampling
and caching behavior, so ratios should not be interpreted as precise universal
speedups. Independent repeated publication runs remain recommended.

Default measurement budget: up to 10 seconds per scenario. A warmed single
call above 30 seconds is recorded as `capped`. Separately, the orchestrator
has a hard wall deadline for each **package job**, including loading,
validation, warmup and child workers. A timed-out job is archived as
`timed_out`, not silently omitted or interpreted as unsupported functionality.
This first workflow does not restart individual scenarios after a job timeout.

## Commands

Prerequisites: Julia, Python, and uv. Mathematica/SNEG is optional.

```sh
# Instantiate the existing Julia environment and the pinned Python lock.
python benchmark/comparison/run.py prepare

# Explicitly update the dependency environments and Python lock.
python benchmark/comparison/run.py prepare --refresh

# Check all implemented algebra families, with no timings.
python benchmark/comparison/run.py run --validate-only

# Run the core suites sequentially into a new immutable run directory.
python benchmark/comparison/run.py run

# Optional installed SNEG; SNEG_PATH may point to the directory with sneg.m.
python benchmark/comparison/run.py run --sneg

# Small timing smoke run, never eligible for publication.
python benchmark/comparison/run.py run --keys jc_heisenberg,fock_reorder_n2 --seconds 0.1

# Validate and render an archived run; --publish also updates tracked outputs.
python benchmark/comparison/run.py report benchmark/comparison/runs/RUN_ID
python benchmark/comparison/run.py report benchmark/comparison/runs/RUN_ID --publish

# Restore archived dependency environments before reproducing a run.
python benchmark/comparison/run.py prepare --from-run benchmark/comparison/runs/RUN_ID
```

`--julia /path/to/julia` before the subcommand selects an executable. `run`
also accepts `--python`, `--timeout`, and `--seconds`. Complete core coverage is
required for publication; SNEG remains optional.
No preparation command updates dependencies unless `--refresh` is explicit.
For frozen reproduction, check out the recorded SQA revision, retain/apply
any archived source changes, and use the recorded runtime versions.

## Archives and publication

Each run stores resolved Julia and Python environments, runner sources, the
local package patch, logs, run ID, contract/hash, source revision/dirty state,
CPU/thread settings, preferences, and timing settings. Setup is a separate
command and never overlaps measurement. Python cap handoff files use a unique
temporary directory shared with workers. Publication checks package completion,
expected keys, script hashes, contract hashes and run IDs before rendering.

Scenario statuses are `passed`, `capped`, `not_implemented`, or `not_run`.
Job statuses additionally distinguish `failed`, `timed_out`, and
`not_installed`. Failures block publication. The current runners do not assert
an `unsupported` status based merely on missing coverage.

The table, version information, PNG/SVG/PDF scaling figures, and relative-time
overview are generated from one archive. Publishing replaces the marked
results block in `docs/src/comparison.md`; no manual table copying is needed.
The overview includes the indexed TC result and excludes the interface row.
These artifacts are supporting research evidence for JOSS. The paper should
summarize the contribution and link the report, rather than reproduce code
examples or the entire benchmark inventory.
