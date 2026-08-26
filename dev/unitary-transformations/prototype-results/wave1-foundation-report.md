# Wave 1 foundation report

Captured: 2026-08-26, Europe/Zurich.

All three spikes used isolated detached worktrees at commit `4a26c868` and changed no production
files. The coordinator reran all three harnesses through separate Julia MCP environments before
accepting their conclusions.

## Closure and Hilbert context

- Exact linear closures covered Fock number/squeeze, two-mode rotation/squeeze, Pauli, Spin,
  ordinary three-level transitions, and one concrete indexed site.
- Closure coordinates should reuse canonical `QTerm` identity and inequality semantics. A copied
  wrapper cost about 3.1 times the construction bytes on the N-level fixture.
- The minimal short-lived basis is `AdjointBasis{H}(space, terms, lookup)`. Action, affine offset,
  worklist, limits, and diagnostics belong to a separate analysis result.
- Explicit Hilbert context is authoritative. An `Op` cannot reconstruct untouched product factors,
  original space names, index ranges, symbolic level labels, or collective level count.
- Free indexed families remain outside the first slice; concrete indexed sites are compatible.

## Gaussian conventions and timed gauges

- Store generators in site-local interleaved order `(a1,a1',a2,a2',...)`, matching current package
  ordering and rule coverage. Construct a split-Nambu permutation only for algorithms that need
  particle/hole blocks.
- The coefficient scanner and independent commutators agreed exactly in both orderings. The
  ladder/quadrature numeric oracle agreed within `2e-14` on a fully coupled complex fixture.
- For package gauge convention `G = im*(dU'/dt)*U`, reconstruct nonscalar terms from body
  velocities `M^-1*Mdot` and `M^-1*ddot`, not spatial velocity `Mdot*M^-1`.
- The scalar gauge is unitary-lift data. Fock and quadrature rotations induce the same canonical
  map while differing by `thetadot/2`; Weyl displacements require their cocycle. Each exact
  strategy must therefore provide a concrete scalar `CNum` during construction.
- The scratch one-pass scanner proved the formula but is not production-worthy: it used 2777
  allocations versus 868 for four direct commutators. The production scanner must accumulate
  directly with coefficient helpers and a precomputed operator lookup.

## Exact response solving

- Retain scalar and explicit 2x2 solvers.
- Split the structural coupling graph before solving larger blocks. On the reproduced symbolic
  block-diagonal 8x8 fixture, two 4x4 solves took about 2.1 ms and 0.47 MB versus 12.9 ms and
  3.25 MB for one whole solve, while cross multiplication certified rational equivalence.
- Faddeev--LeVerrier is acceptable for bounded connected blocks. Validate the division-free
  certificate `A*numerator == det*forcing`; do not expect a post-quotient residual to simplify.
- A structurally nonzero symbolic determinant yields an ordinary quotient valid away from its
  pole. An isolated zero coordinate with zero forcing receives zero. Nonzero forcing on it and all
  connected structurally singular blocks are refused.
- Keep the solver transformation-private until a second real coefficient-level consumer exists.
  Bareiss remains rejected under current composite-division behavior.

## Foundation checkpoint

Wave 2 may proceed. It must pressure-test, rather than reopen without contrary evidence:

1. `QTerm` coordinate reuse and short-lived basis storage;
2. site-interleaved persistent ordering with split algorithm views;
3. strategy-supplied scalar lift alongside shared nonscalar gauge reconstruction;
4. block-first exact solving with explicit small-block paths and division-free certificates;
5. explicit Hilbert context as the generic correctness path.

The disposable spike implementations are not merged. Only the smallest concepts that survive the
Wave 2 materializer, strategy-interface, and API-pressure comparisons may enter production.
