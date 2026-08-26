# Wave 0 baseline

Captured: 2026-08-25, Europe/Zurich.

## Repository

- Branch: `displace`
- Upstream: `origin/displace`
- Commit: `4a26c868a6a76238e46b24805e2f98fde69004d3`
- Commit subject: `displace quadrature phase space`
- Tracked working tree: clean
- Julia: 1.12.6
- Julia threads: 12

Manifest hashes:

```text
45adbfbad6fc1412c3736c754be61b664e10dc779d7639beca4665ecefec4f55  Project.toml
e7e71205c626dcfa798b07500545727065843d0518cf59a88cc56fcbe81e72bd  Manifest.toml
16e27b8bd3854bd5a819fe55f29e1d6a9eab4486c6aac431dd8409772f07e52b  test/Project.toml
f0de10bcd7db904df682487c94093093c1796c4ded1cb699b2313fa2ddab3e3f  test/Manifest.toml
```

## Behavior baseline

The focused `test/arithmetics/unitary_test.jl` suite passed through the Julia MCP:

```text
Exact unitary transformations | 277 passed / 277 total | 39.4 s
```

The frozen fixtures include named Fock and quadrature displacements, rotations, one- and two-mode
squeezes, beamsplitters, Pauli/Spin rotations, N-level matrix rotations, complete timed gauges,
automatic static and bounded displacements, exact refusal cases, inference, concrete storage, and
allocation gates.

## Complete workflow benchmark baseline

Ten samples per workflow, one evaluation per sample. Values are warm medians after package load;
time is nanoseconds.

| Workflow | Time | Memory | Allocations |
| --- | ---: | ---: | ---: |
| Automatic Fock displacement, one tone | 490358 | 175312 | 4624 |
| Automatic Fock displacement, thirty-three sidebands | 9141094 | 2385840 | 63771 |
| Automatic quadrature displacement, one tone | 656582 | 181040 | 4879 |
| Automatic quadrature displacement, thirty-three sidebands | 20568352 | 2200512 | 60533 |
| Fock Gaussian constructor family | 5146392 | 1910216 | 39945 |
| Phase-space, spin, and Pauli family | 1762560 | 782768 | 15165 |
| Static and moving two-level basis | 17375678 | 8518976 | 218861 |
| Static and timed frame composition | 2615642 | 1163312 | 24899 |
| Thirty-three-sideband exact phase pipeline | 8371124 | 4235168 | 82222 |

These are prototype-comparison baselines, not new CI thresholds. A migration decision uses the
plan's complete-workflow 10% time and 20% allocation gates, repeated under the same session and
fixture when comparing two implementations.
