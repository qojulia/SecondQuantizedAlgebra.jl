# Wave 1 exact-response solver prototype report

## Scope and environment

This is a disposable spike for the exact closed-adjoint plan. It changes no production source.
It was run from detached commit `4a26c86` with Julia 1.12.6 and the root project in
`/tmp/sqa-agent-solver-4a26c868`.

The prototype code is in:

- `prototype/exact_response_solver.jl`
- `prototype/run_exact_response_solver.jl`

All reported timings are warm medians of seven runs from one Julia MCP session. They are useful
for relative design decisions, not stable benchmark thresholds.

## Algorithm and signs

For the characteristic polynomial

```math
p(\lambda)=\lambda^n+c_1\lambda^{n-1}+\cdots+c_n,
```

the prototype uses the Faddeev-LeVerrier recurrence

```math
B_0=I,
\qquad
c_k=-\frac{\operatorname{tr}(A B_{k-1})}{k},
\qquad
B_k=A B_{k-1}+c_k I.
```

It then constructs

```math
\det(A)=(-1)^n c_n,
\qquad
\operatorname{adj}(A)=(-1)^{n+1}B_{n-1},
\qquad
x=\frac{\operatorname{adj}(A)f}{\det(A)}.
```

Only division by the exact integer `k` occurs in the recurrence. There are no symbolic pivots.
For a symbolic 2 by 2 matrix the result is exactly

```math
\det(A)=ad-bc,
\qquad
\operatorname{adj}(A)=
\begin{pmatrix}d&-b\\-c&a\end{pmatrix},
```

and both the numerator and final quotient are structurally equal to the existing explicit 2 by 2
formula.

## Correctness evidence

- Symbolic 2 by 2 Faddeev-LeVerrier and explicit formulas are structurally equal.
- Independent floating-point `A \ f` oracles agree for a 2 by 2 and a dense 4 by 4 case, with
  solution errors `1.12e-16` and `1.14e-16`.
- The division-free certificate
  `A * (adj(A) * f) == det(A) * f` is structurally zero for every symbolic fixture.
- The residual formed after division, `A*x-f`, is not structurally zero even for the generic
  symbolic 2 by 2 case. This is a coefficient simplification limitation, not a solver error.
- Solving a block-diagonal 8 by 8 matrix as two 4 by 4 blocks produces quotients that are not
  structurally equal to the whole-block quotients, but cross-multiplication proves their exact
  rational equivalence.

Production validation should therefore retain a division-free determinant/numerator certificate.
It should not require simplification of a completed rational residual.

## Expression size and warm cost

| Fixture | Determinant terms/chars | Adjugate chars | Solution chars | Solve median | Allocations | Certificate median/allocations |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| symbolic 2 by 2, FL | 2 / 9 | 6 | 51 | 177 us | 34,528 B | 5.9 us / 4,672 B |
| symbolic 2 by 2, explicit | 2 / 9 | n/a | 51 | 154 us | 29,504 B | same certificate applies |
| coupled tridiagonal 4 by 4 | 5 / 59 | 214 | 469 | 1.02 ms | 231,504 B | 48.7 us / 19,472 B |
| coupled tridiagonal 8 by 8 | 34 / 651 | 6,866 | 10,288 | 29.4 ms | 7,543,864 B | 473 us / 288,816 B |
| whole block-diagonal 8 by 8 | 14 / 296 | 3,956 | 5,228 | 8.76 ms | 3,199,872 B | 129 us / 132,256 B |
| same matrix split into two 4 by 4 blocks | n/a | n/a | 938 | 1.34 ms | 463,008 B | per-block certificates |

The 8 by 8 tridiagonal result is still tractable for one solve, but repeating it for many harmonic
components would be costly. Structural block decomposition is not optional: on this fixture it is
about 6.5 times faster, uses about 6.9 times fewer allocated bytes, and emits a solution about 5.6
times smaller.

The explicit 2 by 2 path remains worth retaining. Its advantage is modest in this prototype, but
it is simpler, avoids the general recurrence, and is already the common one-mode path.

## Singular systems and bounded-response semantics

The truthful initial boundary is structural, not symbolic-rank based:

| Case | Result |
| --- | --- |
| Nonsingular determinant | Return `adj(A)f/det(A)`, valid away from `det(A)=0` |
| Disconnected zero 1 by 1 block with zero forcing | Return zero for that coordinate |
| Disconnected zero 1 by 1 block with nonzero forcing | Refuse as inconsistent; for a timed harmonic this is resonance with no bounded response |
| Connected structurally singular block, even with demonstrably compatible numeric forcing | Refuse the unresolved block |
| Symbolic determinant not structurally zero | Return the ordinary quotient; do not branch on possible symbolic rank loss |

A coordinate-pivoted particular solution for a connected rank-deficient block would introduce
pivot-domain conditions and make the claimed "zero null-space contribution" basis dependent.
There is no canonical exact choice in the current representation. The first implementation should
only special-case zero-forced structural null blocks exposed by coupling-graph decomposition and
refuse every connected singular block. A future solver may accept additional named structures only
with an explicit, documented particular-solution convention.

## Why Bareiss remains unsuitable

The coefficient probe

```math
\frac{(ad-bc)g}{ad-bc}
```

is stored as `(a*d*g - b*c*g) / (a*d - b*c)`, has 29 printed characters, and is not structurally
equal to `g`. Bareiss would repeatedly rely on precisely these composite exact divisions. Under the
current `CNum` domain it would accumulate rational-expression noise and artificial pivot-domain
conditions instead of maintaining polynomial intermediates. This confirms the plan's rejection of
Bareiss for now.

## Storage, inference, and JET feasibility

- Persistent prototype storage is `Matrix{CNum}`, `Vector{CNum}`, and the concrete immutable
  `FLResult(det::CNum, adjugate::Matrix{CNum})`.
- There are no `Any` containers, closures in stored records, or mode-count type parameters.
- `@inferred` passes for `faddeev_leverrier`, `solve_adjugate`, and the explicit 2 by 2 solver.
- `JET.report_call` and `JET.report_opt`, restricted to the prototype module, report zero findings
  for the recurrence kernel. An unrestricted analysis descends into known Symbolics internals and
  produces hundreds of dependency reports, so production JET should follow the repository's
  existing `target_modules` policy.

The current implementation allocates fresh dense matrices per recurrence step. That is acceptable
for the design spike. A production implementation can reuse two dense buffers after whole-workflow
benchmarks, without changing the API or result types.

## Recommendation

**Adopt with a strict boundary.**

1. Keep scalar and explicit 2 by 2 solvers as specialized fast paths.
2. Decompose the structural coupling graph first and invoke Faddeev-LeVerrier independently on
   each connected block of size greater than two.
3. Return one determinant domain per connected block and retain the division-free
   `A*numerator == det*forcing` certificate during construction or tests.
4. Initially cap generic connected symbolic blocks at a measured size or expression-work budget;
   8 by 8 is feasible but already 29 ms and 7.5 MB for one modest tridiagonal fixture.
5. Accept a structural null block only when it is an isolated zero coordinate with zero forcing,
   choosing zero displacement. Refuse connected singular blocks and nonzero resonant forcing.
6. Keep the solver private to the transformation subsystem. No second concrete coefficient-level
   consumer was found in this spike, so promotion to shared exact linear algebra is premature.
7. Do not implement Bareiss until `CNum` supports exact cancellation of composite polynomial
   divisions.

## Files safe to discard

Everything created by this spike is disposable:

- `prototype/exact_response_solver.jl`
- `prototype/run_exact_response_solver.jl`
- `prototype_report.md`

No production file was edited.

