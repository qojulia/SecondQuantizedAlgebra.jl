# Wave 1B prototype report: multimode Gaussian conventions and gauges

## Question

Which multimode ladder ordering should a closed-adjoint implementation use, does the general
quadratic-affine scanner agree with the package algebra and quadrature representation, and can a
complete timed gauge be reconstructed from the affine map `M,d` and its derivatives alone?

This is a disposable spike. It changes no production source.

## Competing designs

### Ladder ordering

1. Site-local interleaving:
   `zI = (a1,a1',a2,a2',...)`.
2. Split Nambu ordering:
   `zN = (a1,a2,...,a1',a2',...)`.

For two modes their exact commutator matrices and adjoint permutations are

```text
CI = [ 0  1  0  0;      adjoint = (2,1,4,3)
      -1  0  0  0;
       0  0  0  1;
       0  0 -1  0 ]

CN = [ 0  0  1  0;      adjoint = (3,4,1,2)
       0  0  0  1;
      -1  0  0  0;
       0 -1  0  0 ]
```

Split order exposes the familiar particle/hole blocks. Interleaved order matches the package's
existing `sort(Op[...])`, `generators(U)`, site coverage, and `fundamental_operators`-style local
layout. Both can construct adjoint rules from only the annihilation half.

### Gauge reconstruction

1. Reconstruct the entire gauge from the affine path.
2. Reconstruct only its quadratic and linear parts and require an explicit scalar lift/cocycle.

The second design is required. The first is mathematically underdetermined, not merely difficult
to simplify.

## Fixtures and refusal cases

The 12-term two-mode Hamiltonian contains all of:

- unequal number terms;
- complex hopping `g*a1'*a2 + conj(g)*a2'*a1`;
- single-mode anomalous terms `lambda*a1'^2/2 + h.c.`;
- two-mode anomalous terms `mu*a1'*a2' + h.c.`;
- independent complex drives on both modes.

The quadrature fixture contains the equivalent `x^2`, `p^2`, symmetrized `xp`, cross-mode
`x1*x2`, `p1*p2`, `x1*p2`, `p1*x2`, and linear-force terms. Timed fixtures cover Fock and
quadrature displacement, rotation, one-mode squeeze, beamsplitter, two-mode squeeze, inverse, and
composition.

The scratch scanner deterministically refuses an operator outside the supplied basis and a cubic
Hamiltonian. Gauge reconstruction refuses a singular/noncanonical commutator basis.

## Independent oracle

The scan uses the coefficient identity

```math
[z_i z_j,z_k]=z_i C_{jk}+C_{ik}z_j
```

and is compared exactly with independent calls to package `commutator`. In split order the first
two hand-derived action rows are

```math
i[H,a_1]=-i(\omega_1 a_1+g a_2+\lambda a_1^\dagger+\mu a_2^\dagger+\eta),
```

```math
i[H,a_2]=-i(g^* a_1+\omega_2a_2+\mu a_1^\dagger+\xi).
```

An independent numerical coordinate matrix uses

```math
x=(a+a^\dagger)/\sqrt 2,
\qquad
p=i(a^\dagger-a)/\sqrt 2,
\qquad
a=(x+ip)/\sqrt 2.
```

At the recorded nontrivial complex parameter point, `Aq = T*Aa*T^-1` and `bq = T*ba` agree to
`2e-14`.

The exact ladder-to-quadrature term table used by the fixture is:

| Ladder term | Quadrature term |
| --- | --- |
| `omega*a'*a` | `omega*(x^2+p^2-1)/2` |
| `lambda*a'^2/2 + h.c.` | `real(lambda)*(x^2-p^2)/2 + imag(lambda)*(xp+px)/2` |
| `eta*a' + h.c.` | `sqrt(2)*(real(eta)*x + imag(eta)*p)` |
| `g*a1'*a2 + h.c.` | `real(g)*(x1*x2+p1*p2) + imag(g)*(p1*x2-x1*p2)` |
| `mu*a1'*a2' + h.c.` | `real(mu)*(x1*x2-p1*p2) + imag(mu)*(x1*p2+p1*x2)` |

The two-mode number offset is `-(omega1+omega2)/2`. In the package's canonical `x`-before-`p`
storage, the symmetrized local cross term additionally exposes `-im*imag(lambda)/2`; that scalar
is paired with the `xp` term and does not make the Hermitian Hamiltonian non-Hermitian.

## Correctness result

`prototype/run_gaussian_conventions.jl` passes 51 tests after the two refusal checks are included
(the immediately preceding captured run was 49/49 before adding those two direct refusal tests).
It establishes:

- exact equality of one-pass and direct-commutator affine actions in both orderings;
- the two hand-derived action rows and complex forcing;
- permutation equivalence of both ladder orderings;
- numerical ladder/quadrature equivalence;
- identical rule dictionaries materialized from either ordering;
- the timed velocity equation for every named Gaussian constructor:
  `d(U'*z*U)/dt = -im*[gauge,U'*z*U]`;
- the same equation for inverse and composed transforms;
- equality of the reconstructed and existing quadratic/linear gauges.

With package convention `G = im*(dU'/dt)*U`, and with the stored gauge expressed in the original
generator labels, the correct body-coordinate reconstruction is

```math
M^{-1}\dot M=iCK,
\qquad
M^{-1}\dot d=iC\ell,
```

for `G = z'Kz/2 + ell'z + s`. The initially tempting spatial velocity
`dot(M)*M^-1` gives wrong moving-squeeze and quadrature signs.

The remaining exact scalar residuals after reconstructing `K,ell` were:

| Constructor | Required scalar lift |
| --- | --- |
| Fock displacement `alpha=u*t+im*v*t^2` | `u*v*t^2` |
| Quadrature displacement `dx=u*t, dp=v*t^2` | `u*v*t^2/2` |
| Fock rotation `theta*t` | `theta/2` |
| Quadrature rotation | `0` |
| One-mode squeeze | `0` |
| Beamsplitter | `0` |
| Two-mode squeeze | `0` |
| Quadrature squeeze | `0` |

The displacement values are the Weyl cocycle `(dx*dot(dp)-dp*dot(dx))/2` in the appropriate
shifts. Most importantly, Fock `exp(-im*theta*n)` and quadrature
`exp(-im*theta*(x^2+p^2)/2)` induce the same canonical rotation but their gauges differ by
`dot(theta)/2`. Thus `M,d,dot(M),dot(d)` cannot determine the complete gauge. An explicit unitary
lift/phase convention is unavoidable.

## Inference and JET result

`@inferred` passes for the commutator matrix, direct action, scanned action, and rule materializer;
their return storage is concrete `Matrix{Coeff}`, `Vector{Coeff}`, and `Dict{Op,QAdd}`.

A broad `JET.report_call` on the scratch `commutator_matrix(::Vector{Op})` entered unrestricted
Symbolics internals and produced 599 possible upstream errors, so it is not a useful acceptance
signal for this disposable generic signature. Production JET checks should target the bounded
coefficient scanner and fixed materializer methods, as the existing suite does, rather than a
generic `Vector{Op}` scratch entry point.

## Runtime, allocations, and expression size

Warm medians on Julia 1.12.6 from 21 samples:

| Operation | Time | Allocations | Bytes |
| --- | ---: | ---: | ---: |
| Materialize 33-mode interleaved rules | 51.1 us | 929 | 128184 |
| Materialize 33-mode split rules | 53.4 us | 929 | 128184 |
| Four direct package commutators | 38.0 us | 868 | 58480 |
| Naive scratch single Hamiltonian scan | 211 us | 2777 | 139120 |

The 12 Hamiltonian terms produce 18 nonzero action/forcing coefficients. Ordering has no material
storage or allocation effect. The scratch scanner proves the algebraic identity but is slower
because it constructs intermediate `QAdd`s and repeatedly searches the basis. Do not port it.
The production scanner must accumulate `CNum` entries directly with a precomputed `Op => Int`
lookup (or equivalent site arithmetic), as the design plan already requires.

## Production concepts and lines introduced

Production lines introduced: zero.

Scratch code is 474 lines, deliberately including fixtures, independent oracles, diagnostics, and
measurement scaffolding. Only these concepts should be ported:

1. persistent site-local interleaved generator order;
2. an optional explicit permutation/view into split Nambu blocks for Gaussian solver strategies;
3. direct quadratic-affine coefficient scanning from the commutator matrix;
4. non-central gauge reconstruction from `M^-1*dot(M)` and `M^-1*dot(d)`;
5. a strategy-supplied scalar lift/cocycle.

## Failure modes

- A complete gauge inferred from affine data silently chooses a unitary phase and cannot reproduce
  both existing rotation conventions.
- Using `dot(M)*M^-1` reconstructs a spatial generator, whereas the stored package gauge requires
  the body generator `M^-1*dot(M)`.
- Generic `LinearAlgebra` multiplication on `Matrix{Coeff}` fails because `zero(::Coeff)` is not
  defined; production code needs the package's explicit `CNum` operations.
- Plain Symbolics simplification did not eliminate one hyperbolic residual; applying the
  transform's certified `cosh^2-sinh^2=1` relation did.
- A naive single scan can allocate more than direct commutators unless it accumulates directly.
- Persistent split ordering would require needless permutation at every existing site/rule API
  boundary; persistent interleaving merely needs a split view for algorithms that use BdG blocks.

## Recommendation

### Adopt

- Use site-local interleaving as the persistent `AdjointBasis` generator order because it matches
  current Hilbert-space and `UnitaryTransform` behavior.
- Let Gaussian solvers request a split-Nambu permutation/view; do not store a second basis.
- Share quadratic and linear timed-gauge reconstruction using body velocities.
- Make the scalar lift an explicit output of each exact strategy. Weyl displacement, named
  rotations, squeezes, and exact exponentiation strategies must state their lift convention.

### Reject

- Reject a universal complete-gauge builder whose inputs are only `M,d,dot(M),dot(d)`.
- Reject the scratch scanner implementation and generic matrix arithmetic on `Coeff`.
- Reject changing public generator ordering merely to expose Nambu blocks.

### Investigate in the materializer wave

- The smallest concrete representation of the scalar lift: preferably a computed `CNum` supplied
  alongside the certified map, not a callback or metadata object retained in `UnitaryTransform`.
- Whether homogeneous exact-flow strategies standardize on the metaplectic/Weyl lift and named
  Fock rotation adds its known `dot(theta)/2` phase correction at construction.

## Files safe to discard

- `prototype/gaussian_conventions.jl`
- `prototype/run_gaussian_conventions.jl`

The entire isolated worktree may be removed after the coordinator copies this report. No
production file or main worktree was edited.

