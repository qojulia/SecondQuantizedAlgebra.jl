# JOSS readiness audit

Checked on 17 October 2026 against the [JOSS review checklist](https://joss.readthedocs.io/en/latest/review_checklist.html)
and [review criteria](https://joss.readthedocs.io/en/latest/review_criteria.html).
The software source is commit `c4c2dbae2bed34dd2cbb81dc055e3285d4a566b1`
on `joss-release`, inherited from `type-stability`; manuscript and onboarding
changes are still uncommitted. This is a preparation audit, not an editorial
acceptance decision.

## Documentation and functionality

| Criterion | Evidence and status |
| :-- | :-- |
| Open repository and license | Public repository with a complete MIT `LICENSE`. |
| Installation | Installed registered SQA v0.14.0 and QuantumOpticsBase in a fresh temporary project on Julia 1.13.1. The README commutator matched its analytic expression and numerical conversion produced the expected 10-by-10 matrix. This used the registered package rather than a development path. |
| Requirements | `Project.toml` declares Julia 1.10 and dependencies. README now states the minimum Julia version, automatic dependency installation, and the optional numerical backends. |
| Need and audience | README and paper identify symbolic derivations of quantum-optical and many-body models and reuse by research-software developers. |
| Manual and API | Nine manual chapters cover spaces/operators, algebra, coefficients, indexed systems, formal functions, transformations, averages/numerics, and open systems. The fresh full build completed, including doctests, `checkdocs = :exports`, cross-references, and HTML rendering. |
| Research examples | All nine existing Literate examples completed in the full build. Their source files were not changed. No further example work is planned at the maintainer's request. |
| Published site | Stable and development home pages, the tutorial, numerics, comparison, and contributor guide return HTTP 200. All 26 navigation-linked development pages return HTTP 200 and have no missing same-page fragment targets. Browser-based visual inspection was unavailable; these are HTML and HTTP checks. |
| Automated tests | A fresh local test subprocess reports 5,039 passing assertions across 73 test files on Julia 1.13.1. The MCP parent timed out, but the subprocess completed with `SUCCESS`; exact-source CI independently passed on minimum and current Julia. JET is excluded by the ordinary runner and checked separately. |
| JET | Fresh local analysis and exact-source [JET CI](https://github.com/qojulia/SecondQuantizedAlgebra.jl/actions/runs/39025610219) passed: five package/extension correctness checks plus 588 passing per-feature checks and two expected broken checks for known backend dispatch. This is not a claim that every possible call is inferable. |
| Contribution, bug reports, support | Root `CONTRIBUTING.md` links the developer guide and explicitly directs all three activities to GitHub. README now links this guide and explains `make test`, `make jet`, and `make docs`. |

Exact-source CI also passed [tests](https://github.com/qojulia/SecondQuantizedAlgebra.jl/actions/runs/39025610173)
and [documentation](https://github.com/qojulia/SecondQuantizedAlgebra.jl/actions/runs/39025610002).
Those runs validate the underlying software source, not the unpushed manuscript
or README changes. The documentation log confirms Literate execution, doctests,
and HTML rendering. Its only error messages are post-build cache-cleanup HTTP
403s; they do not indicate documentation failures.

## Remaining release and submission work

1. **Complete the two manuscript author actions.** Funding applicability and the
   complete AI-use history cannot be established from public affiliations or
   this session. Authors must review the assisted material and confirm the
   corresponding disclosure. Approve authorship, affiliation details, and the
   final manuscript before submission. ORCIDs are optional.
2. **Select the reviewed software version and resolve its downstream story.**
   QuantumCumulants currently permits SQA 0.12/0.13 and QuantumInputOutput permits
   0.10/0.11/0.12. Their green [downstream jobs](https://github.com/qojulia/SecondQuantizedAlgebra.jl/actions/runs/39025610110)
   skipped on resolver incompatibilities; FloquetExpansions actually ran and
   passed. Existing research use remains valid, but running those projects
   against this release has not been demonstrated. Updating other repositories
   is separate work. JOSS does not require every downstream package to support
   the newest release, but the paper and demonstrations must accurately identify
   which version is used.
3. **Assess the performance claims before relying on them.** The documentation's
   comparison numbers are explicitly a July 2026 snapshot for SQA 0.10.0 and
   named competitor versions; they were not rerun for this audit. The
   [benchmark workflow for the preparation commit](https://github.com/qojulia/SecondQuantizedAlgebra.jl/actions/runs/39025610129)
   failed its configured performance threshold. This is not a failing algebra
   test or proof of a source regression; it needs controlled measurement before
   an explanation. The paper makes no quantitative speedup claim. Keep the
   historical scope clear, or revalidate the comparison for the reviewed release.
4. **Publish the prepared material.** User commits and pushes are still needed,
   followed by submitting the actual JOSS form.
   No submission or release has been made.
5. **Archive the accepted version at the appropriate stage.** After review,
   create the agreed tag and software archive DOI and add its reference. This
   acceptance-stage requirement is distinct from the current manuscript tasks.

The repository has sustained development since June 2025, accessible early
PRs and releases, and three concrete research integrations. The scientific
provenance and research-use evidence are collected in [the preparation notes](README.md#research-use-evidence).

## Fresh local build results

The full documentation build completed on Julia 1.13.1 with draft mode off,
two Literate workers, doctests, exported-API checks, and link checking enabled.
All nine existing examples completed and 26 HTML pages were generated.
There were no documentation errors; precompilation emitted stale-pidfile and
upstream ModelingToolkitBase `@nospecialize` warnings. A static scan of the
built HTML found only the expected `siteinfo.js` and `versions.js` references
that Documenter supplies during deployment, rather than in a local build.

The initial docs MCP connection timed out while Julia continued building; the
completed log contains `JOSS AUDIT DOCS COMPLETE`. A separate test-environment
MCP session could not start because its helper expected `TestEnv`. Standalone
JET was therefore run in the already-resolved `Pkg.test` temporary environment
against the same checkout. It completed with `JOSS AUDIT JET COMPLETE`, five
passing package/extension checks, and 588 passing per-feature checks with two
expected broken checks. No package dependencies were added for this fallback.

Local transient logs are `/tmp/sqa-joss-docs.log`, `/tmp/sqa-joss-jet.log`,
`/tmp/sqa-joss-tests.log`, and `/tmp/sqa-joss-install.log`. Julia's generated
absolute development path was removed from the tracked docs project afterward.
No Julia source, example, or test files were changed.
