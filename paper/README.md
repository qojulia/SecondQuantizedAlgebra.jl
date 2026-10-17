# JOSS submission preparation

The manuscript is [`paper.md`](paper.md), with references in
[`paper.bib`](paper.bib). This is a draft for author review, not a submitted
application. The author order follows the maintainer's instructions:
Orjan Ameye, Christoph Hotter, David Plankensteiner.

## Submission details

| Field | Prepared value |
| :-- | :-- |
| Journal | Journal of Open Source Software |
| Title | SecondQuantizedAlgebra.jl: Reusable symbolic operator algebra for quantum optics and many-body physics |
| Repository | https://github.com/qojulia/SecondQuantizedAlgebra.jl |
| Preparation branch | `joss-release` |
| Paper path | `paper/paper.md` |
| License | MIT, in the repository's root `LICENSE` |
| Language | Julia |
| Documentation | https://qojulia.github.io/SecondQuantizedAlgebra.jl/ |
| Package version at preparation | 0.14.0; the reviewed release is still to be selected |
| Archive DOI | Pending the reviewed release and Zenodo/figshare deposit |

Select the editor and reviewer suggestions in the live submission form after
checking expertise, availability, and conflicts. Declare related publications,
including the QuantumCumulants software paper and the quantum-dot input-output
preprint cited in the manuscript. The JOSS paper describes the algebra package
and its reuse; those papers describe the predecessor framework and a downstream
science application.

## Author details and funding

Affiliations were checked against public sources on 17 October 2026:

| Author | Draft affiliation | Evidence |
| :-- | :-- | :-- |
| Orjan Ameye | Department of Physics, University of Konstanz, Germany | [Personal page](https://oameye.github.io/) and [group directory](https://www.zilberberg.uni-konstanz.de/home/index.php?id=group) |
| Christoph Hotter | Department of Physics, Massachusetts Institute of Technology, USA | Maintainer confirmation; [GitHub profile](https://github.com/ChristophHotter); [August 2026 preprint](https://arxiv.org/html/2608.28470v1) lists the MIT physics affiliation |
| David Plankensteiner | QuEra Computing Inc., USA | [GitHub profile](https://github.com/david-pl) lists Quera Computing Inc.; confirm his preferred affiliation for this work |

ORCIDs are omitted until verified against the authors' own records. Do not
substitute an identifier belonging to a coauthor with a similar search result.
Confirm that all authors agree to authorship, order, affiliations, and the
submission, including David's contribution to the inherited algebra.

Funding cannot be assigned to this software from affiliation alone. The
input-output preprint acknowledges FWF FG5, COE1, J4865 and the Carlsberg
Foundation's QCooL project. These are candidates to check with Christoph,
not established funding acknowledgements for this JOSS submission. No grants
have been copied into the manuscript without confirmation. Complete the
acknowledgements and AI disclosure before submission; the manuscript marks
both remaining author actions explicitly.

## Research-use evidence

These sources were inspected on 17 October 2026. Compatibility ranges are a
snapshot, not a request to update other repositories as part of this draft.

| Workflow | Evidence | Meaning for this submission |
| :-- | :-- | :-- |
| QuantumCumulants | [Published framework](https://doi.org/10.22331/q-2022-01-04-617); [dependency declaration](https://github.com/qojulia/QuantumCumulants.jl/blob/98702fe3534a0e7d9cf02e0c33544e62b949152d/Project.toml) | Published scientific provenance plus current reuse. Older QuantumCumulants papers are not evidence of using the later standalone package. Current declared SQA compatibility is 0.12/0.13. |
| QuantumInputOutput | [README](https://github.com/qojulia/QuantumInputOutput.jl/blob/ef58639bf01303ff170291e4290838a251bb9b6b/README.md); [dependencies](https://github.com/qojulia/QuantumInputOutput.jl/blob/ef58639bf01303ff170291e4290838a251bb9b6b/Project.toml); [SUPER preprint](https://arxiv.org/html/2608.28470v1) | Symbolic input-output modeling and a documented research application: the preprint cites QuantumInputOutput documentation and reproducible materials. Current declared SQA compatibility is 0.10/0.11/0.12. |
| FloquetExpansions | [Dependencies](https://github.com/oameye/FloquetExpansions.jl/blob/6ec610ac9b36c2300f44d57400d03fd5956e39c6/Project.toml); [research workflow](https://github.com/oameye/FloquetExpansions.jl/blob/6ec610ac9b36c2300f44d57400d03fd5956e39c6/docs/src/index.md) | Active research integration for high-frequency effective generators and micromotion, with SQA 0.14 compatibility. No publication count is claimed. |

The [downstream CI workflow](../.github/workflows/Downstream.yml) targets all
three projects and treats dependency-resolution incompatibilities as successful
skips. In the [run for the preparation source commit](https://github.com/qojulia/SecondQuantizedAlgebra.jl/actions/runs/39025610110),
FloquetExpansions actually ran and passed its tests; QuantumCumulants and
QuantumInputOutput skipped because their compatibility bounds excluded the
candidate. A green aggregate status does not establish that those two projects
work with this release. Resolve or explicitly document that release relationship
before presenting them as demonstrations of the candidate version.

The repository was created on 3 June 2025 and is currently public. Local Git
history spans June 2025 to October 2026, with commits in 13 separate months.
Recent public releases include v0.12.0 (10 September), v0.13.0 (5 October),
and v0.14.0 (15 October 2026). Creation and commit dates alone do not establish
when the repository became public. Early accessible development records include
[the initial documentation PR](https://github.com/qojulia/SecondQuantizedAlgebra.jl/pull/4)
(3 June 2025) and [v0.1.0](https://github.com/qojulia/SecondQuantizedAlgebra.jl/releases/tag/v0.1.0)
(7 June 2025). These supplement the sustained commit and release history for
the six-month public-development requirement. Some issues were transferred from
QuantumCumulants and have older creation dates; do not use those dates to date
the standalone repository.

## Before submission

- [ ] Authors approve the paper, author order, and affiliations; optionally supply verified ORCIDs.
- [ ] Complete funding acknowledgements and the AI usage disclosure, including earlier development and human validation.
- [ ] Add any further representative research papers the authors want cited; verify their relationship to the standalone package.
- [ ] Select the submission commit and version, accounting for changes inherited from the `type-stability` branch on which `joss-release` was created.
- [ ] Run installation and representative documented examples from a fresh environment for that commit.
- [ ] Run `make test`, `make jet`, and `make docs` for the release candidate; record results and resolve failures.
- [ ] Check downstream results against the candidate and resolve or explain compatibility restrictions before relying on those projects as release demonstrations.
- [x] Check public development-history evidence and that contribution/support links are accessible; see the dated audit.
- [ ] Compile and review the final PDF; remove every `AUTHOR ACTION` note and update the date.
- [ ] Have the user commit and push the branch.
- [ ] Fill in the [JOSS submission form](https://joss.theoj.org/papers/new), declaring related publications and conflicts.

## Building the paper

For a local build from the repository root:

```sh
podman run --rm --userns=keep-id \
  --volume "$PWD/paper:/data" \
  --env JOURNAL=joss \
  docker.io/openjournals/inara:latest
```

Docker can use the equivalent command from the
[JOSS build instructions](https://joss.readthedocs.io/en/latest/paper.html#checking-that-your-paper-compiles).
Generated PDFs and JATS files are ignored by Git. The bibliography's article
metadata was checked against Crossref; repository references identify the
software directly.

Preparation checks on 17 October 2026: the draft compiled successfully with
`openjournals/inara` (image ID `341ae8f00b95`), and all four PDF pages were
rendered and visually inspected. The manuscript has approximately 1,080 words;
all required sections and ten citation keys were checked. The preview is also saved as
`output/pdf/SecondQuantizedAlgebra-JOSS-draft.pdf` for local review. Standard
JOSS draft DOI, editor, and submission-date placeholders are supplied by the
renderer and will be replaced during editorial processing. The subsequent [readiness audit](audit.md) records installation,
tests, JET, documentation, and downstream checks. The preparation changes do not
alter Julia source or tests.

## During and after review

JOSS review is public and iterative. Authors should normally respond within
two weeks and complete requested changes within four to six weeks, or agree
a different schedule with the editor. Under the current
[AI policy](https://joss.readthedocs.io/en/latest/submitting.html#ai-usage-policy),
AI may assist with software, documentation, and manuscript preparation with
disclosure, but may not write conversational interactions with editors or
reviewers except for translation. Authors must handle that correspondence.

After successful review, make the agreed tagged release, archive the matching
source with Zenodo or figshare, and provide its DOI and version in the review
thread. The software archive DOI and the JOSS article DOI are distinct.
Add the archive reference to the paper and update software citation metadata
once the actual identifiers exist. Commits, pushes, tags, deposits, and the
submission itself remain user actions; repository policy prohibits this agent
from committing or pushing.

Requirements checked against the current [submission guidelines](https://joss.readthedocs.io/en/latest/submitting.html),
[paper format](https://joss.readthedocs.io/en/latest/paper.html), and
[review checklist](https://joss.readthedocs.io/en/latest/review_checklist.html).
