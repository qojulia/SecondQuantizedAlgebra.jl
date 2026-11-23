# Docs cleanup checklist (`docs-cleanup`)

Line numbers refer to the state of the branch after the typo pass.

## Done

- [x] Typos and wording in the examples and manual
- [x] British → American spelling (43 replacements)

## 1. Physics and numerical errors

- [x] **1.1** [dicke.jl:101-102](examples/dicke.jl#L101-L102): reverse the Holstein–Primakoff mapping to S₊ ≈ √(2S) b†, S₋ ≈ √(2S) b (expansion around S_z = −S, noted against Wikipedia's S_z = +S form)
- [x] **1.2** [dicke.jl:180](examples/dicke.jl#L180): change `\mp` to `\pm` so that ε₋ is the lower branch
- [ ] **1.3** [dicke.jl:205-206](examples/dicke.jl#L205-L206): point the `Squeeze` cross-link to the parametric oscillator example instead of optomechanics
- [x] **1.4** [dicke.jl:302-304](examples/dicke.jl#L302-L304): soften "below λ_c they coincide exactly"
- [x] **1.5** [optomechanics.jl:244-245](examples/optomechanics.jl#L244-L245): change ε(1), ε(2) to ≈ 1.1832 and 1.3416
- [x] **1.6** [optomechanics.jl:75-83](examples/optomechanics.jl#L75-L83): fix the "ã = a" justification (only ã†ã = a†a holds), and change "three substitutions" to two
- [x] **1.7** [optomechanics.jl:192-193, 246-250](examples/optomechanics.jl#L192): describe it as the single-mode analogue of the Bogoliubov example, and fix the phonon-counting logic
- [x] **1.8** [schrieffer_wolff.jl:131-139](examples/schrieffer_wolff.jl#L131-L139): fix χ_eff for the excited branch (the extra |n,e⟩ check was left out on purpose)
- [x] **1.9** [parametric_oscillator.jl:165-169](examples/parametric_oscillator.jl#L165-L169): define M and det, or derive the result from ε = √(δ²−κ²)
- [ ] **1.11** [two_atom.jl](examples/two_atom.jl):
  - [ ] **1.11a** lines 71-72: say that distinct sites commute exactly
  - [ ] **1.11b** lines 191-194: describe the actual commutator output, or apply `expand_completeness`
  - [ ] **1.11c** line 240: the single-atom splitting is 2g, not g
  - [ ] **1.11d** lines 115-118, 150: rename the "singlet" states
  - [ ] **1.11e** line 169: rename the clashing S₋
  - [ ] **1.11f** lines 85-88: fix the 1/√2 normalization
  - [ ] **1.11g** line 102: define Γ
- [x] **1.12** [truncated_wigner.jl:257](examples/truncated_wigner.jl#L257): the plot title should say t = 400; line 259 is a stale Makie comment that splits the code block
- [x] **1.13** [bogoliubov.jl:193-196](examples/bogoliubov.jl#L193-L196): define u, v before line 135
- [ ] **1.14** [spaces_and_operators.md:225](docs/src/manual/spaces_and_operators.md#L225): the LaTeX example shows counter-rotating terms; it probably needs `Transition(h, 1, 2, 2)`
- [ ] **1.15** [algebra.md:109](docs/src/manual/algebra.md#L109): the Pauli commutator is 2i times the third component
- [ ] **1.16** [README.md:23](README.md#L23): remove the stale "Extensible for custom operator types"

## 2. Integrate the sections appended by recent PRs

- [ ] **2.1** [coefficients.md:190-202](docs/src/manual/coefficients.md#L190-L202): delete the orphan `pr305` paragraph, which repeats line 46
- [ ] **2.2** [coefficients.md:204-238](docs/src/manual/coefficients.md#L204-L238), harmonic projection:
  - [ ] **2.2a** rename the `pr311` block
  - [ ] **2.2b** drop the unused `ω t` and the unneeded `import`
  - [ ] **2.2c** move it after unitary transformations, or link there
- [ ] **2.3** [coefficients.md:240-250](docs/src/manual/coefficients.md#L240-L250): move `map_coefficients` and `coefficient` to "Taking expressions apart", with examples
- [ ] **2.4** [many_body.md:257-260](docs/src/manual/many_body.md#L257-L260): move `get_sum_indices` to "Inspecting indices"
- [ ] **2.5** [unitary_transformations.md:194-202](docs/src/manual/unitary_transformations.md#L194-L202): fold `displacement` into the Displacements section, with an example
- [ ] **2.6** [open_systems.md:334-345](docs/src/manual/open_systems.md#L334-L345): move "Secular projection" after "Rotating frame", and drop the unneeded `import`
- [ ] **2.7** [averages_numerics.md:145-205](docs/src/manual/averages_numerics.md#L145-L205): rewrite "Bases and levels" with executed `@example` blocks, removing the duplication with line 129 and the QuantumToolbox-callback detour
- [ ] **2.8** [extending.md:134-140](docs/src/dev/extending.md#L134-L140): move the accumulation paragraph up next to its table row

## 3. Formatting consistency

- [ ] **3.1** Join hard-wrapped paragraphs into one line per paragraph, as on the rest of the manual
- [ ] **3.2** Remove the double blank lines (coefficients, many_body, unitary_transformations, open_systems)
- [ ] **3.3** Decide on `CollapsedDocStrings`: set it on every page or on none (currently only open_systems)
- [ ] **3.4** Use hyphens or en dashes consistently in names (Jaynes–Cummings, Holstein–Primakoff, Bose–Hubbard, …)

## 4. Examples as a set

- [ ] **4.1** Align the nine `# # Title` headings with [pages.jl](docs/pages.jl)
- [ ] **4.2** Expand [tavis_cummings.jl](examples/tavis_cummings.jl) with an intro, sections, output, a plot and a closing paragraph, or merge it with two_atom's N-atom part
- [ ] **4.3** Add closing paragraphs to schrieffer_wolff and differentiable_control
- [ ] **4.4** two_atom: use `to_numeric` instead of building the operators by hand three times, and fix the "no bookkeeping" claims (lines 131, 151, 322)
- [ ] **4.5** Add citations through `refs.bib`: Schrieffer–Wolff, Emary–Brandes, Holstein–Primakoff, Aspelmeyer RMP, Polkovnikov
- [ ] **4.6** Unify plot sizes, label style (LaTeXStrings), title case, and `# hide` vs `#hide`
- [ ] **4.7** Use symbol levels (`:g, :e`) everywhere, including schrieffer_wolff
- [ ] **4.8** Tone down "story", "The star of this example", and the heavy bolding
- [ ] **4.9** differentiable_control: remove the unused `using Symbolics`, deduplicate the plotting code, and fix the "near unity" comment

## 5. Small items

- [ ] **5.1** [comparison.md](docs/src/comparison.md): rerun or relabel the results (they say SQA 0.10.0; the package is now 0.15.0), and change the sneg link to https
- [ ] **5.2** [many_body.md:183](docs/src/manual/many_body.md#L183): add the machine and version to the timing table, or drop the absolute numbers
- [ ] **5.3** Delete the unused [docs/src/assets/eager_vs_lazy.svg](docs/src/assets/eager_vs_lazy.svg)
- [ ] **5.4** [index.md](docs/src/index.md):
  - [ ] **5.4a** replace the hot-linked SciML icon
  - [ ] **5.4b** move the citation next to QuantumCumulants
  - [ ] **5.4c** rename `:atoms` to `:atom` (also in the README)
- [ ] **5.5** README: use consistent heading levels, put `pkg>` in a `julia-repl` block, and drop the redundant `simplify`
- [ ] **5.6** [getting_started.md:7](docs/src/manual/getting_started.md#L7): remove "We use units with ℏ=1."
- [ ] **5.7** [operator_functions.md:135](docs/src/manual/operator_functions.md#L135): link the Schrieffer–Wolff example

## 6. Before the PR

- [ ] **6.1** `make docs` builds cleanly (examples run, doctests pass, no broken `@ref`)
- [ ] **6.2** Check the changed plots and the formulas in the rendered pages
- [ ] **6.3** Add a changelog entry, if docs changes need one under your CI check
