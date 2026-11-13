# Julia quantum algebra comparison for JOSS

Checked 2026-10-26 against primary sources and the live local SQA checkout. No new timings or paper changes. This report distinguishes current implementations from what another package could implement.

## Conclusion

SQA's strongest case is a particular combination: an eagerly reduced, indexed operator-polynomial algebra for quantum dynamics; expectation values represented as Symbolics scalars for existing QuantumCumulants/ModelingToolkit workflows; direct numerical conversion; and measured efficiency for substantial analytical calculations. Neither Julia, Symbolics compatibility, numerical translation, nor ecosystem ambitions are unique by themselves.

Against QuantumAlgebra, emphasize the established QuantumCumulants integration and the representation boundary between operator algebra and Symbolics equation systems, alongside demonstrated workload-specific performance. Against QuantumSymbolics, emphasize the specialized canonical polynomial and indexed-sum machinery for many-body/dynamical calculations. QuantumSymbolics already has a substantial ecosystem role and more numerical representation families.

## Versions and scope

- Local SQA is 0.14.0, with uncommitted work. Evidence below describes this checkout, not every registered release.
- Latest registered/tagged QuantumAlgebra is **1.7.0**. Verified in [General Versions.toml](https://github.com/JuliaRegistries/General/blob/master/Q/QuantumAlgebra/Versions.toml) and remote tags. The benchmark's 1.7.0 is therefore current registered. Upstream main inspected at `c7d46fcef5019c9a439abd1dcd0f1fa9989819ea` has Project version 1.7.1 but no corresponding published tag/General entry. Its release notes describe further canonicalization and coefficient optimizations. Qudits, spins, and so(n) are already in 1.7.0; u(n) and qudit construction from a matrix are explicitly unreleased. [Release notes](https://github.com/jfeist/QuantumAlgebra.jl/blob/c7d46fcef5019c9a439abd1dcd0f1fa9989819ea/docs/src/release_notes.md#L3-L47)
- Latest registered/tagged QuantumSymbolics is **0.4.20**. Inspected tagged Project and main at `434c8c64fc5eeb95cf1e2e79ead5f26f3eed8c33`. Its rendered docs can lag source: the source includes QuantumToolbox and Gabs extensions although introductory docs discuss chiefly QuantumOptics and QuantumClifford. [Tagged Project](https://github.com/QuantumSavory/QuantumSymbolics.jl/blob/v0.4.20/Project.toml)

## Comparison matrix

| Aspect | SQA local 0.14 | QuantumAlgebra 1.7 | QuantumSymbolics 0.4.20/current source |
| --- | --- | --- | --- |
| Principal scope | Canonical operator algebra for quantum optics/many-body dynamics | Closely overlapping operator algebra, expectation/cumulant manipulations, equation generation | General symbolic quantum objects and translation between numerical formalisms, especially quantum information |
| Built-in objects | Boson, finite-level transitions, Pauli, general spin, quadratures, collective transitions | Boson, fermion, TLS, qudit, spin, so(n); u(n) on main | Generic bras/kets/operators/superoperators; Pauli/gates/channels; Fock/coherent/squeezed/thermal states and Gaussian operations |
| Polynomial arithmetic | Distribution and local algebra reduction on construction | Normal form explicit by default, optionally automatic; specialized direct commutator | Constructors collect additive terms and flatten products; physical rewrites and distribution through explicit qsimplify/qexpand |
| Symbolic site sums | Index carries subsystem and range; exclusions and diagonal splitting; user index identity retained | Formal indexed sums, Kronecker deltas, canonical dummy-index renumbering; arbitrary index tuples | No bound-index, analytic-sum, or exclusion-constraint machinery found in current source/docs |
| Symbolics relationship | Scalar coefficients and expectation expressions enter Symbolics; time-dependent averages retain operator metadata | Symbolics scalar prefactors supported by extension; ExpVal/Corr remain in QuExpr; Julia Expr generation | Symbolics/SymbolicUtils/TermInterface hard dependencies, substitution and rewrite integration |
| Dynamics | Operator/superoperator/Lindbladian algebra, frames, GKSL extraction; moment closure delegated to QC | Built-in Lindblad Heisenberg equations and truncated expectation/cumulant systems | Generic superoperator/Kraus/channel expressions and translation; no analogous moment-equation generator found |
| Numeric backends bundled | QuantumOpticsBase, QuantumToolbox/SciMLOperators, including parameterized/time-dependent and generator conversion | No bundled operator-to-QO/QToolbox conversion found; emits Julia Expr for scalar expectations/correlators | QuantumOpticsBase, QuantumToolbox, QuantumClifford, Gabs, mixed Clifford/optics extension |
| Extensibility trade-off | New operator kind requires internal changes; external numeric backend interface documented | Internal operator-kind rules; arbitrary Number/CAS scalar coefficient types | External symbolic object subtypes, user rewrite chains, representation-specific express_nolookup methods |

The table describes inspected implementations. Missing an integration is not evidence that it is impossible or would be impractical to contribute.

## Substantive distinctions and evidence

### 1. QuantumAlgebra is a serious close alternative

It already provides formal sums, normal ordering, vacuum expectations, expectation-to-cumulant expansions, Lindblad Heisenberg drifts, hierarchy generation/truncation, and code emission. Its direct commutator avoids generating terms that cancel in AB-BA. Do not describe it as merely a commutator calculator or as disconnected from numerical workflows. [Module overview](https://github.com/jfeist/QuantumAlgebra.jl/blob/c7d46fcef5019c9a439abd1dcd0f1fa9989819ea/src/QuantumAlgebra.jl#L1-L9), [equation-system implementation](https://github.com/jfeist/QuantumAlgebra.jl/blob/c7d46fcef5019c9a439abd1dcd0f1fa9989819ea/src/eqsys.jl#L12-L36), [commutator documentation](https://github.com/jfeist/QuantumAlgebra.jl/blob/v1.7.0/README.md)

The precise distinction is the symbolic equation boundary. QA's Symbolics extension maps simplify/substitute over scalar coefficients. Its expectations and correlators are native quantum terms. Its julia_expression produces Julia AST for scalar expressions and explicitly refuses terms with bare operators. SQA's average produces ordinary SymbolicUtils/Symbolics scalar expressions; make_time_dependent gives them symbolic function metadata suitable for the established QC path. This is a difference in supplied interfaces, not an assertion that QA could never feed ModelingToolkit. [QA Symbolics extension](https://github.com/jfeist/QuantumAlgebra.jl/blob/v1.7.0/ext/QuantumAlgebraSymbolicsExt.jl), [QA code emission](https://github.com/jfeist/QuantumAlgebra.jl/blob/c7d46fcef5019c9a439abd1dcd0f1fa9989819ea/src/convert_to_expression.jl#L14-L24); SQA `src/average.jl:192-226,285-305`, `benchmark/comparison/demonstrations/interoperability.jl:33-50`.

There are also meaningful representation choices. QA represents a qudit ground projector by completeness and canonicalizes dummy sums by renumbering. SQA deliberately retains all finite-level projectors until explicit completeness expansion and keeps user-provided indices, rejecting ambiguous binder collisions. These choices help preserve model/moment identity and control some expression growth, while QA's choices identify more equivalent spellings automatically. Neither convention is uniformly superior. [QA qudit documentation](https://github.com/jfeist/QuantumAlgebra.jl/blob/c7d46fcef5019c9a439abd1dcd0f1fa9989819ea/README.md#L158-L188), [QA dummy sums](https://github.com/jfeist/QuantumAlgebra.jl/blob/c7d46fcef5019c9a439abd1dcd0f1fa9989819ea/src/sum_indices.jl#L1-L12); SQA `docs/src/dev/architecture.md:51-59,65-79`.

QA supports fermions and fixed-spin finite-representation identities. SQA's current operator kinds contain no fermionic family, and its SpinSpace describes unspecified angular momentum. Avoid claims of categorically broader algebra coverage by SQA. SQA `src/operators/op.jl:17`, `docs/src/manual/spaces_and_operators.md:123`; [QA release notes](https://github.com/jfeist/QuantumAlgebra.jl/blob/c7d46fcef5019c9a439abd1dcd0f1fa9989819ea/docs/src/release_notes.md#L88-L110).

### 2. QuantumSymbolics already is an ecosystem symbolic frontend

Its QuantumSavory role and conversion facilities are central, not incidental. It supports quantum states and gates as well as operators and channels, and can express objects in stabilizer and Gaussian formalisms in addition to conventional quantum-optics operators. It shares Symbolics, SymbolicUtils, TermInterface, and QuantumInterface with SQA. Therefore “SQA is unique because it connects symbolic and numerical Julia tools” would be false. [QuantumSavory symbolic frontend reference](https://qs.quantumsavory.org/dev/symbolics/), [conversion documentation](https://doc.quantumsavory.org/QuantumSymbolics.jl/stable/express/), [dependencies/extensions](https://github.com/QuantumSavory/QuantumSymbolics.jl/blob/v0.4.20/Project.toml).

Its algebra architecture differs. Product constructors preserve an expression tree after scalar collection/flattening. Manual rewriting covers Pauli identities, simple bosonic commutators, state actions, and selected Gaussian identities; qexpand distributes products and expands commutators. I found no general bosonic operator-word normal-order algorithm or bound-site-sum/inequality system. A custom implementation remains possible. SQA supplies those facilities as its core arithmetic, with physical subsystem/range data and automatic diagonal splitting. [QS multiplication constructor](https://github.com/QuantumSavory/QuantumSymbolics.jl/blob/434c8c64fc5eeb95cf1e2e79ead5f26f3eed8c33/src/QSymbolicsBase/basic_ops_homogeneous.jl#L158-L180), [QS rewrite rules](https://github.com/QuantumSavory/QuantumSymbolics.jl/blob/v0.4.20/src/QSymbolicsBase/rules.jl); SQA `docs/src/dev/architecture.md:30,51-79`, `docs/src/manual/many_body.md:22-106`.

SQA also supplies transformations of time-dependent Hamiltonians and Lindblad generators, and symbolic conversion of sandwich maps to GKSL coordinates without a Hilbert-space cutoff. QS supplies generic superoperator application and Kraus-channel expressions. No equivalent dedicated dynamical frame/GKSL API was found in QS source. Do not say QS lacks superoperators. SQA `docs/src/manual/unitary_transformations.md:31-72`, `docs/src/manual/open_systems.md:142-208`; [QS superoperator source](https://github.com/QuantumSavory/QuantumSymbolics.jl/blob/434c8c64fc5eeb95cf1e2e79ead5f26f3eed8c33/src/QSymbolicsBase/basic_superops.jl).

### 3. Build-versus-contribute rationale

The defensible historical argument is preserving and extracting working QuantumCumulants infrastructure, not claiming the alternatives are intrinsically incapable. Replacing its algebra would require translating operator identities, index semantics, expectation representations, and downstream equation-generation interfaces. Extracting it makes that investment reusable independently of cumulant truncation, and allows performance and numeric-backend development at that layer. SQA `README.md:25`, `docs/src/dev/architecture.md:18-30`.

This migration-cost/continuity argument is a reasoned inference from the differing representations, not a measured migration study. The author can establish that it motivated the original decision. It should not be presented as proof that contributing to another package was technically impossible.

The verified QC-to-SciML demonstration currently requires a pinned QC development revision; released QC/default branch does not resolve with SQA 0.14. Preserve this limitation when asserting current interoperability. `benchmark/comparison/demonstrations/README.md:70-74`, `interoperability.jl:33-50,85-86`.

## Performance claims

The refreshed kernel dataset compares SQA 0.14.0 with registered QA 1.7.0. It shows workload dependence: TC indexed commutator 22.65 microseconds SQA versus 6.49 microseconds QA; larger chain construction also favors QA. Expectation construction has different outputs and is not ratio-ranked. `benchmark/comparison/results/comparison_table.md:2,14-18` and the full JSONs.

The independent exact-coefficient checks and the same graded BCH algorithm through twentieth order provide a strong, narrow analytical example: 0.684 s SQA versus 21.834 s QA on this workstation. This transforms JC using a known generator; it does not solve arbitrary Schrieffer-Wolff/Feshbach/Floquet derivations. No measurements compare QS, so do not rank its performance from its tree-rewrite design. `benchmark/comparison/demonstrations/README.md:18-58,105-115`, `high_order.jl:35-53`.

## Candidate framing for the Julia portion of State of the field

Within Julia, QuantumAlgebra.jl provides closely related tools for indexed operator algebra, normal ordering, and expectation and cumulant equations, including support for Symbolics.jl coefficients. QuantumSymbolics.jl provides symbolic states, operators, gates, and channels with translation to several numerical representations in the QuantumSavory ecosystem. SecondQuantizedAlgebra.jl develops the algebra originating in QuantumCumulants.jl as an independent foundation for analytical and dynamical calculations. Extracting this established infrastructure retains the operator and expectation representations used by QuantumCumulants and its ModelingToolkit workflow, while making them available independently of cumulant expansions. Its eagerly reduced polynomial arithmetic and symbolic site-sum machinery address repeated operator calculations, and numerical backends connect the resulting models to direct simulations. Reproducible comparisons provide workload-specific evidence for its computational efficiency.

This is proposed framing, not final copy. Performance details can go in impact evidence, while projector/index representation trade-offs belong primarily in Software design. Cite related packages even where their goals overlap strongly.

## Claims to avoid

- “Only SQA integrates Symbolics/SciML/the quantum ecosystem.” Distinguish actual supplied workflows from general compatibility.
- “QA lacks cumulants, finite-level systems, spins, symbolic sums, or symbolic coefficients.” All present, subject to version notes above.
- “QS lacks quantum-optics or QuantumToolbox conversion.” Both present in current source.
- “QS is merely a circuit tool” or “cannot do continuous-variable algebra.” It includes bosonic objects and selected rules/representations.
- “SQA is faster” without workload/version/statistic; “high-order methods are unavailable in QA”; or any performance comparison with QS without measurements.
- “SQA is more extensible” without specifying the dimension: external numeric backends are supported, but new operator kinds require internal changes.
- “The ecosystem is already unified on SQA.” Existing tested integrations and the ambition for further adoption are different claims.
