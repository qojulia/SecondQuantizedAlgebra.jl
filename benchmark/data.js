window.BENCHMARK_DATA = {
  "lastUpdate": 1781339748409,
  "repoUrl": "https://github.com/qojulia/SecondQuantizedAlgebra.jl",
  "entries": {
    "Benchmark Results": [
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "3ce0affe7b361f362e542796fc54849758c85bec",
          "message": "SQA rewrite (#99)\n\n* feat: add abstract type hierarchy (QField, QSym, QTerm)\n\nMove old implementation to src/legacy/ and test/legacy/ for reference.\n\n* feat: add FockSpace and ProductSpace with concrete Tuple storage\n\n* feat: add Destroy and Create operators with concrete fields\n\n* feat: add QMul with lazy multiplication, canonical ordering\n\n* feat: add QAdd with type-stable addition and distribution\n\n* feat: add @qnumbers macro and TermInterface integration\n\nRename legacy test files to .jl.bak to exclude from ParallelTestRunner.\n\n* feat: add normal_order() and simplify()\n\ncanonical_lt now only sorts by space_index, preserving operator order\nwithin a space so normal_order can detect and apply commutation relations.\n\n* feat: add REPL printing and LaTeX rendering\n\n* feat: add to_numeric and numeric_average via QuantumOpticsBase\n\n* test: add code quality and integration tests\n\nClean up stale imports. Aqua, JET, ExplicitImports, CheckConcreteStructs\nall passing. Integration tests cover full workflow.\n\n* chore: remove unused imports from module header\n\nRemoved SciMLBase, MacroTools, LaTeXStrings, Combinatorics imports\nthat are not used by the new implementation. Deps kept in Project.toml\nfor future feature additions.\n\n* refactor: split code quality tests into separate files\n\nAqua, ExplicitImports, JET, and CheckConcreteStructs each in their own\ntest file for parallel execution. Concrete test auto-discovers structs.\n\n* fix: JET clean, Aqua strict, fix QAdd negation type stability\n\n- Use JET.test_package with target_modules to filter upstream reports\n- Aqua.test_all without exceptions (no ambiguities/piracies allowed)\n- Fix QAdd negation to produce Vector{QMul{T}} not Vector{QMul}\n- Linter formatting cleanup across source files\n\n* refactor: update canonical ordering and simplify operator handling in QAdd and QMul\n\n* test printing\n\n* other hilbertspaces\n\n* Implement PhaseSpace, Position, and Momentum operators in SecondQuantizedAlgebra\n\n* force type-stabilty and allocations tests\n\n* commutator\n\n* cluster spaces\n\n* indexing\n\n* add order_by_index support\n\n* benchmarks\n\n* make sort stable\n\n* some perf improvements\n\n* all structs concrete\n\n* remove skipped concrete structs from tests\n\n* format\n\n* average\n\n* numeric_average, fundamental_operators, find_operators, conj, adjoint\n\n* use Combinatorics as it was already used anyway\n\n* clean up\n\n* port tests\n\n* fixing bugs\n\n* more fixes\n\n* more tests\n\n* LaTeX rendering and add tests for symbolic variables and indexed expressions\n\n* dict based QAdd. Sick performance improvements\n\n* Implement symbolic substitution functionality and enhance NLevelSpace with symbolic levels support\n\n* Enhance indexing functionality: implement insert_index and update numeric conversions with ranges support\n\n* Implement eager diagonal splitting in Σ and QAdd multiplication; update tests for index handling and constraints\n\n* Add legacy tests for nested sum equivalences and symbolic N variants; address issues #221 and #223\n\n* remove done todo's\n\n* Limit changelog to last 3 benchmark runs\n\n* switch to runic\n\n* remove QMul\n\n* remove legacy tests\n\n* get rid of reduce\n\n* get rid of mtk weakdep\n\n* fix tests on lts\n\n* update docs\n\n* docstrings\n\n* add extend-exclude configuration to .typos.toml\n\n* fix: correct formatting in commutator documentation\n\n* format\n\n* remove claude skills\n\n* Schrieffer-Wolff example\n\n* fix docs\n\n* add precompilation workload\n\n* update PrecompileTools compatibility and fix typo in precompile comments\n\n* fix benchmark directory references in workflow configuration\n\n* fix printing\n\n* narrow type signatures in methods\n\n* prune some tests\n\n* update SW example to use `numeric_average`\n\n* fix sbstitute type piracy\n\n* Superradiant laser commutators\n\n* remove cluster (space)\n\n* remove spin\n\n* add Destroy(h::ProductSpace, name::Symbol) convenience function\n\n* anticommutator\n\n* remove expand_sums\n\n* canonical order\n\n* Port PR #103\n\n* Add overloads for numeric_average to handle SymbolicUtils.BasicSymbolic\n\n* add QField to docs\n\n* Add OrderingConvention documentation to API.md\n\n* refactor: no_equal on QTerm level\n\n* Update TODOs.md\n\n* feat: integrate ScopedValues for improved ordering management\n\n* OrderedTerm instead of Tuple{CNum, Vector{QSym}}\n\n* refactor: restructure QTerm and QAdd types, update documentation and tests\n\n* refactor: streamline numeric_average function\n\n* refactor: optimize change_index function to check for variable presence before substitution\n\n* refactor: simplify algebraic reduction functions\n\n* update docs\n\n* eager canonical form *\n\n* QTerm and QAdd handling for physical operator tracking\n\n- Updated normal_order.jl to track physical operators during normal ordering.\n- Modified qadd_arithmetic.jl to ensure physical operators are preserved in multiplication and accumulation.\n- Introduced _phys_sort! function to maintain physical order of operators.\n- Enhanced _apply_diag_terms! and _emit_diagonal! to utilize physical operators.\n- Updated simplify.jl to simplify prefactors while preserving physical operator information.\n- Adjusted substitute.jl to handle physical operators during substitution.\n- Improved cnum.jl with enhanced zero-checking for symbolic numbers.\n- Refined index.jl to ensure physical operators are correctly indexed during changes.\n- Updated qadd.jl to maintain physical operator order in QAdd operations.\n- Enhanced qterm.jl to include physical operators in QTerm structure and equality checks.\n- Modified tests in indexing_test.jl and integration_test.jl to validate new physical operator handling.\n\n* refactor: streamline physical operator tracking in arithmetic and expression functions\n\n* update devdocs\n\n* remove lazy ordering\n\n* refactor with a pipeline desing\n\n* docs and migration guide\n\n* clean up tests\n\n* update public API docstrings\n\n* proper Hermitian conjugation helpers\n\n* update more docs\n\n* fix docs\n\n* type stability pass\n\n* more type-stabilty\n\n* clarify usage of numeric_average for symbolic scalar expressions and operator expressions\n\n* more type stability\n\n* refine JET report_opt tests with structured tolerance buckets for operator families\n\n* fix docs\n\n* add examples\n\n* format\n\n* refactor _latex_prefactor to improve handling of symbolic and real components\n\n* apply codecoverage findings: dead code and coverage gaps\n\n* fix docs\n\n* improve LaTeX rendering for compound names\n\n* cleanup examples\n\n* fix conj on average\n\n* better error for index index multiplication\n\n* feat: proper avg printing in terminal\n\n* fix latex printing\n\n* refactor multiplication operations to absorb pinned summation indices\n\n* chistoph bugs and more test\n\n* chrisoph bug and more tests\n\n* fix: collapse sums\n\n* add  (i::Index)(k::Integer)\n\n* feat: enforce dead-NE invariant in QAdd constructor\n\nStrip non-equal pairs that reference no operator index, coefficient\nindex, or enclosing sum scope at QAdd construction. Pairs failing all\nthree tests encode nothing observable and would only obstruct dict\ndedup, so the inner constructor now routes through _prune_dead_ne,\nwhich walks each term, drops dead pairs, and merges colliding keys via\n_addto_key!. Idempotent on already-clean input.\n\nWith the invariant in place, average(::QAdd) drops the isempty(term.ne)\nguard when collapsing a singleton no-scope term to its bare operator:\nany surviving NE on such a term is observationally inert (no second op\nto constrain, no sum range to filter), and wrapping it in a _single_qadd\nonly obstructs downstream dict lookups.\n\nDocument the invariant in devdocs.md and link it from symbolic_sums.md.\n\n* fix jet\n\n* fix: NE constraint becomes contradictory under the rename",
          "timestamp": "2026-06-01T10:23:13-04:00",
          "tree_id": "55b2ac78c55a8f8a569ae721e1cb127f361020a9",
          "url": "https://github.com/qojulia/SecondQuantizedAlgebra.jl/commit/3ce0affe7b361f362e542796fc54849758c85bec"
        },
        "date": 1780324427795,
        "tool": "julia",
        "benches": [
          {
            "name": "Commutator/Nested JC/depth=1",
            "value": 34452,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19040\nallocs=422\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=2",
            "value": 274353,
            "unit": "ns",
            "extra": "gctime=0\nmemory=92960\nallocs=2334\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=3",
            "value": 779052,
            "unit": "ns",
            "extra": "gctime=0\nmemory=239600\nallocs=6107\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=4",
            "value": 1599216,
            "unit": "ns",
            "extra": "gctime=0\nmemory=482320\nallocs=12259\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=5",
            "value": 2803863,
            "unit": "ns",
            "extra": "gctime=0\nmemory=838912\nallocs=21243\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, V]",
            "value": 167281,
            "unit": "ns",
            "extra": "gctime=0\nmemory=58800\nallocs=1402\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, [S, H0]]",
            "value": 527413,
            "unit": "ns",
            "extra": "gctime=0\nmemory=172336\nallocs=4189\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_Dicke, S_j]",
            "value": 109810.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=73248\nallocs=1365\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_JC, σ_j]",
            "value": 96616,
            "unit": "ns",
            "extra": "gctime=0\nmemory=68160\nallocs=1111\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/double-sum spin-spin",
            "value": 15974,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8336\nallocs=126\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/indexed JC H",
            "value": 30696,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13536\nallocs=262\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/double Σ_ij(J_ij*S_i*S_j)",
            "value": 12238,
            "unit": "ns",
            "extra": "gctime=0\nmemory=15888\nallocs=179\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/single Σ_i(σ_i*σ_j)",
            "value": 10095,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12208\nallocs=139\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=2",
            "value": 2600.5555555555557,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2928\nallocs=59\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=3",
            "value": 5137.714285714285,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5248\nallocs=119\nparams={\"evals\":7,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=4",
            "value": 8572.666666666666,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8448\nallocs=201\nparams={\"evals\":3,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=5",
            "value": 13050,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12528\nallocs=305\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Ground state/3-level rewrite",
            "value": 2977.8888888888887,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3152\nallocs=47\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Multi-mode/2-mode 6-op chain",
            "value": 7321,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6608\nallocs=154\nparams={\"evals\":5,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H",
            "value": 4121.25,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4048\nallocs=82\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H²",
            "value": 801115,
            "unit": "ns",
            "extra": "gctime=0\nmemory=318624\nallocs=6773\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H",
            "value": 4506.714285714285,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4512\nallocs=95\nparams={\"evals\":7,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H²",
            "value": 1274648,
            "unit": "ns",
            "extra": "gctime=0\nmemory=512824\nallocs=11403\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H",
            "value": 91197,
            "unit": "ns",
            "extra": "gctime=0\nmemory=46496\nallocs=887\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H²",
            "value": 1459095,
            "unit": "ns",
            "extra": "gctime=0\nmemory=599104\nallocs=12867\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "4e4e6d7854e238f40766917c521721074e4ed4a9",
          "message": "feat: LaTeX rendering for averages (#157)",
          "timestamp": "2026-06-03T18:56:14+02:00",
          "tree_id": "c65dfd8033567826a727cec293c3c02b6a725627",
          "url": "https://github.com/qojulia/SecondQuantizedAlgebra.jl/commit/4e4e6d7854e238f40766917c521721074e4ed4a9"
        },
        "date": 1780506274564,
        "tool": "julia",
        "benches": [
          {
            "name": "Commutator/Nested JC/depth=1",
            "value": 32809,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19968\nallocs=441\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=2",
            "value": 285050,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106432\nallocs=2884\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=3",
            "value": 827373,
            "unit": "ns",
            "extra": "gctime=0\nmemory=280272\nallocs=7817\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=4",
            "value": 1701372,
            "unit": "ns",
            "extra": "gctime=0\nmemory=566352\nallocs=15838\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=5",
            "value": 2963100,
            "unit": "ns",
            "extra": "gctime=0\nmemory=984336\nallocs=27500\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, V]",
            "value": 175745,
            "unit": "ns",
            "extra": "gctime=0\nmemory=62704\nallocs=1773\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, [S, H0]]",
            "value": 567336,
            "unit": "ns",
            "extra": "gctime=0\nmemory=190656\nallocs=5556\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_Dicke, S_j]",
            "value": 107512.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75296\nallocs=1454\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_JC, σ_j]",
            "value": 95675,
            "unit": "ns",
            "extra": "gctime=0\nmemory=70096\nallocs=1227\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/double-sum spin-spin",
            "value": 17566,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8544\nallocs=137\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/indexed JC H",
            "value": 33871,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13888\nallocs=280\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/double Σ_ij(J_ij*S_i*S_j)",
            "value": 12228,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16072\nallocs=199\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/single Σ_i(σ_i*σ_j)",
            "value": 10055,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12312\nallocs=152\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=2",
            "value": 2435.8888888888887,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2928\nallocs=59\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=3",
            "value": 4740,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5248\nallocs=119\nparams={\"evals\":7,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=4",
            "value": 7866.75,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8448\nallocs=201\nparams={\"evals\":4,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=5",
            "value": 11788,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12528\nallocs=305\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Ground state/3-level rewrite",
            "value": 2862.1111111111113,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3168\nallocs=46\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Multi-mode/2-mode 6-op chain",
            "value": 6768.2,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6608\nallocs=154\nparams={\"evals\":5,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H",
            "value": 4054.875,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4048\nallocs=82\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H²",
            "value": 860929,
            "unit": "ns",
            "extra": "gctime=0\nmemory=338752\nallocs=7843\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H",
            "value": 4368.142857142857,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4512\nallocs=95\nparams={\"evals\":7,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H²",
            "value": 1328209,
            "unit": "ns",
            "extra": "gctime=0\nmemory=542648\nallocs=13120\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H",
            "value": 97101.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=48544\nallocs=1001\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H²",
            "value": 1551814,
            "unit": "ns",
            "extra": "gctime=0\nmemory=621360\nallocs=15156\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "49699333+dependabot[bot]@users.noreply.github.com",
            "name": "dependabot[bot]",
            "username": "dependabot[bot]"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "27e80cdebc5bbd13c5877ef524619aa103ceb909",
          "message": "build(deps): bump julia-actions/setup-julia from 2 to 3 (#158)\n\nBumps [julia-actions/setup-julia](https://github.com/julia-actions/setup-julia) from 2 to 3.\n- [Release notes](https://github.com/julia-actions/setup-julia/releases)\n- [Commits](https://github.com/julia-actions/setup-julia/compare/v2...v3)\n\n---\nupdated-dependencies:\n- dependency-name: julia-actions/setup-julia\n  dependency-version: '3'\n  dependency-type: direct:production\n  update-type: version-update:semver-major\n...\n\nSigned-off-by: dependabot[bot] <support@github.com>\nCo-authored-by: dependabot[bot] <49699333+dependabot[bot]@users.noreply.github.com>",
          "timestamp": "2026-06-07T19:14:05+02:00",
          "tree_id": "17376355d1f61b2d8019a7691c1396239a940376",
          "url": "https://github.com/qojulia/SecondQuantizedAlgebra.jl/commit/27e80cdebc5bbd13c5877ef524619aa103ceb909"
        },
        "date": 1780852918751,
        "tool": "julia",
        "benches": [
          {
            "name": "Commutator/Nested JC/depth=1",
            "value": 31677,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19216\nallocs=397\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=2",
            "value": 259874,
            "unit": "ns",
            "extra": "gctime=0\nmemory=95856\nallocs=2282\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=3",
            "value": 743439,
            "unit": "ns",
            "extra": "gctime=0\nmemory=248832\nallocs=6038\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=4",
            "value": 1536111,
            "unit": "ns",
            "extra": "gctime=0\nmemory=501568\nallocs=12179\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=5",
            "value": 2673711.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=871824\nallocs=21151\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, V]",
            "value": 159205,
            "unit": "ns",
            "extra": "gctime=0\nmemory=56272\nallocs=1403\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, [S, H0]]",
            "value": 504806,
            "unit": "ns",
            "extra": "gctime=0\nmemory=166480\nallocs=4196\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_Dicke, S_j]",
            "value": 104485,
            "unit": "ns",
            "extra": "gctime=0\nmemory=73264\nallocs=1336\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_JC, σ_j]",
            "value": 92651.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=67920\nallocs=1097\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/double-sum spin-spin",
            "value": 15998.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8336\nallocs=126\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/indexed JC H",
            "value": 30665,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13536\nallocs=262\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/double Σ_ij(J_ij*S_i*S_j)",
            "value": 12388,
            "unit": "ns",
            "extra": "gctime=0\nmemory=15784\nallocs=181\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/single Σ_i(σ_i*σ_j)",
            "value": 10145,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12104\nallocs=139\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=2",
            "value": 2437,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2928\nallocs=59\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=3",
            "value": 4641.142857142857,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5248\nallocs=119\nparams={\"evals\":7,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=4",
            "value": 7631.375,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8448\nallocs=201\nparams={\"evals\":4,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=5",
            "value": 11407,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12528\nallocs=305\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Ground state/3-level rewrite",
            "value": 2850.8888888888887,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3152\nallocs=45\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Multi-mode/2-mode 6-op chain",
            "value": 5815.333333333333,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6608\nallocs=154\nparams={\"evals\":6,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H",
            "value": 3885.75,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4048\nallocs=82\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H²",
            "value": 748360,
            "unit": "ns",
            "extra": "gctime=0\nmemory=316192\nallocs=6736\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H",
            "value": 4244.857142857143,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4512\nallocs=95\nparams={\"evals\":7,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H²",
            "value": 1164560,
            "unit": "ns",
            "extra": "gctime=0\nmemory=502504\nallocs=11188\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H",
            "value": 85025,
            "unit": "ns",
            "extra": "gctime=0\nmemory=45920\nallocs=875\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H²",
            "value": 1266531,
            "unit": "ns",
            "extra": "gctime=0\nmemory=563808\nallocs=12458\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "4186da6d9e4c8a4051d77774954fca8d72c8fbb3",
          "message": "propegate the inequility contrraint (#161)\n\n* propegate the inequility contrraint\n\n* update changelog\n\n* format changelog",
          "timestamp": "2026-06-07T19:14:29+02:00",
          "tree_id": "76b426cd14c376cf8fba5a49dbc70d26ecb11a12",
          "url": "https://github.com/qojulia/SecondQuantizedAlgebra.jl/commit/4186da6d9e4c8a4051d77774954fca8d72c8fbb3"
        },
        "date": 1780853038291,
        "tool": "julia",
        "benches": [
          {
            "name": "Commutator/Nested JC/depth=1",
            "value": 32780.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19968\nallocs=441\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=2",
            "value": 275893,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106432\nallocs=2884\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=3",
            "value": 782619,
            "unit": "ns",
            "extra": "gctime=0\nmemory=280272\nallocs=7817\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=4",
            "value": 1607708.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=566352\nallocs=15838\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=5",
            "value": 2821510,
            "unit": "ns",
            "extra": "gctime=0\nmemory=984336\nallocs=27500\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, V]",
            "value": 168548.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=62704\nallocs=1773\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, [S, H0]]",
            "value": 541439,
            "unit": "ns",
            "extra": "gctime=0\nmemory=190656\nallocs=5556\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_Dicke, S_j]",
            "value": 102643.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75296\nallocs=1454\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_JC, σ_j]",
            "value": 91565,
            "unit": "ns",
            "extra": "gctime=0\nmemory=70096\nallocs=1227\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/double-sum spin-spin",
            "value": 15366,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8528\nallocs=136\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/indexed JC H",
            "value": 29845,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13856\nallocs=278\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/double Σ_ij(J_ij*S_i*S_j)",
            "value": 11332,
            "unit": "ns",
            "extra": "gctime=0\nmemory=15496\nallocs=187\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/single Σ_i(σ_i*σ_j)",
            "value": 9456.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12120\nallocs=148\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=2",
            "value": 2391.222222222222,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2928\nallocs=59\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=3",
            "value": 4668,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5248\nallocs=119\nparams={\"evals\":7,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=4",
            "value": 7607,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8448\nallocs=201\nparams={\"evals\":4,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=5",
            "value": 11266,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12528\nallocs=305\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Ground state/3-level rewrite",
            "value": 2934.3888888888887,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3168\nallocs=46\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Multi-mode/2-mode 6-op chain",
            "value": 6426.9,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6608\nallocs=154\nparams={\"evals\":5,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H",
            "value": 4073,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4048\nallocs=82\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H²",
            "value": 770394,
            "unit": "ns",
            "extra": "gctime=0\nmemory=331808\nallocs=7692\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H",
            "value": 4343,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4512\nallocs=95\nparams={\"evals\":7,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H²",
            "value": 1193250,
            "unit": "ns",
            "extra": "gctime=0\nmemory=529160\nallocs=12830\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H",
            "value": 86457,
            "unit": "ns",
            "extra": "gctime=0\nmemory=47664\nallocs=982\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H²",
            "value": 1291368,
            "unit": "ns",
            "extra": "gctime=0\nmemory=600048\nallocs=14699\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "1d963dc5875299c8aa5f1ab9a2d3b020c63a9260",
          "message": "perf: for commuting commutators and symbolic rational coefficients (#162)\n\n* perf: for commuting commutators and symbolic rational coefficients\n\n* format changelog",
          "timestamp": "2026-06-09T17:38:15+02:00",
          "tree_id": "8fa34a72f3e8df97dff4c4b7041b943a21d03200",
          "url": "https://github.com/qojulia/SecondQuantizedAlgebra.jl/commit/1d963dc5875299c8aa5f1ab9a2d3b020c63a9260"
        },
        "date": 1781020111087,
        "tool": "julia",
        "benches": [
          {
            "name": "Commutator/Nested JC/depth=1",
            "value": 26530,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26032\nallocs=379\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=2",
            "value": 290086,
            "unit": "ns",
            "extra": "gctime=0\nmemory=162096\nallocs=3122\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=3",
            "value": 872382.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=428944\nallocs=8841\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=4",
            "value": 1826118,
            "unit": "ns",
            "extra": "gctime=0\nmemory=847568\nallocs=18153\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=5",
            "value": 3231092,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1452992\nallocs=31988\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, V]",
            "value": 204377,
            "unit": "ns",
            "extra": "gctime=0\nmemory=87832\nallocs=2139\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, [S, H0]]",
            "value": 446271,
            "unit": "ns",
            "extra": "gctime=0\nmemory=180832\nallocs=4548\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_Dicke, S_j]",
            "value": 103285,
            "unit": "ns",
            "extra": "gctime=0\nmemory=81648\nallocs=1371\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_JC, σ_j]",
            "value": 89615,
            "unit": "ns",
            "extra": "gctime=0\nmemory=76384\nallocs=1144\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/double-sum spin-spin",
            "value": 16715,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8480\nallocs=135\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/indexed JC H",
            "value": 31647,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13760\nallocs=276\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/double Σ_ij(J_ij*S_i*S_j)",
            "value": 12689,
            "unit": "ns",
            "extra": "gctime=0\nmemory=15496\nallocs=187\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/single Σ_i(σ_i*σ_j)",
            "value": 10078.333333333334,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12120\nallocs=148\nparams={\"evals\":3,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=2",
            "value": 2713,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2928\nallocs=59\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=3",
            "value": 5153.571428571428,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5248\nallocs=119\nparams={\"evals\":7,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=4",
            "value": 8523,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8448\nallocs=201\nparams={\"evals\":3,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=5",
            "value": 12800,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12528\nallocs=305\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Ground state/3-level rewrite",
            "value": 2888.777777777778,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3168\nallocs=46\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Multi-mode/2-mode 6-op chain",
            "value": 6890.4,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6608\nallocs=154\nparams={\"evals\":5,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H",
            "value": 4177.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4048\nallocs=82\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H²",
            "value": 776183,
            "unit": "ns",
            "extra": "gctime=0\nmemory=331808\nallocs=7692\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H",
            "value": 4724.285714285715,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4512\nallocs=95\nparams={\"evals\":7,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H²",
            "value": 1202584,
            "unit": "ns",
            "extra": "gctime=0\nmemory=529160\nallocs=12830\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H",
            "value": 86850,
            "unit": "ns",
            "extra": "gctime=0\nmemory=47664\nallocs=982\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H²",
            "value": 1314948.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=600048\nallocs=14699\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "49699333+dependabot[bot]@users.noreply.github.com",
            "name": "dependabot[bot]",
            "username": "dependabot[bot]"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "08fb182eb59d45e74d9ae1d08724fd9ae2b14660",
          "message": "build(deps): bump codecov/codecov-action from 6 to 7 (#166)\n\nBumps [codecov/codecov-action](https://github.com/codecov/codecov-action) from 6 to 7.\n- [Release notes](https://github.com/codecov/codecov-action/releases)\n- [Changelog](https://github.com/codecov/codecov-action/blob/main/CHANGELOG.md)\n- [Commits](https://github.com/codecov/codecov-action/compare/v6...v7)\n\n---\nupdated-dependencies:\n- dependency-name: codecov/codecov-action\n  dependency-version: '7'\n  dependency-type: direct:production\n  update-type: version-update:semver-major\n...\n\nSigned-off-by: dependabot[bot] <support@github.com>\nCo-authored-by: dependabot[bot] <49699333+dependabot[bot]@users.noreply.github.com>",
          "timestamp": "2026-06-09T23:10:21+02:00",
          "tree_id": "0e2eea904b9c628f7736c03ac22b3110fa824889",
          "url": "https://github.com/qojulia/SecondQuantizedAlgebra.jl/commit/08fb182eb59d45e74d9ae1d08724fd9ae2b14660"
        },
        "date": 1781039628002,
        "tool": "julia",
        "benches": [
          {
            "name": "Commutator/Nested JC/depth=1",
            "value": 26419,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26032\nallocs=379\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=2",
            "value": 311630,
            "unit": "ns",
            "extra": "gctime=0\nmemory=162096\nallocs=3122\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=3",
            "value": 936032,
            "unit": "ns",
            "extra": "gctime=0\nmemory=428944\nallocs=8843\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=4",
            "value": 1944079.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=847568\nallocs=18158\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=5",
            "value": 3425849,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1452992\nallocs=31987\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, V]",
            "value": 226502,
            "unit": "ns",
            "extra": "gctime=0\nmemory=87888\nallocs=2141\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, [S, H0]]",
            "value": 502355,
            "unit": "ns",
            "extra": "gctime=0\nmemory=180864\nallocs=4548\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_Dicke, S_j]",
            "value": 108332,
            "unit": "ns",
            "extra": "gctime=0\nmemory=81648\nallocs=1371\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_JC, σ_j]",
            "value": 100968,
            "unit": "ns",
            "extra": "gctime=0\nmemory=76384\nallocs=1144\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/double-sum spin-spin",
            "value": 16641,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8480\nallocs=135\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/indexed JC H",
            "value": 31779,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13760\nallocs=276\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/double Σ_ij(J_ij*S_i*S_j)",
            "value": 12142,
            "unit": "ns",
            "extra": "gctime=0\nmemory=15496\nallocs=187\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/single Σ_i(σ_i*σ_j)",
            "value": 10389,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12120\nallocs=148\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=2",
            "value": 2613.777777777778,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2928\nallocs=59\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=3",
            "value": 5052.142857142857,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5248\nallocs=119\nparams={\"evals\":7,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=4",
            "value": 8339,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8448\nallocs=201\nparams={\"evals\":3,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=5",
            "value": 12383,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12528\nallocs=305\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Ground state/3-level rewrite",
            "value": 2748.4444444444443,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3168\nallocs=46\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Multi-mode/2-mode 6-op chain",
            "value": 6844.8,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6608\nallocs=154\nparams={\"evals\":5,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H",
            "value": 3957.375,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4048\nallocs=82\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H²",
            "value": 904073,
            "unit": "ns",
            "extra": "gctime=0\nmemory=331808\nallocs=7692\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H",
            "value": 4465.428571428572,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4512\nallocs=95\nparams={\"evals\":7,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H²",
            "value": 1366465,
            "unit": "ns",
            "extra": "gctime=0\nmemory=529160\nallocs=12830\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H",
            "value": 102080,
            "unit": "ns",
            "extra": "gctime=0\nmemory=47664\nallocs=982\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H²",
            "value": 1457223,
            "unit": "ns",
            "extra": "gctime=0\nmemory=600048\nallocs=14699\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "c34c002ebfd218cf469990637dd1dbeade015a33",
          "message": "Average as mtk unknown (#168)\n\n* perf: for commuting commutators and symbolic rational coefficients\n\n* format changelog\n\n* refactor: change average symtype to number and add make_time_dependent\n\n* add changelog",
          "timestamp": "2026-06-12T09:26:59+02:00",
          "tree_id": "f956098e8ebef21c322c2716732af0e41bd278bd",
          "url": "https://github.com/qojulia/SecondQuantizedAlgebra.jl/commit/c34c002ebfd218cf469990637dd1dbeade015a33"
        },
        "date": 1781249680817,
        "tool": "julia",
        "benches": [
          {
            "name": "Commutator/Nested JC/depth=1",
            "value": 26660,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26032\nallocs=379\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=2",
            "value": 313298.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=162096\nallocs=3122\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=3",
            "value": 943160,
            "unit": "ns",
            "extra": "gctime=0\nmemory=428944\nallocs=8845\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=4",
            "value": 1961677.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=847568\nallocs=18158\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=5",
            "value": 3470196,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1452992\nallocs=31987\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, V]",
            "value": 229638,
            "unit": "ns",
            "extra": "gctime=0\nmemory=87832\nallocs=2139\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, [S, H0]]",
            "value": 509361,
            "unit": "ns",
            "extra": "gctime=0\nmemory=180832\nallocs=4548\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_Dicke, S_j]",
            "value": 108522,
            "unit": "ns",
            "extra": "gctime=0\nmemory=81648\nallocs=1371\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_JC, σ_j]",
            "value": 100863.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=76384\nallocs=1144\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/double-sum spin-spin",
            "value": 22372,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9184\nallocs=153\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/indexed JC H",
            "value": 42550,
            "unit": "ns",
            "extra": "gctime=0\nmemory=15168\nallocs=312\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/double Σ_ij(J_ij*S_i*S_j)",
            "value": 12033,
            "unit": "ns",
            "extra": "gctime=0\nmemory=15496\nallocs=187\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/single Σ_i(σ_i*σ_j)",
            "value": 10399,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12120\nallocs=148\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=2",
            "value": 2518,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2928\nallocs=59\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=3",
            "value": 4841.857142857143,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5248\nallocs=119\nparams={\"evals\":7,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=4",
            "value": 7982.25,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8448\nallocs=201\nparams={\"evals\":4,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=5",
            "value": 11892,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12528\nallocs=305\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Ground state/3-level rewrite",
            "value": 2775.222222222222,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3168\nallocs=46\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Multi-mode/2-mode 6-op chain",
            "value": 6657.4,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6608\nallocs=154\nparams={\"evals\":5,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H",
            "value": 3995,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4048\nallocs=82\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H²",
            "value": 1002521,
            "unit": "ns",
            "extra": "gctime=0\nmemory=345248\nallocs=8028\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H",
            "value": 4408.285714285715,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4512\nallocs=95\nparams={\"evals\":7,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H²",
            "value": 1502383,
            "unit": "ns",
            "extra": "gctime=0\nmemory=547080\nallocs=13278\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H",
            "value": 110647,
            "unit": "ns",
            "extra": "gctime=0\nmemory=48944\nallocs=1014\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H²",
            "value": 1576132,
            "unit": "ns",
            "extra": "gctime=0\nmemory=613488\nallocs=15035\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "c34c002ebfd218cf469990637dd1dbeade015a33",
          "message": "Average as mtk unknown (#168)\n\n* perf: for commuting commutators and symbolic rational coefficients\n\n* format changelog\n\n* refactor: change average symtype to number and add make_time_dependent\n\n* add changelog",
          "timestamp": "2026-06-12T09:26:59+02:00",
          "tree_id": "f956098e8ebef21c322c2716732af0e41bd278bd",
          "url": "https://github.com/qojulia/SecondQuantizedAlgebra.jl/commit/c34c002ebfd218cf469990637dd1dbeade015a33"
        },
        "date": 1781250098372,
        "tool": "julia",
        "benches": [
          {
            "name": "Commutator/Nested JC/depth=1",
            "value": 26869.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26032\nallocs=379\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=2",
            "value": 294113,
            "unit": "ns",
            "extra": "gctime=0\nmemory=162096\nallocs=3122\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=3",
            "value": 878378,
            "unit": "ns",
            "extra": "gctime=0\nmemory=428944\nallocs=8843\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=4",
            "value": 1840029,
            "unit": "ns",
            "extra": "gctime=0\nmemory=847568\nallocs=18153\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=5",
            "value": 3253632.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1452992\nallocs=31986\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, V]",
            "value": 204782,
            "unit": "ns",
            "extra": "gctime=0\nmemory=87832\nallocs=2139\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, [S, H0]]",
            "value": 447469,
            "unit": "ns",
            "extra": "gctime=0\nmemory=180864\nallocs=4548\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_Dicke, S_j]",
            "value": 101690,
            "unit": "ns",
            "extra": "gctime=0\nmemory=81648\nallocs=1371\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_JC, σ_j]",
            "value": 90283,
            "unit": "ns",
            "extra": "gctime=0\nmemory=76384\nallocs=1144\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/double-sum spin-spin",
            "value": 19989,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9184\nallocs=153\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/indexed JC H",
            "value": 38477,
            "unit": "ns",
            "extra": "gctime=0\nmemory=15168\nallocs=312\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/double Σ_ij(J_ij*S_i*S_j)",
            "value": 12038,
            "unit": "ns",
            "extra": "gctime=0\nmemory=15496\nallocs=187\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/single Σ_i(σ_i*σ_j)",
            "value": 10255,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12120\nallocs=148\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=2",
            "value": 2670.5555555555557,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2928\nallocs=59\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=3",
            "value": 5123.214285714286,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5248\nallocs=119\nparams={\"evals\":7,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=4",
            "value": 8442.333333333334,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8448\nallocs=201\nparams={\"evals\":3,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=5",
            "value": 12708,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12528\nallocs=305\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Ground state/3-level rewrite",
            "value": 2963.222222222222,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3168\nallocs=46\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Multi-mode/2-mode 6-op chain",
            "value": 6960.2,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6608\nallocs=154\nparams={\"evals\":5,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H",
            "value": 4133.285714285715,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4048\nallocs=82\nparams={\"evals\":7,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H²",
            "value": 825827,
            "unit": "ns",
            "extra": "gctime=0\nmemory=345248\nallocs=8028\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H",
            "value": 4561,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4512\nallocs=95\nparams={\"evals\":7,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H²",
            "value": 1269279,
            "unit": "ns",
            "extra": "gctime=0\nmemory=547080\nallocs=13278\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H",
            "value": 91244,
            "unit": "ns",
            "extra": "gctime=0\nmemory=48944\nallocs=1014\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H²",
            "value": 1370318,
            "unit": "ns",
            "extra": "gctime=0\nmemory=613488\nallocs=15035\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "e4d0371b36fde1c78f11ffe9939dbcdb17685640",
          "message": "feat: add order_key total ordering for operators(#169)\n\n* feat: add order_key total ordering for operators\n\nAdd public order_key / term_order_key / qadd_order_key: a total,\nidentity-faithful structural ordering for operators, products and sums,\nbuilt from one per-type order_key method (no generic QSym fallback, so a\nnew operator type without a key is a loud MethodError). The key ties two\noperators exactly when they are isequal, giving downstream packages a\nreproducible way to pick canonical representatives and compare expressions\nwithout round-tripping through show or hash. Includes API/devdocs/changelog\nand a consistency + total-order test.\n\n* update minot\n\n* add term_order_key / qadd_order_key tests",
          "timestamp": "2026-06-13T10:32:25+02:00",
          "tree_id": "f955e67967d42394e021a0b556b7e0046a354975",
          "url": "https://github.com/qojulia/SecondQuantizedAlgebra.jl/commit/e4d0371b36fde1c78f11ffe9939dbcdb17685640"
        },
        "date": 1781339747834,
        "tool": "julia",
        "benches": [
          {
            "name": "Commutator/Nested JC/depth=1",
            "value": 26389,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26032\nallocs=379\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=2",
            "value": 288735,
            "unit": "ns",
            "extra": "gctime=0\nmemory=162096\nallocs=3122\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=3",
            "value": 869321,
            "unit": "ns",
            "extra": "gctime=0\nmemory=428944\nallocs=8845\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=4",
            "value": 1813844,
            "unit": "ns",
            "extra": "gctime=0\nmemory=847568\nallocs=18154\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=5",
            "value": 3217072,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1452992\nallocs=31983\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, V]",
            "value": 203740,
            "unit": "ns",
            "extra": "gctime=0\nmemory=87856\nallocs=2141\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, [S, H0]]",
            "value": 445487,
            "unit": "ns",
            "extra": "gctime=0\nmemory=180832\nallocs=4548\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_Dicke, S_j]",
            "value": 100808,
            "unit": "ns",
            "extra": "gctime=0\nmemory=81648\nallocs=1371\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_JC, σ_j]",
            "value": 88691,
            "unit": "ns",
            "extra": "gctime=0\nmemory=76384\nallocs=1144\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/double-sum spin-spin",
            "value": 20671,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9184\nallocs=153\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/indexed JC H",
            "value": 37715,
            "unit": "ns",
            "extra": "gctime=0\nmemory=15168\nallocs=312\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/double Σ_ij(J_ij*S_i*S_j)",
            "value": 12037,
            "unit": "ns",
            "extra": "gctime=0\nmemory=15496\nallocs=187\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/single Σ_i(σ_i*σ_j)",
            "value": 10195,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12120\nallocs=148\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=2",
            "value": 2518.1111111111113,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2928\nallocs=59\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=3",
            "value": 4822.857142857143,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5248\nallocs=119\nparams={\"evals\":7,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=4",
            "value": 7974.25,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8448\nallocs=201\nparams={\"evals\":4,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=5",
            "value": 11977,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12528\nallocs=305\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Ground state/3-level rewrite",
            "value": 2932.1111111111113,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3168\nallocs=46\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Multi-mode/2-mode 6-op chain",
            "value": 6665.8,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6608\nallocs=154\nparams={\"evals\":5,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H",
            "value": 4027.125,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4048\nallocs=82\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H²",
            "value": 818367,
            "unit": "ns",
            "extra": "gctime=0\nmemory=345248\nallocs=8028\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H",
            "value": 4412.214285714286,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4512\nallocs=95\nparams={\"evals\":7,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H²",
            "value": 1261689,
            "unit": "ns",
            "extra": "gctime=0\nmemory=547080\nallocs=13278\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H",
            "value": 90553.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=48944\nallocs=1014\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H²",
            "value": 1365694,
            "unit": "ns",
            "extra": "gctime=0\nmemory=613488\nallocs=15035\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          }
        ]
      }
    ]
  }
}