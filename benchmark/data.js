window.BENCHMARK_DATA = {
  "lastUpdate": 1784738688279,
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
          "id": "fe7d49ca0b9366b2e3b8698df1672ffa6b3a4cf8",
          "message": "feat: add index_slot (#170)",
          "timestamp": "2026-06-14T11:33:04+02:00",
          "tree_id": "44166998ef8efa5da5c3e23332d83a491fec927e",
          "url": "https://github.com/qojulia/SecondQuantizedAlgebra.jl/commit/fe7d49ca0b9366b2e3b8698df1672ffa6b3a4cf8"
        },
        "date": 1781429782835,
        "tool": "julia",
        "benches": [
          {
            "name": "Commutator/Nested JC/depth=1",
            "value": 26119,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26032\nallocs=379\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=2",
            "value": 305764,
            "unit": "ns",
            "extra": "gctime=0\nmemory=162096\nallocs=3122\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=3",
            "value": 916130.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=428944\nallocs=8845\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=4",
            "value": 1916037,
            "unit": "ns",
            "extra": "gctime=0\nmemory=847568\nallocs=18153\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=5",
            "value": 3361724,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1452992\nallocs=31982\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, V]",
            "value": 221245,
            "unit": "ns",
            "extra": "gctime=0\nmemory=87832\nallocs=2139\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, [S, H0]]",
            "value": 494869,
            "unit": "ns",
            "extra": "gctime=0\nmemory=180832\nallocs=4548\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_Dicke, S_j]",
            "value": 106059,
            "unit": "ns",
            "extra": "gctime=0\nmemory=81648\nallocs=1371\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_JC, σ_j]",
            "value": 99156,
            "unit": "ns",
            "extra": "gctime=0\nmemory=76384\nallocs=1144\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/double-sum spin-spin",
            "value": 21670,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9184\nallocs=153\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/indexed JC H",
            "value": 42099,
            "unit": "ns",
            "extra": "gctime=0\nmemory=15168\nallocs=312\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/double Σ_ij(J_ij*S_i*S_j)",
            "value": 11912,
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
            "value": 2693.8888888888887,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2928\nallocs=59\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=3",
            "value": 5278.333333333333,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5248\nallocs=119\nparams={\"evals\":6,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=4",
            "value": 8806.333333333334,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8448\nallocs=201\nparams={\"evals\":3,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=5",
            "value": 13145,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12528\nallocs=305\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Ground state/3-level rewrite",
            "value": 2750.777777777778,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3168\nallocs=46\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Multi-mode/2-mode 6-op chain",
            "value": 7217.4,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6608\nallocs=154\nparams={\"evals\":5,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H",
            "value": 4097.714285714285,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4048\nallocs=82\nparams={\"evals\":7,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H²",
            "value": 967119.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=345248\nallocs=8028\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H",
            "value": 4504.142857142857,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4512\nallocs=95\nparams={\"evals\":7,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H²",
            "value": 1455093,
            "unit": "ns",
            "extra": "gctime=0\nmemory=547080\nallocs=13278\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H",
            "value": 106540,
            "unit": "ns",
            "extra": "gctime=0\nmemory=48944\nallocs=1014\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H²",
            "value": 1513053.5,
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
          "id": "37800934f5dd0bb1253b7f57039c6cf78e4e3823",
          "message": "fix: Render lifted time-dependent averages as ⟨op⟩(t) (#173)\n\n* Render lifted time-dependent averages as ⟨op⟩(t)\n\nAdd AverageOperator metadata display hooks (REPL show_metadata + latexify\n_toexpr_metadata) so a lifted average prints as ⟨op⟩(t) instead of leaking the\nstructurally-injective _avg_… variable name. The SumIndices REPL hook defers to\nthe new hook for lifted vars.\n\n* Bump to v0.6.3, drop inline comments, broaden display tests\n\n* Update root Changelog.md (enforced by CI) for v0.6.3 and backfill v0.6.2\n\n* Generate docs/src/changelog.md on build; untrack it\n\nThe docs changelog is a copy of the root Changelog.md. Regenerate it on every\ndoc build (not just CI) and gitignore it so the tracked copy can no longer drift.\nservedocs skips the regenerated file to avoid a LiveServer rebuild loop.",
          "timestamp": "2026-06-16T11:48:27+02:00",
          "tree_id": "4c5551ab08a462dc7dbdf5dfcb6e80dd499465f7",
          "url": "https://github.com/qojulia/SecondQuantizedAlgebra.jl/commit/37800934f5dd0bb1253b7f57039c6cf78e4e3823"
        },
        "date": 1781603754001,
        "tool": "julia",
        "benches": [
          {
            "name": "Commutator/Nested JC/depth=1",
            "value": 26069,
            "unit": "ns",
            "extra": "gctime=0\nmemory=25552\nallocs=349\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=2",
            "value": 301409.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=152912\nallocs=2548\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=3",
            "value": 909503,
            "unit": "ns",
            "extra": "gctime=0\nmemory=400752\nallocs=7083\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=4",
            "value": 1888997,
            "unit": "ns",
            "extra": "gctime=0\nmemory=788368\nallocs=14455\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=5",
            "value": 3336218,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1348768\nallocs=25471\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, V]",
            "value": 212688,
            "unit": "ns",
            "extra": "gctime=0\nmemory=80736\nallocs=1694\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, [S, H0]]",
            "value": 469520,
            "unit": "ns",
            "extra": "gctime=0\nmemory=163104\nallocs=3438\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_Dicke, S_j]",
            "value": 109295,
            "unit": "ns",
            "extra": "gctime=0\nmemory=79984\nallocs=1267\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_JC, σ_j]",
            "value": 97633,
            "unit": "ns",
            "extra": "gctime=0\nmemory=74528\nallocs=1028\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/double-sum spin-spin",
            "value": 19156,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9328\nallocs=145\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/indexed JC H",
            "value": 37140,
            "unit": "ns",
            "extra": "gctime=0\nmemory=15520\nallocs=300\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/double Σ_ij(J_ij*S_i*S_j)",
            "value": 12102,
            "unit": "ns",
            "extra": "gctime=0\nmemory=15240\nallocs=171\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/single Σ_i(σ_i*σ_j)",
            "value": 9985.333333333334,
            "unit": "ns",
            "extra": "gctime=0\nmemory=11912\nallocs=135\nparams={\"evals\":3,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=2",
            "value": 2656,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2928\nallocs=59\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=3",
            "value": 5206.857142857143,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5248\nallocs=119\nparams={\"evals\":7,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=4",
            "value": 8639.333333333334,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8448\nallocs=201\nparams={\"evals\":3,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=5",
            "value": 12964,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12528\nallocs=305\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Ground state/3-level rewrite",
            "value": 2768.5555555555557,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3152\nallocs=45\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Multi-mode/2-mode 6-op chain",
            "value": 7223.6,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6608\nallocs=154\nparams={\"evals\":5,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H",
            "value": 4100.428571428572,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4048\nallocs=82\nparams={\"evals\":7,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H²",
            "value": 954247,
            "unit": "ns",
            "extra": "gctime=0\nmemory=333120\nallocs=7099\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H",
            "value": 4554.285714285715,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4512\nallocs=95\nparams={\"evals\":7,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H²",
            "value": 1424776,
            "unit": "ns",
            "extra": "gctime=0\nmemory=525384\nallocs=11675\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H",
            "value": 105232.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=47584\nallocs=910\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H²",
            "value": 1466835,
            "unit": "ns",
            "extra": "gctime=0\nmemory=581152\nallocs=12824\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
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
          "id": "061991ded5653351900128a68542786a49392d2d",
          "message": "fix: `change_index` now correctly zeros `DoubleIndexedVariable` (#174)\n\n* fix: `change_index` now correctly zeros `DoubleIndexedVariable`\n\n* make changelog",
          "timestamp": "2026-06-16T14:43:11+02:00",
          "tree_id": "5e7edc8b8ac4fc428729d925482f55023ac666fa",
          "url": "https://github.com/qojulia/SecondQuantizedAlgebra.jl/commit/061991ded5653351900128a68542786a49392d2d"
        },
        "date": 1781614354413,
        "tool": "julia",
        "benches": [
          {
            "name": "Commutator/Nested JC/depth=1",
            "value": 24993,
            "unit": "ns",
            "extra": "gctime=0\nmemory=25904\nallocs=360\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=2",
            "value": 273663.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=154672\nallocs=2603\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=3",
            "value": 815788,
            "unit": "ns",
            "extra": "gctime=0\nmemory=405104\nallocs=7219\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=4",
            "value": 1703580.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=796560\nallocs=14711\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=5",
            "value": 3008257.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1362176\nallocs=25890\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, V]",
            "value": 194512.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=81280\nallocs=1711\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, [S, H0]]",
            "value": 419621,
            "unit": "ns",
            "extra": "gctime=0\nmemory=164160\nallocs=3471\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_Dicke, S_j]",
            "value": 95957,
            "unit": "ns",
            "extra": "gctime=0\nmemory=80752\nallocs=1291\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_JC, σ_j]",
            "value": 82188,
            "unit": "ns",
            "extra": "gctime=0\nmemory=77184\nallocs=1063\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/double-sum spin-spin",
            "value": 18481,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9376\nallocs=146\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/indexed JC H",
            "value": 35666,
            "unit": "ns",
            "extra": "gctime=0\nmemory=15616\nallocs=302\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/double Σ_ij(J_ij*S_i*S_j)",
            "value": 11057,
            "unit": "ns",
            "extra": "gctime=0\nmemory=15768\nallocs=175\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/single Σ_i(σ_i*σ_j)",
            "value": 9118.666666666666,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12488\nallocs=141\nparams={\"evals\":3,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=2",
            "value": 2375.3333333333335,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2928\nallocs=59\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=3",
            "value": 4640.857142857143,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5248\nallocs=119\nparams={\"evals\":7,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=4",
            "value": 7602.125,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8448\nallocs=201\nparams={\"evals\":4,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=5",
            "value": 11255,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12528\nallocs=305\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Ground state/3-level rewrite",
            "value": 2982,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3184\nallocs=46\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Multi-mode/2-mode 6-op chain",
            "value": 6290.2,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6608\nallocs=154\nparams={\"evals\":5,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H",
            "value": 3910.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4048\nallocs=82\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H²",
            "value": 863022,
            "unit": "ns",
            "extra": "gctime=0\nmemory=333184\nallocs=7101\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H",
            "value": 4233,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4512\nallocs=95\nparams={\"evals\":7,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H²",
            "value": 1318381,
            "unit": "ns",
            "extra": "gctime=0\nmemory=525448\nallocs=11677\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H",
            "value": 95292,
            "unit": "ns",
            "extra": "gctime=0\nmemory=47584\nallocs=910\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H²",
            "value": 1394232,
            "unit": "ns",
            "extra": "gctime=0\nmemory=581152\nallocs=12824\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
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
          "id": "8c6676622425f9a99c6d96470599213f5e61652e",
          "message": "feat: support @variables g::Number (#177)\n\n* feat: support @variables g::Number\n\n* add changelog\n\n* small optimizations\n\n* fix ci\n\n* fix docs ci\n\n* fail docs on missing exported docstring api\n\n* another attempt",
          "timestamp": "2026-06-20T21:19:36+02:00",
          "tree_id": "9b59d4174eae97ec8a75c8a4938cd23ab7c6de37",
          "url": "https://github.com/qojulia/SecondQuantizedAlgebra.jl/commit/8c6676622425f9a99c6d96470599213f5e61652e"
        },
        "date": 1781983801854,
        "tool": "julia",
        "benches": [
          {
            "name": "Commutator/Nested JC/depth=1",
            "value": 26329,
            "unit": "ns",
            "extra": "gctime=0\nmemory=24720\nallocs=334\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=2",
            "value": 294431,
            "unit": "ns",
            "extra": "gctime=0\nmemory=144800\nallocs=2362\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=3",
            "value": 876752,
            "unit": "ns",
            "extra": "gctime=0\nmemory=376928\nallocs=6510\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=4",
            "value": 1831590,
            "unit": "ns",
            "extra": "gctime=0\nmemory=739216\nallocs=13250\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=5",
            "value": 3254022,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1263328\nallocs=23356\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, V]",
            "value": 216906,
            "unit": "ns",
            "extra": "gctime=0\nmemory=81040\nallocs=1590\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, [S, H0]]",
            "value": 471536,
            "unit": "ns",
            "extra": "gctime=0\nmemory=161200\nallocs=3169\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_Dicke, S_j]",
            "value": 109595,
            "unit": "ns",
            "extra": "gctime=0\nmemory=77168\nallocs=1215\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_JC, σ_j]",
            "value": 102602,
            "unit": "ns",
            "extra": "gctime=0\nmemory=74912\nallocs=1015\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/double-sum spin-spin",
            "value": 21480,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9328\nallocs=145\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/indexed JC H",
            "value": 41117,
            "unit": "ns",
            "extra": "gctime=0\nmemory=15520\nallocs=300\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/double Σ_ij(J_ij*S_i*S_j)",
            "value": 12774,
            "unit": "ns",
            "extra": "gctime=0\nmemory=15736\nallocs=173\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/single Σ_i(σ_i*σ_j)",
            "value": 11221,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12392\nallocs=139\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=2",
            "value": 2557,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2928\nallocs=59\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=3",
            "value": 4926.285714285715,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5248\nallocs=119\nparams={\"evals\":7,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=4",
            "value": 8135.25,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8448\nallocs=201\nparams={\"evals\":4,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=5",
            "value": 12254,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12528\nallocs=305\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Ground state/3-level rewrite",
            "value": 2870.8888888888887,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3088\nallocs=44\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Multi-mode/2-mode 6-op chain",
            "value": 6742.6,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6608\nallocs=154\nparams={\"evals\":5,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H",
            "value": 4042.625,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4048\nallocs=82\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H²",
            "value": 935531,
            "unit": "ns",
            "extra": "gctime=0\nmemory=328128\nallocs=6969\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H",
            "value": 4432.428571428572,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4512\nallocs=95\nparams={\"evals\":7,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H²",
            "value": 1401808,
            "unit": "ns",
            "extra": "gctime=0\nmemory=521512\nallocs=11578\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H",
            "value": 104636,
            "unit": "ns",
            "extra": "gctime=0\nmemory=47216\nallocs=901\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H²",
            "value": 1451228,
            "unit": "ns",
            "extra": "gctime=0\nmemory=603104\nallocs=12758\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
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
          "id": "9184fd79e75165390d59a21a3c7825591f6efb12",
          "message": "Averaging an indexed sum (#179)\n\n* fix: Averaging an indexed sum\n\n* remove  verbose comments\n\n* fix env\n\n* fix: order-insensitive QAdd equality on bound sum indices\n\nQAdd isequal/hash compared .indices order-sensitively, but bound\nsummation indices are a set (Σ_iΣ_j ≡ Σ_jΣ_i) with no canonical order,\nso undo_average(average(Σ_i a_i + Σ_j a_j')) == s4 was seed-dependent.\nCompare/hash .indices as a set. Also document is_indexed_sum in API.md.",
          "timestamp": "2026-06-21T11:09:28+02:00",
          "tree_id": "912780a9c07fa62cf9af844e129fc63901109786",
          "url": "https://github.com/qojulia/SecondQuantizedAlgebra.jl/commit/9184fd79e75165390d59a21a3c7825591f6efb12"
        },
        "date": 1782033378152,
        "tool": "julia",
        "benches": [
          {
            "name": "Commutator/Nested JC/depth=1",
            "value": 25317,
            "unit": "ns",
            "extra": "gctime=0\nmemory=24720\nallocs=334\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=2",
            "value": 285083,
            "unit": "ns",
            "extra": "gctime=0\nmemory=144800\nallocs=2362\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=3",
            "value": 867511,
            "unit": "ns",
            "extra": "gctime=0\nmemory=376928\nallocs=6510\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=4",
            "value": 1811525.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=739216\nallocs=13250\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=5",
            "value": 3214707,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1263328\nallocs=23356\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, V]",
            "value": 219049,
            "unit": "ns",
            "extra": "gctime=0\nmemory=81040\nallocs=1590\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, [S, H0]]",
            "value": 481871.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=161200\nallocs=3169\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_Dicke, S_j]",
            "value": 103133,
            "unit": "ns",
            "extra": "gctime=0\nmemory=77168\nallocs=1215\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_JC, σ_j]",
            "value": 94216,
            "unit": "ns",
            "extra": "gctime=0\nmemory=74912\nallocs=1015\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/double-sum spin-spin",
            "value": 19868,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9328\nallocs=145\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/indexed JC H",
            "value": 38412,
            "unit": "ns",
            "extra": "gctime=0\nmemory=15520\nallocs=300\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/double Σ_ij(J_ij*S_i*S_j)",
            "value": 12073,
            "unit": "ns",
            "extra": "gctime=0\nmemory=15736\nallocs=173\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/single Σ_i(σ_i*σ_j)",
            "value": 10410,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12392\nallocs=139\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=2",
            "value": 2464.5555555555557,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2928\nallocs=59\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=3",
            "value": 4674.428571428572,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5248\nallocs=119\nparams={\"evals\":7,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=4",
            "value": 7680.625,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8448\nallocs=201\nparams={\"evals\":4,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=5",
            "value": 11471,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12528\nallocs=305\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Ground state/3-level rewrite",
            "value": 2791.8888888888887,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3088\nallocs=44\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Multi-mode/2-mode 6-op chain",
            "value": 5894.333333333333,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6608\nallocs=154\nparams={\"evals\":6,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H",
            "value": 4017.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4048\nallocs=82\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H²",
            "value": 896555,
            "unit": "ns",
            "extra": "gctime=0\nmemory=328128\nallocs=6969\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H",
            "value": 4349.571428571428,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4512\nallocs=95\nparams={\"evals\":7,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H²",
            "value": 1338521,
            "unit": "ns",
            "extra": "gctime=0\nmemory=521512\nallocs=11578\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H",
            "value": 99846,
            "unit": "ns",
            "extra": "gctime=0\nmemory=47216\nallocs=901\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H²",
            "value": 1410015,
            "unit": "ns",
            "extra": "gctime=0\nmemory=603104\nallocs=12758\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
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
          "id": "7b7310dffcffa63113d3c0891cca0457d72c2925",
          "message": "perf: Native numeric and parameter-polynomial coefficients (#183)\n\n* Tier 1: native-numeric coefficient fast path\n\nReplace the `CNum = Complex{Num}` coefficient type with a concrete `Coeff`\nstruct: an inline `ComplexF64` fast path plus a `Complex{Num}` symbolic tail.\nConcrete numeric coefficients (the common case) now use machine arithmetic and\nnever touch SymbolicUtils hashconsing, while genuinely symbolic coefficients\ndelegate to Symbolics through the tail.\n\nDetails:\n- `Coeff(z::ComplexF64, tail::Union{Nothing, Complex{Num}})`; `tail === nothing`\n  marks the native fast path. Native `_mul_cnum`/`_add_cnum`/etc. are\n  allocation-free (measured 0 B/op vs ~2.7 kB/op for `Complex{Num}`).\n- Exactness gate (`_native_complex`): values that do not round-trip through\n  `ComplexF64` (e.g. `1//3`, bignums) stay symbolic, so no silent truncation.\n- Const folding (`_numeric_value`): `substitute`/`simplify` results that reduce\n  to a number collapse back to the native path.\n- Signed-zero normalization in `_native` keeps structurally equal coefficients\n  `isequal` and identically hashed (e.g. `conj(2)` vs `2`).\n- Boundary materialization via `to_num`/`real`/`imag`/`_to_complex`; public\n  `prefactor`/`getindex` keep returning `Complex{Num}`.\n- `Coeff` supports scalar arithmetic and mixed comparison with `Number` so it\n  stays a drop-in for downstream code that iterates coefficients.\n\nFull test suite green (2636 tests, incl. Aqua); adds test/coeff_test.jl.\n\n* Tier 2: native parameter-polynomial coefficients\n\nAdd a sparse parameter-polynomial (Poly of Monomial) as a third Coeff\nform alongside the native ComplexF64 and Complex{Num} tails. Products and\nsums of named parameters now stay native: factor identity uses objectid\nand ===, exploiting SymbolicUtils hashconsing so structurally equal\nsymbols are the same object. Scalars are ComplexF64 and exponents Int, so\nthe common symbolic coefficient never routes through hashconsing.\n\nRecognition is a total Coeff-algebra interpreter (_rec) with no nothing\nsentinels; an unrecognized subterm becomes a symbolic-tail Coeff that the\narithmetic absorbs. The native fast path is marked by a zero-size Native\nsingleton rather than nothing.\n\nPolynomial arithmetic is fully type-stable (verified zero dispatch via\n@report_opt). _term_to_num builds powers by repeated Num multiplication\nbecause Num^Int infers to Any in SymbolicUtils, which otherwise degraded\nto_num and real(::Coeff) throughout.\n\nBenchmarks vs base: numeric power 2.1x, symbolic H^4 2.65x, many-mode H^2\n(M=8) 3.54x, nested commutator 2.36x. Tests: 2454 functional plus 583\nquality (Aqua, ExplicitImports, JET) green.\n\nNote: factor canonicalization relies on SymbolicUtils hashconsing being\nenabled (the default).\n\n* format\n\n* fix benchmark and docs CI\n\n* fix change perf change_index regression\n\n* fix regressed printing and test\n\n* implement code review",
          "timestamp": "2026-06-21T15:48:35+02:00",
          "tree_id": "dba7e532e52794d18af10491dbc129a4896b3072",
          "url": "https://github.com/qojulia/SecondQuantizedAlgebra.jl/commit/7b7310dffcffa63113d3c0891cca0457d72c2925"
        },
        "date": 1782050119643,
        "tool": "julia",
        "benches": [
          {
            "name": "Commutator/Nested JC/depth=1",
            "value": 10446,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23200\nallocs=232\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=2",
            "value": 85567,
            "unit": "ns",
            "extra": "gctime=0\nmemory=143712\nallocs=1575\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=3",
            "value": 264030,
            "unit": "ns",
            "extra": "gctime=0\nmemory=379712\nallocs=4467\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=4",
            "value": 589170,
            "unit": "ns",
            "extra": "gctime=0\nmemory=753776\nallocs=9344\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=5",
            "value": 1146092,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1312208\nallocs=16968\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, V]",
            "value": 234727,
            "unit": "ns",
            "extra": "gctime=0\nmemory=97280\nallocs=1936\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, [S, H0]]",
            "value": 476403,
            "unit": "ns",
            "extra": "gctime=0\nmemory=185264\nallocs=3624\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_Dicke, S_j]",
            "value": 97464,
            "unit": "ns",
            "extra": "gctime=0\nmemory=85136\nallocs=1305\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_JC, σ_j]",
            "value": 98576,
            "unit": "ns",
            "extra": "gctime=0\nmemory=81952\nallocs=1137\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/double-sum spin-spin",
            "value": 18888,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9632\nallocs=148\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/indexed JC H",
            "value": 35212,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16000\nallocs=306\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/double Σ_ij(J_ij*S_i*S_j)",
            "value": 12448,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16504\nallocs=176\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/single Σ_i(σ_i*σ_j)",
            "value": 11577,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13048\nallocs=145\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=2",
            "value": 2362.3333333333335,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3168\nallocs=59\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=3",
            "value": 4452.285714285715,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5648\nallocs=119\nparams={\"evals\":7,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=4",
            "value": 7383.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9072\nallocs=201\nparams={\"evals\":4,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=5",
            "value": 10996,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13440\nallocs=305\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Ground state/3-level rewrite",
            "value": 1228.8,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2752\nallocs=29\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Multi-mode/2-mode 6-op chain",
            "value": 5518.166666666667,
            "unit": "ns",
            "extra": "gctime=0\nmemory=7104\nallocs=154\nparams={\"evals\":6,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H",
            "value": 3318.625,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4368\nallocs=82\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H²",
            "value": 12529,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12096\nallocs=291\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H",
            "value": 3666.625,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4864\nallocs=95\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H²",
            "value": 24015,
            "unit": "ns",
            "extra": "gctime=0\nmemory=27984\nallocs=632\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H",
            "value": 3711.75,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4384\nallocs=89\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H²",
            "value": 8172,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8384\nallocs=203\nparams={\"evals\":4,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
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
          "id": "4f05f2032ebfa94aaf5e32b2d6a3063881f19a7d",
          "message": "perf: Cached / cheaper QTerm hashing (#184)\n\n* Cached / cheaper QTerm hashing\n\n* add changelog",
          "timestamp": "2026-06-21T21:48:04+02:00",
          "tree_id": "f0672f6b2918d37d18457cae574773dcfc04fcae",
          "url": "https://github.com/qojulia/SecondQuantizedAlgebra.jl/commit/4f05f2032ebfa94aaf5e32b2d6a3063881f19a7d"
        },
        "date": 1782071455638,
        "tool": "julia",
        "benches": [
          {
            "name": "Commutator/Nested JC/depth=1",
            "value": 7143.4,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22432\nallocs=150\nparams={\"evals\":5,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=2",
            "value": 61105,
            "unit": "ns",
            "extra": "gctime=0\nmemory=130720\nallocs=956\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=3",
            "value": 191669,
            "unit": "ns",
            "extra": "gctime=0\nmemory=333904\nallocs=2628\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=4",
            "value": 438850,
            "unit": "ns",
            "extra": "gctime=0\nmemory=645408\nallocs=5381\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=5",
            "value": 861280,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1099104\nallocs=9626\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, V]",
            "value": 257882,
            "unit": "ns",
            "extra": "gctime=0\nmemory=92048\nallocs=1707\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, [S, H0]]",
            "value": 522845.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=173696\nallocs=3147\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_Dicke, S_j]",
            "value": 88356,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75824\nallocs=989\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_JC, σ_j]",
            "value": 102802,
            "unit": "ns",
            "extra": "gctime=0\nmemory=78144\nallocs=949\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/double-sum spin-spin",
            "value": 18555,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9392\nallocs=137\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/indexed JC H",
            "value": 33973,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14608\nallocs=255\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/double Σ_ij(J_ij*S_i*S_j)",
            "value": 12103,
            "unit": "ns",
            "extra": "gctime=0\nmemory=15832\nallocs=146\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/single Σ_i(σ_i*σ_j)",
            "value": 13536,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14304\nallocs=144\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=2",
            "value": 1033.9333333333334,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1520\nallocs=21\nparams={\"evals\":15,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=3",
            "value": 1760.3,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1808\nallocs=34\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=4",
            "value": 2579.3333333333335,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2176\nallocs=51\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=5",
            "value": 4238.5625,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2624\nallocs=72\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Ground state/3-level rewrite",
            "value": 1001.4347826086956,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2912\nallocs=25\nparams={\"evals\":23,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Multi-mode/2-mode 6-op chain",
            "value": 2535.777777777778,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1952\nallocs=41\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H",
            "value": 1509.9,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2880\nallocs=29\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H²",
            "value": 4866.285714285715,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4192\nallocs=71\nparams={\"evals\":7,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H",
            "value": 1647.1,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2944\nallocs=33\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H²",
            "value": 8492.333333333334,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12320\nallocs=129\nparams={\"evals\":3,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H",
            "value": 1777.4,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3168\nallocs=35\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H²",
            "value": 3196,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4000\nallocs=59\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
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
          "id": "10bf6e6eb1a08ec9e9be94c5faf5dd7284333373",
          "message": "fix: promote_symtype averages (#185)",
          "timestamp": "2026-06-21T22:14:35+02:00",
          "tree_id": "17ec8519a067572e68b2dae5571fcb1fc7d631ee",
          "url": "https://github.com/qojulia/SecondQuantizedAlgebra.jl/commit/10bf6e6eb1a08ec9e9be94c5faf5dd7284333373"
        },
        "date": 1782073053564,
        "tool": "julia",
        "benches": [
          {
            "name": "Commutator/Nested JC/depth=1",
            "value": 7511.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22432\nallocs=150\nparams={\"evals\":4,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=2",
            "value": 64926.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=130720\nallocs=956\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=3",
            "value": 197500.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=333904\nallocs=2628\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=4",
            "value": 450876,
            "unit": "ns",
            "extra": "gctime=0\nmemory=645408\nallocs=5381\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=5",
            "value": 888967,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1099104\nallocs=9626\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, V]",
            "value": 250109,
            "unit": "ns",
            "extra": "gctime=0\nmemory=92048\nallocs=1707\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, [S, H0]]",
            "value": 513994,
            "unit": "ns",
            "extra": "gctime=0\nmemory=173696\nallocs=3147\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_Dicke, S_j]",
            "value": 87104,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75824\nallocs=989\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_JC, σ_j]",
            "value": 103614,
            "unit": "ns",
            "extra": "gctime=0\nmemory=78144\nallocs=949\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/double-sum spin-spin",
            "value": 19817,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9392\nallocs=137\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/indexed JC H",
            "value": 36218,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14608\nallocs=255\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/double Σ_ij(J_ij*S_i*S_j)",
            "value": 11862,
            "unit": "ns",
            "extra": "gctime=0\nmemory=15832\nallocs=146\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/single Σ_i(σ_i*σ_j)",
            "value": 13875,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14304\nallocs=144\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=2",
            "value": 1101,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1520\nallocs=21\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=3",
            "value": 2179,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1808\nallocs=34\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=4",
            "value": 2676.1111111111113,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2176\nallocs=51\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=5",
            "value": 4264.1875,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2624\nallocs=72\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Ground state/3-level rewrite",
            "value": 1067,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2912\nallocs=25\nparams={\"evals\":14,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Multi-mode/2-mode 6-op chain",
            "value": 2042.8,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1952\nallocs=41\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H",
            "value": 1523.8,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2880\nallocs=29\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H²",
            "value": 4618.714285714285,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4192\nallocs=71\nparams={\"evals\":7,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H",
            "value": 1650.1,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2944\nallocs=33\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H²",
            "value": 7576.75,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12320\nallocs=129\nparams={\"evals\":4,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H",
            "value": 1872.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3168\nallocs=35\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H²",
            "value": 3455.25,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4000\nallocs=59\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
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
          "id": "725b5437b5d87289ad2fc0a5ad9e0fe967573927",
          "message": "build: bump patch (#186)",
          "timestamp": "2026-06-21T22:16:13+02:00",
          "tree_id": "1b906f9a54917654b3fd115f8439329a44e7b0b2",
          "url": "https://github.com/qojulia/SecondQuantizedAlgebra.jl/commit/725b5437b5d87289ad2fc0a5ad9e0fe967573927"
        },
        "date": 1782073158010,
        "tool": "julia",
        "benches": [
          {
            "name": "Commutator/Nested JC/depth=1",
            "value": 7869.25,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22432\nallocs=150\nparams={\"evals\":4,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=2",
            "value": 63925,
            "unit": "ns",
            "extra": "gctime=0\nmemory=130720\nallocs=956\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=3",
            "value": 191543.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=333904\nallocs=2628\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=4",
            "value": 434844,
            "unit": "ns",
            "extra": "gctime=0\nmemory=645408\nallocs=5381\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=5",
            "value": 837376,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1099104\nallocs=9626\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, V]",
            "value": 225354,
            "unit": "ns",
            "extra": "gctime=0\nmemory=92048\nallocs=1707\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, [S, H0]]",
            "value": 451960,
            "unit": "ns",
            "extra": "gctime=0\nmemory=173696\nallocs=3147\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_Dicke, S_j]",
            "value": 83915,
            "unit": "ns",
            "extra": "gctime=0\nmemory=75824\nallocs=989\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_JC, σ_j]",
            "value": 91025,
            "unit": "ns",
            "extra": "gctime=0\nmemory=78144\nallocs=949\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/double-sum spin-spin",
            "value": 18568,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9392\nallocs=137\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/indexed JC H",
            "value": 34071,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14608\nallocs=255\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/double Σ_ij(J_ij*S_i*S_j)",
            "value": 11587,
            "unit": "ns",
            "extra": "gctime=0\nmemory=15832\nallocs=146\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/single Σ_i(σ_i*σ_j)",
            "value": 12589,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14304\nallocs=144\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=2",
            "value": 1031.5172413793102,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1520\nallocs=21\nparams={\"evals\":29,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=3",
            "value": 2250.3,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1808\nallocs=34\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=4",
            "value": 2537.1111111111113,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2176\nallocs=51\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=5",
            "value": 4231.25,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2624\nallocs=72\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Ground state/3-level rewrite",
            "value": 1157.4285714285713,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2912\nallocs=25\nparams={\"evals\":14,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Multi-mode/2-mode 6-op chain",
            "value": 2272.222222222222,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1952\nallocs=41\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H",
            "value": 1589.4,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2880\nallocs=29\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H²",
            "value": 4582.571428571428,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4192\nallocs=71\nparams={\"evals\":7,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H",
            "value": 1672.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2944\nallocs=33\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H²",
            "value": 8910,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12320\nallocs=129\nparams={\"evals\":3,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H",
            "value": 1948.9,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3168\nallocs=35\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H²",
            "value": 3317.375,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4000\nallocs=59\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
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
          "id": "7c2d8fce3cb463139beaace824715cc871a504ab",
          "message": "perf(breaking): the operator types collapsed into one concrete (#187)\n\n* perf(breaking): the seven operator types collapsed into one concrete `Op`.\n\n* Parameter-polynomial coefficient arithmetic\n\n* remove redundant poly add\n\n* fix ci",
          "timestamp": "2026-06-22T14:42:28+02:00",
          "tree_id": "aa2a8be57b1ef5aa5fee37e734808f421cf997af",
          "url": "https://github.com/qojulia/SecondQuantizedAlgebra.jl/commit/7c2d8fce3cb463139beaace824715cc871a504ab"
        },
        "date": 1782132451998,
        "tool": "julia",
        "benches": [
          {
            "name": "Commutator/Nested JC/depth=1",
            "value": 6671.1,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22976\nallocs=118\nparams={\"evals\":5,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=2",
            "value": 45052.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=120928\nallocs=651\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=3",
            "value": 131097,
            "unit": "ns",
            "extra": "gctime=0\nmemory=308864\nallocs=1698\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=4",
            "value": 288574.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=597328\nallocs=3353\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=5",
            "value": 559412,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1021040\nallocs=5876\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, V]",
            "value": 220502,
            "unit": "ns",
            "extra": "gctime=0\nmemory=91376\nallocs=1575\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, [S, H0]]",
            "value": 443157,
            "unit": "ns",
            "extra": "gctime=0\nmemory=174448\nallocs=2904\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_Dicke, S_j]",
            "value": 73781,
            "unit": "ns",
            "extra": "gctime=0\nmemory=71696\nallocs=711\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_JC, σ_j]",
            "value": 84757,
            "unit": "ns",
            "extra": "gctime=0\nmemory=72000\nallocs=731\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/double-sum spin-spin",
            "value": 17597,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9328\nallocs=127\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/indexed JC H",
            "value": 32338,
            "unit": "ns",
            "extra": "gctime=0\nmemory=15216\nallocs=239\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/double Σ_ij(J_ij*S_i*S_j)",
            "value": 10516,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14856\nallocs=104\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/single Σ_i(σ_i*σ_j)",
            "value": 12108.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12960\nallocs=109\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=2",
            "value": 903.5873015873016,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1760\nallocs=8\nparams={\"evals\":63,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=3",
            "value": 1391.1,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2256\nallocs=9\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=4",
            "value": 2811.3,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2928\nallocs=10\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=5",
            "value": 2740.8888888888887,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3744\nallocs=11\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Ground state/3-level rewrite",
            "value": 2007,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2880\nallocs=16\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Multi-mode/2-mode 6-op chain",
            "value": 1604.4,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2560\nallocs=9\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H",
            "value": 1557.4,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3680\nallocs=17\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H²",
            "value": 3814.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6976\nallocs=27\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H",
            "value": 1686.55,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3840\nallocs=17\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H²",
            "value": 7012.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=17248\nallocs=41\nparams={\"evals\":5,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H",
            "value": 1929.9,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4032\nallocs=21\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H²",
            "value": 3377.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5952\nallocs=29\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
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
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "distinct": true,
          "id": "cd8d75acb62a46d02ac303173ae7b5d586ac0d54",
          "message": "perf(breaking): the operator types collapsed into one concrete (#187)\n\n* perf(breaking): the seven operator types collapsed into one concrete `Op`.\n\n* Parameter-polynomial coefficient arithmetic\n\n* remove redundant poly add\n\n* fix ci",
          "timestamp": "2026-06-22T18:53:25+02:00",
          "tree_id": "aeeec071d9e339afc9500f1c2152fa18444f3d42",
          "url": "https://github.com/qojulia/SecondQuantizedAlgebra.jl/commit/cd8d75acb62a46d02ac303173ae7b5d586ac0d54"
        },
        "date": 1782160357037,
        "tool": "julia",
        "benches": [
          {
            "name": "Commutator/Nested JC/depth=1",
            "value": 6338,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22976\nallocs=118\nparams={\"evals\":6,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=2",
            "value": 44397,
            "unit": "ns",
            "extra": "gctime=0\nmemory=122336\nallocs=651\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=3",
            "value": 132751,
            "unit": "ns",
            "extra": "gctime=0\nmemory=314112\nallocs=1698\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=4",
            "value": 294726,
            "unit": "ns",
            "extra": "gctime=0\nmemory=610384\nallocs=3353\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=5",
            "value": 575726.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1049504\nallocs=5876\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, V]",
            "value": 228326,
            "unit": "ns",
            "extra": "gctime=0\nmemory=94192\nallocs=1737\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, [S, H0]]",
            "value": 458218.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=179024\nallocs=3172\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_Dicke, S_j]",
            "value": 78719,
            "unit": "ns",
            "extra": "gctime=0\nmemory=73808\nallocs=825\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_JC, σ_j]",
            "value": 90137,
            "unit": "ns",
            "extra": "gctime=0\nmemory=79008\nallocs=934\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/double-sum spin-spin",
            "value": 18037,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9456\nallocs=135\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/indexed JC H",
            "value": 1684.6,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4048\nallocs=19\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/double Σ_ij(J_ij*S_i*S_j)",
            "value": 10196,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14920\nallocs=108\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/single Σ_i(σ_i*σ_j)",
            "value": 22765,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16992\nallocs=232\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=2",
            "value": 861.1868131868132,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1760\nallocs=8\nparams={\"evals\":91,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=3",
            "value": 1390.1,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2256\nallocs=9\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=4",
            "value": 2750.1499999999996,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2928\nallocs=10\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=5",
            "value": 2678.4444444444443,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3744\nallocs=11\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Ground state/3-level rewrite",
            "value": 1090.84,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2880\nallocs=16\nparams={\"evals\":25,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Multi-mode/2-mode 6-op chain",
            "value": 1565.4,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2560\nallocs=9\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H",
            "value": 1494.75,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3680\nallocs=17\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H²",
            "value": 3678,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6976\nallocs=27\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H",
            "value": 1628.95,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3840\nallocs=17\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H²",
            "value": 6652.2,
            "unit": "ns",
            "extra": "gctime=0\nmemory=17248\nallocs=41\nparams={\"evals\":5,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H",
            "value": 1856.8,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4032\nallocs=21\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H²",
            "value": 3277.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5952\nallocs=29\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
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
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "distinct": true,
          "id": "67d61a7660b98c9f2dc834d620cd7d04ee30f12f",
          "message": "fix: restore type stability and perf in Coeff radical recognition\n\nThe v0.7.2 rational-exponent work introduced two regressions caught by CI on\nmain, both now fixed:\n\n- Type instability (JET): _term_to_num's general fractional branch used\n  Num ^ Rational, which infers as Any (the same reason the integer path uses a\n  repeated-multiply loop), making _term_to_num / to_num return-unstable. Assert\n  ::Num on that one materialization step.\n\n- Benchmark regression (~1.9x on indexed-sum construction): the broadened\n  _is_atom recognized an indexed variable g(i) as a Poly atom, so the index\n  machinery paid a materialize/re-recognize round-trip per index op. Require the\n  call's operation to be a real Function, so g(i) (a callable symbol) stays a\n  symbolic leaf as before, while real/imag/sqrt/exp stay recognized.\n\nAlso restructure the ^ recognizer to narrow the exponent with isa guards\ninstead of calling a helper on an Any value (which forced runtime dispatch\nflagged by JET), and drop float-exponent recognition (only exact Rational\nexponents fold now; sqrt/g^(1//2) are unaffected). Removes _rat_exp.",
          "timestamp": "2026-06-22T23:15:52+02:00",
          "tree_id": "e721f68a425c7c8372e95b80dd0ab752d25dddb5",
          "url": "https://github.com/qojulia/SecondQuantizedAlgebra.jl/commit/67d61a7660b98c9f2dc834d620cd7d04ee30f12f"
        },
        "date": 1782163138661,
        "tool": "julia",
        "benches": [
          {
            "name": "Commutator/Nested JC/depth=1",
            "value": 5785.083333333333,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22976\nallocs=118\nparams={\"evals\":6,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=2",
            "value": 41639,
            "unit": "ns",
            "extra": "gctime=0\nmemory=122336\nallocs=651\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=3",
            "value": 127918,
            "unit": "ns",
            "extra": "gctime=0\nmemory=314112\nallocs=1698\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=4",
            "value": 288838,
            "unit": "ns",
            "extra": "gctime=0\nmemory=610384\nallocs=3353\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=5",
            "value": 572160,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1049392\nallocs=5876\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, V]",
            "value": 233970,
            "unit": "ns",
            "extra": "gctime=0\nmemory=91376\nallocs=1575\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, [S, H0]]",
            "value": 473040,
            "unit": "ns",
            "extra": "gctime=0\nmemory=174448\nallocs=2904\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_Dicke, S_j]",
            "value": 73468,
            "unit": "ns",
            "extra": "gctime=0\nmemory=71792\nallocs=711\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_JC, σ_j]",
            "value": 92113,
            "unit": "ns",
            "extra": "gctime=0\nmemory=72032\nallocs=731\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/double-sum spin-spin",
            "value": 18284,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9328\nallocs=127\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/indexed JC H",
            "value": 33924,
            "unit": "ns",
            "extra": "gctime=0\nmemory=15216\nallocs=239\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/double Σ_ij(J_ij*S_i*S_j)",
            "value": 11101,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14856\nallocs=104\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/single Σ_i(σ_i*σ_j)",
            "value": 11452,
            "unit": "ns",
            "extra": "gctime=0\nmemory=11944\nallocs=96\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=2",
            "value": 795.1709401709402,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1760\nallocs=8\nparams={\"evals\":117,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=3",
            "value": 1269.4,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2256\nallocs=9\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=4",
            "value": 2465.1,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2928\nallocs=10\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=5",
            "value": 2507,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3744\nallocs=11\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Ground state/3-level rewrite",
            "value": 971.8222222222222,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2880\nallocs=16\nparams={\"evals\":45,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Multi-mode/2-mode 6-op chain",
            "value": 1462.8,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2560\nallocs=9\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H",
            "value": 1352.6,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3680\nallocs=17\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H²",
            "value": 3280,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6976\nallocs=27\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H",
            "value": 1428.7,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3840\nallocs=17\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H²",
            "value": 6174.916666666666,
            "unit": "ns",
            "extra": "gctime=0\nmemory=17248\nallocs=41\nparams={\"evals\":6,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H",
            "value": 1734.25,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4032\nallocs=21\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H²",
            "value": 2883.1111111111113,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5952\nallocs=29\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
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
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "distinct": true,
          "id": "b3d15f09a0a1af2f6213889a99f987dd2255433d",
          "message": "fix: restore type stability and perf after rational-exponent Coeff work\n\nThe v0.7.2 rational-exponent work tripped two CI checks on main; both fixed,\nwith the broadened-recognition performance win kept intact.\n\nType stability (JET report_opt):\n- _term_to_num's general fractional branch used Num ^ Rational, which infers as\n  Any (same reason the integer path uses a repeated-multiply loop), making\n  _term_to_num / to_num return-unstable. Assert ::Num on that one step.\n- The ^ recognizer called a helper on the Any-typed exponent, forcing runtime\n  dispatch. Narrow with isa guards instead (Integer / Rational{Int}); drop\n  float-exponent recognition (only exact Rational exponents fold; sqrt and\n  g^(1//2) are unaffected). Removes _rat_exp.\n\nPerformance (Indexing benchmarks):\n- Keep recognizing an indexed variable g(i) as a Poly atom: this makes\n  simplify of indexed Hamiltonians ~25x faster (the broadened-_is_atom win).\n- It regressed sum construction ~3x because _depends_on_index_term materialized\n  the Poly coefficient to a Complex{Num} just to scan its variables. Add a Poly\n  fast path that scans the monomial atoms directly (wrapping each in Num so\n  get_variables stays type-stable), restoring sum-construction to baseline.",
          "timestamp": "2026-06-22T23:37:12+02:00",
          "tree_id": "a76ee8cf6a699f914d1d15b29133cb9783a90648",
          "url": "https://github.com/qojulia/SecondQuantizedAlgebra.jl/commit/b3d15f09a0a1af2f6213889a99f987dd2255433d"
        },
        "date": 1782164415116,
        "tool": "julia",
        "benches": [
          {
            "name": "Commutator/Nested JC/depth=1",
            "value": 5827.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22976\nallocs=118\nparams={\"evals\":6,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=2",
            "value": 41457,
            "unit": "ns",
            "extra": "gctime=0\nmemory=122336\nallocs=651\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=3",
            "value": 127698.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=314112\nallocs=1698\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=4",
            "value": 287256,
            "unit": "ns",
            "extra": "gctime=0\nmemory=610384\nallocs=3353\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=5",
            "value": 560545,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1049392\nallocs=5876\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, V]",
            "value": 240789,
            "unit": "ns",
            "extra": "gctime=0\nmemory=91376\nallocs=1575\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, [S, H0]]",
            "value": 494853,
            "unit": "ns",
            "extra": "gctime=0\nmemory=174448\nallocs=2904\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_Dicke, S_j]",
            "value": 49757.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=64816\nallocs=545\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_JC, σ_j]",
            "value": 70250.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=71264\nallocs=674\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/double-sum spin-spin",
            "value": 19085,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9344\nallocs=128\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/indexed JC H",
            "value": 1486.8,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4048\nallocs=19\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/double Σ_ij(J_ij*S_i*S_j)",
            "value": 10991,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14856\nallocs=104\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/single Σ_i(σ_i*σ_j)",
            "value": 14407,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13776\nallocs=131\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=2",
            "value": 786.697247706422,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1760\nallocs=8\nparams={\"evals\":109,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=3",
            "value": 1833.45,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2256\nallocs=9\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=4",
            "value": 2521.7,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2928\nallocs=10\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=5",
            "value": 2496.8888888888887,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3744\nallocs=11\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Ground state/3-level rewrite",
            "value": 964.1025641025641,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2880\nallocs=16\nparams={\"evals\":39,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Multi-mode/2-mode 6-op chain",
            "value": 1451.7,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2560\nallocs=9\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H",
            "value": 1348.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3680\nallocs=17\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H²",
            "value": 3289.875,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6976\nallocs=27\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H",
            "value": 1415.7,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3840\nallocs=17\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H²",
            "value": 6031.333333333333,
            "unit": "ns",
            "extra": "gctime=0\nmemory=17248\nallocs=41\nparams={\"evals\":6,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H",
            "value": 1694.15,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4032\nallocs=21\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H²",
            "value": 2856.4444444444443,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5952\nallocs=29\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
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
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "distinct": true,
          "id": "98a5c6b096d9dbf4c8f76912d7e0fdc2516d0182",
          "message": "docs: trim verbose comments in the Coeff/Poly tier\n\nCondense the multi-line rationale blocks in cnum.jl and monomial.jl to concise\none-liners. Comment-only; no behavior change.",
          "timestamp": "2026-06-22T23:47:51+02:00",
          "tree_id": "bd97b829a82bd2ea14e6c19d8930e7dc2ca90e31",
          "url": "https://github.com/qojulia/SecondQuantizedAlgebra.jl/commit/98a5c6b096d9dbf4c8f76912d7e0fdc2516d0182"
        },
        "date": 1782165064773,
        "tool": "julia",
        "benches": [
          {
            "name": "Commutator/Nested JC/depth=1",
            "value": 6814.4,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22976\nallocs=118\nparams={\"evals\":5,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=2",
            "value": 46390,
            "unit": "ns",
            "extra": "gctime=0\nmemory=122336\nallocs=651\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=3",
            "value": 139962,
            "unit": "ns",
            "extra": "gctime=0\nmemory=314112\nallocs=1698\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=4",
            "value": 302826,
            "unit": "ns",
            "extra": "gctime=0\nmemory=610384\nallocs=3353\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=5",
            "value": 591552,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1049392\nallocs=5876\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, V]",
            "value": 222230.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=91376\nallocs=1575\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, [S, H0]]",
            "value": 447034,
            "unit": "ns",
            "extra": "gctime=0\nmemory=174448\nallocs=2904\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_Dicke, S_j]",
            "value": 52399,
            "unit": "ns",
            "extra": "gctime=0\nmemory=64816\nallocs=545\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_JC, σ_j]",
            "value": 67171,
            "unit": "ns",
            "extra": "gctime=0\nmemory=71264\nallocs=674\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/double-sum spin-spin",
            "value": 19650,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9344\nallocs=128\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/indexed JC H",
            "value": 2202.3,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4424\nallocs=25\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/double Σ_ij(J_ij*S_i*S_j)",
            "value": 10535,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14856\nallocs=104\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/single Σ_i(σ_i*σ_j)",
            "value": 12869,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13048\nallocs=122\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=2",
            "value": 876.8048780487804,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1760\nallocs=8\nparams={\"evals\":82,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=3",
            "value": 2023.1,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2256\nallocs=9\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=4",
            "value": 2740.1,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2928\nallocs=10\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=5",
            "value": 2763.6111111111113,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3744\nallocs=11\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Ground state/3-level rewrite",
            "value": 1238.75,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2880\nallocs=16\nparams={\"evals\":16,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Multi-mode/2-mode 6-op chain",
            "value": 1611.4,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2560\nallocs=9\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H",
            "value": 1550.85,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3680\nallocs=17\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H²",
            "value": 3746.875,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6976\nallocs=27\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H",
            "value": 1700,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3840\nallocs=17\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H²",
            "value": 6953.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=17248\nallocs=41\nparams={\"evals\":5,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H",
            "value": 2037.1,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4032\nallocs=21\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H²",
            "value": 3290.4444444444443,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5952\nallocs=29\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
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
          "id": "7194c3835c1fc2f4eeda4e7be870900689865df9",
          "message": "Isbit operators (#189)\n\n* perf: make `Op` stack allocated\n\n* format\n\n* Global intern tables for Sym\n\n* fix code review",
          "timestamp": "2026-06-23T10:37:15+02:00",
          "tree_id": "f8ae7c5b7b10a2a7c82626ccbddb8f5d449faa4b",
          "url": "https://github.com/qojulia/SecondQuantizedAlgebra.jl/commit/7194c3835c1fc2f4eeda4e7be870900689865df9"
        },
        "date": 1782204060000,
        "tool": "julia",
        "benches": [
          {
            "name": "Commutator/Nested JC/depth=1",
            "value": 5555.333333333333,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21840\nallocs=117\nparams={\"evals\":6,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=2",
            "value": 40085,
            "unit": "ns",
            "extra": "gctime=0\nmemory=115232\nallocs=650\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=3",
            "value": 121191.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=293504\nallocs=1697\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=4",
            "value": 269759,
            "unit": "ns",
            "extra": "gctime=0\nmemory=566880\nallocs=3352\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=5",
            "value": 530922,
            "unit": "ns",
            "extra": "gctime=0\nmemory=970640\nallocs=5875\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, V]",
            "value": 243154.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=90144\nallocs=1613\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, [S, H0]]",
            "value": 497600,
            "unit": "ns",
            "extra": "gctime=0\nmemory=170816\nallocs=2954\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_Dicke, S_j]",
            "value": 48020,
            "unit": "ns",
            "extra": "gctime=0\nmemory=54880\nallocs=563\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_JC, σ_j]",
            "value": 68127,
            "unit": "ns",
            "extra": "gctime=0\nmemory=62208\nallocs=692\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/double-sum spin-spin",
            "value": 20188,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8752\nallocs=128\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/indexed JC H",
            "value": 1340.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3488\nallocs=19\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/double Σ_ij(J_ij*S_i*S_j)",
            "value": 10139,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12984\nallocs=104\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/single Σ_i(σ_i*σ_j)",
            "value": 12233,
            "unit": "ns",
            "extra": "gctime=0\nmemory=11656\nallocs=125\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=2",
            "value": 694.0396341463414,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1568\nallocs=8\nparams={\"evals\":164,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=3",
            "value": 1056.7058823529412,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1904\nallocs=9\nparams={\"evals\":17,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=4",
            "value": 1436.7,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2352\nallocs=10\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=5",
            "value": 1985.7,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2848\nallocs=11\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Ground state/3-level rewrite",
            "value": 894.9464285714286,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2816\nallocs=16\nparams={\"evals\":84,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Multi-mode/2-mode 6-op chain",
            "value": 1589,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2080\nallocs=9\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H",
            "value": 1259.4,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3264\nallocs=17\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H²",
            "value": 2874.222222222222,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5472\nallocs=27\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H",
            "value": 1445.2,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328\nallocs=17\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H²",
            "value": 6649.916666666666,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14432\nallocs=41\nparams={\"evals\":6,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H",
            "value": 2334.3500000000004,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=21\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H²",
            "value": 2661.6666666666665,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4992\nallocs=29\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
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
          "id": "efef89c40b0577795bd33ad9d3fed84443636d7b",
          "message": "feat: MutableArithmetics accumulator (#190)\n\n* feat: MutableArithmetics accumulator\n\n* add accumaltion docs\n\n* fix changelog\n\n* fix docs\n\n* fix coverage",
          "timestamp": "2026-06-23T15:44:19+02:00",
          "tree_id": "661b59e886ca9081d219c3d4726c7d5920811793",
          "url": "https://github.com/qojulia/SecondQuantizedAlgebra.jl/commit/efef89c40b0577795bd33ad9d3fed84443636d7b"
        },
        "date": 1782222702800,
        "tool": "julia",
        "benches": [
          {
            "name": "Accumulation/Many-mode H/foldl M=16",
            "value": 13054,
            "unit": "ns",
            "extra": "gctime=0\nmemory=59472\nallocs=213\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/foldl M=24",
            "value": 28674,
            "unit": "ns",
            "extra": "gctime=0\nmemory=131408\nallocs=436\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/foldl M=8",
            "value": 3614.5555555555557,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12432\nallocs=63\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=16",
            "value": 3021.6,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6816\nallocs=24\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=24",
            "value": 2767.8888888888887,
            "unit": "ns",
            "extra": "gctime=0\nmemory=7968\nallocs=32\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=8",
            "value": 809.4113475177305,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=13\nparams={\"evals\":141,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Same-site/foldl",
            "value": 7143.4,
            "unit": "ns",
            "extra": "gctime=0\nmemory=30912\nallocs=138\nparams={\"evals\":5,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Same-site/sum",
            "value": 3785.3888888888887,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4576\nallocs=29\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=1",
            "value": 5247.642857142857,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21840\nallocs=117\nparams={\"evals\":7,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=2",
            "value": 38792,
            "unit": "ns",
            "extra": "gctime=0\nmemory=115232\nallocs=650\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=3",
            "value": 118500,
            "unit": "ns",
            "extra": "gctime=0\nmemory=293504\nallocs=1697\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=4",
            "value": 264502,
            "unit": "ns",
            "extra": "gctime=0\nmemory=566880\nallocs=3352\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=5",
            "value": 521881,
            "unit": "ns",
            "extra": "gctime=0\nmemory=970496\nallocs=5875\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, V]",
            "value": 239786,
            "unit": "ns",
            "extra": "gctime=0\nmemory=88352\nallocs=1552\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, [S, H0]]",
            "value": 496574,
            "unit": "ns",
            "extra": "gctime=0\nmemory=167840\nallocs=2861\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_Dicke, S_j]",
            "value": 47648,
            "unit": "ns",
            "extra": "gctime=0\nmemory=54304\nallocs=545\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_JC, σ_j]",
            "value": 68097,
            "unit": "ns",
            "extra": "gctime=0\nmemory=61376\nallocs=658\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/double-sum spin-spin",
            "value": 20318,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8688\nallocs=124\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/indexed JC H",
            "value": 1295.4,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3488\nallocs=19\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/double Σ_ij(J_ij*S_i*S_j)",
            "value": 9828,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12936\nallocs=101\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/single Σ_i(σ_i*σ_j)",
            "value": 13485,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12224\nallocs=127\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=2",
            "value": 655.4416666666666,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1568\nallocs=8\nparams={\"evals\":180,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=3",
            "value": 932.7586206896551,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1904\nallocs=9\nparams={\"evals\":58,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=4",
            "value": 1287.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2352\nallocs=10\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=5",
            "value": 2360.9,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2848\nallocs=11\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Ground state/3-level rewrite",
            "value": 857.8352941176471,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2816\nallocs=16\nparams={\"evals\":85,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Multi-mode/2-mode 6-op chain",
            "value": 1109.9285714285713,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2080\nallocs=9\nparams={\"evals\":14,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H",
            "value": 1223.3,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3264\nallocs=17\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H²",
            "value": 2752.9444444444443,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5472\nallocs=27\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H",
            "value": 2184,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328\nallocs=17\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H²",
            "value": 5061.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14432\nallocs=41\nparams={\"evals\":7,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H",
            "value": 2283.7,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=21\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H²",
            "value": 2537,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4992\nallocs=29\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
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
          "id": "8d905d64ff4fd66e0451327f07899dd6d95d063e",
          "message": "build(deps): bump actions/checkout from 6 to 7 (#191)\n\nBumps [actions/checkout](https://github.com/actions/checkout) from 6 to 7.\n- [Release notes](https://github.com/actions/checkout/releases)\n- [Changelog](https://github.com/actions/checkout/blob/main/CHANGELOG.md)\n- [Commits](https://github.com/actions/checkout/compare/v6...v7)\n\n---\nupdated-dependencies:\n- dependency-name: actions/checkout\n  dependency-version: '7'\n  dependency-type: direct:production\n  update-type: version-update:semver-major\n...\n\nSigned-off-by: dependabot[bot] <support@github.com>\nCo-authored-by: dependabot[bot] <49699333+dependabot[bot]@users.noreply.github.com>",
          "timestamp": "2026-06-23T23:20:23+02:00",
          "tree_id": "05d18086a8131592334058d86cd1d1ef3b935cb5",
          "url": "https://github.com/qojulia/SecondQuantizedAlgebra.jl/commit/8d905d64ff4fd66e0451327f07899dd6d95d063e"
        },
        "date": 1782250205154,
        "tool": "julia",
        "benches": [
          {
            "name": "Accumulation/Many-mode H/foldl M=16",
            "value": 14371.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=59472\nallocs=213\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/foldl M=24",
            "value": 29893.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=131408\nallocs=436\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/foldl M=8",
            "value": 3839.0555555555557,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12432\nallocs=63\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=16",
            "value": 3028.3,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6816\nallocs=24\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=24",
            "value": 3617.5555555555557,
            "unit": "ns",
            "extra": "gctime=0\nmemory=7968\nallocs=32\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=8",
            "value": 903.9568965517242,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=13\nparams={\"evals\":116,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Same-site/foldl",
            "value": 8520.375,
            "unit": "ns",
            "extra": "gctime=0\nmemory=30912\nallocs=138\nparams={\"evals\":4,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Same-site/sum",
            "value": 4090.277777777778,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4576\nallocs=29\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=1",
            "value": 5328.666666666667,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21840\nallocs=117\nparams={\"evals\":6,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=2",
            "value": 37634.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=115232\nallocs=650\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=3",
            "value": 112194.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=293504\nallocs=1697\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=4",
            "value": 256552,
            "unit": "ns",
            "extra": "gctime=0\nmemory=566880\nallocs=3352\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=5",
            "value": 505272,
            "unit": "ns",
            "extra": "gctime=0\nmemory=970704\nallocs=5875\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, V]",
            "value": 212196,
            "unit": "ns",
            "extra": "gctime=0\nmemory=88928\nallocs=1575\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, [S, H0]]",
            "value": 421713,
            "unit": "ns",
            "extra": "gctime=0\nmemory=169216\nallocs=2904\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_Dicke, S_j]",
            "value": 45068,
            "unit": "ns",
            "extra": "gctime=0\nmemory=54304\nallocs=545\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_JC, σ_j]",
            "value": 59667,
            "unit": "ns",
            "extra": "gctime=0\nmemory=61632\nallocs=674\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/double-sum spin-spin",
            "value": 20019.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8816\nallocs=130\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/indexed JC H",
            "value": 1340.35,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3488\nallocs=19\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/double Σ_ij(J_ij*S_i*S_j)",
            "value": 8444.333333333334,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12984\nallocs=104\nparams={\"evals\":3,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/single Σ_i(σ_i*σ_j)",
            "value": 10818,
            "unit": "ns",
            "extra": "gctime=0\nmemory=11560\nallocs=122\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=2",
            "value": 661.1988636363636,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1568\nallocs=8\nparams={\"evals\":176,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=3",
            "value": 898.0526315789473,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1904\nallocs=9\nparams={\"evals\":57,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=4",
            "value": 1223.7,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2352\nallocs=10\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=5",
            "value": 2004.45,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2848\nallocs=11\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Ground state/3-level rewrite",
            "value": 893.5138888888889,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2816\nallocs=16\nparams={\"evals\":72,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Multi-mode/2-mode 6-op chain",
            "value": 1045.0882352941176,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2080\nallocs=9\nparams={\"evals\":17,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H",
            "value": 1228.85,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3264\nallocs=17\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H²",
            "value": 2644.222222222222,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5472\nallocs=27\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H",
            "value": 1837.15,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328\nallocs=17\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H²",
            "value": 4984.071428571428,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14432\nallocs=41\nparams={\"evals\":7,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H",
            "value": 2082.3,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=21\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H²",
            "value": 2531.7777777777774,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4992\nallocs=29\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
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
          "id": "cbecd73ddd4476d82b289246f1e053a08b140323",
          "message": "perf: stream QAdd commutators into one dict (no-index fast path) (#192)\n\n* perf: stream QAdd commutators into one dict (no-index fast path)\n\ncommutator(::QAdd, ::QSym) and commutator(::QAdd, ::QAdd) on operators\nwithout summation indices now emit +ta·tb and -tb·ta straight into one\nshared result dict via _emit_product!, instead of wrapping each term of\nthe left operand in a temporary single-term QAdd and building a*b, b*a,\nand their difference as separate QAdds before merging. Operators that\ncarry summation indices keep the existing path (its per-term\n_absorb_pinned_sums index-scope bookkeeping is unchanged).\n\nByte-identical to a*b - b*a across Fock/NLevel/Pauli/Spin operands and\nself-commutator residuals. Measured on Σ_{i,j} g aᵢ† aⱼ: QAdd×QSym about\n1.9x-2.2x faster, QAdd×QAdd about 2.2x faster, both about 5x less memory.\n\nBumps to v0.8.2 with a changelog entry.\n\n* try to make benchmarks robust",
          "timestamp": "2026-06-25T09:53:39+02:00",
          "tree_id": "5488f4430bf70893dcd6fd668c609c8ed118a648",
          "url": "https://github.com/qojulia/SecondQuantizedAlgebra.jl/commit/cbecd73ddd4476d82b289246f1e053a08b140323"
        },
        "date": 1782374626276,
        "tool": "julia",
        "benches": [
          {
            "name": "Accumulation/Many-mode H/foldl M=16",
            "value": 11631,
            "unit": "ns",
            "extra": "gctime=0\nmemory=59472\nallocs=213\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/foldl M=24",
            "value": 26079,
            "unit": "ns",
            "extra": "gctime=0\nmemory=131408\nallocs=436\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/foldl M=8",
            "value": 2929.8888888888887,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12432\nallocs=63\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=16",
            "value": 1715.2,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6816\nallocs=24\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=24",
            "value": 2370,
            "unit": "ns",
            "extra": "gctime=0\nmemory=7968\nallocs=32\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=8",
            "value": 758.3879310344828,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=13\nparams={\"evals\":116,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Same-site/foldl",
            "value": 6354,
            "unit": "ns",
            "extra": "gctime=0\nmemory=30912\nallocs=138\nparams={\"evals\":5,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Same-site/sum",
            "value": 2569.222222222222,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4576\nallocs=29\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=1",
            "value": 1673.1,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4128\nallocs=39\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=2",
            "value": 18425,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26912\nallocs=280\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=3",
            "value": 67535,
            "unit": "ns",
            "extra": "gctime=0\nmemory=84464\nallocs=839\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=4",
            "value": 169186,
            "unit": "ns",
            "extra": "gctime=0\nmemory=184784\nallocs=1806\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=5",
            "value": 378015,
            "unit": "ns",
            "extra": "gctime=0\nmemory=361616\nallocs=3440\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, V]",
            "value": 201656,
            "unit": "ns",
            "extra": "gctime=0\nmemory=55248\nallocs=1283\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, [S, H0]]",
            "value": 366362,
            "unit": "ns",
            "extra": "gctime=0\nmemory=95872\nallocs=2167\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_Dicke, S_j]",
            "value": 44753,
            "unit": "ns",
            "extra": "gctime=0\nmemory=54304\nallocs=545\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_JC, σ_j]",
            "value": 64280,
            "unit": "ns",
            "extra": "gctime=0\nmemory=61632\nallocs=674\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/double-sum spin-spin",
            "value": 16260,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8800\nallocs=129\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/indexed JC H",
            "value": 1203.3,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3488\nallocs=19\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/double Σ_ij(J_ij*S_i*S_j)",
            "value": 8726,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12984\nallocs=104\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/single Σ_i(σ_i*σ_j)",
            "value": 11001,
            "unit": "ns",
            "extra": "gctime=0\nmemory=11560\nallocs=122\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=2",
            "value": 608.4943820224719,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1568\nallocs=8\nparams={\"evals\":178,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=3",
            "value": 877.4583333333334,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1904\nallocs=9\nparams={\"evals\":48,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=4",
            "value": 1188.2,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2352\nallocs=10\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=5",
            "value": 1582,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2848\nallocs=11\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Ground state/3-level rewrite",
            "value": 818.4516129032259,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2816\nallocs=16\nparams={\"evals\":62,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Multi-mode/2-mode 6-op chain",
            "value": 963.8,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2080\nallocs=9\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H",
            "value": 1131,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3264\nallocs=17\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H²",
            "value": 2524.6666666666665,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5472\nallocs=27\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H",
            "value": 1166.1,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328\nallocs=17\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H²",
            "value": 4541.428571428572,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14432\nallocs=41\nparams={\"evals\":7,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H",
            "value": 1469.8,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=21\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H²",
            "value": 2434.5555555555557,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4992\nallocs=29\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
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
          "id": "12976da14300a7535912e34533ae3b2cff3190fd",
          "message": "fix: fold conj(conj(x)) and conj of a constant in coefficient algebra (#193)\n\n* fix: fold conj(conj(x)) and conj of a constant in coefficient algebra\n\nAdjoint and coefficient recognition left double conjugates and conj of numeric\nconstants unfolded:\n\n* `qadjoint(conj(x))` built `conj(conj(x))` instead of `x`. `inner_adjoint`\n  already folded this case; `qadjoint` did not.\n* the operator-level coefficient conjugation (`_conj_atom` / `_sym_conj`)\n  wrapped an already-conjugated atom in a second `conj`.\n* `_rec` kept `conj(<number>)` symbolic, e.g. `conj(0.0)` left after substituting\n  a complex parameter to a real value.\n\nThese unfolded forms never simplify, so a coefficient that is mathematically\nzero (for instance a complex coupling substituted to 0) stays in the expression\nand survives downstream as dead time-dependent terms that get re-evaluated every\nODE step after translation.\n\nFold `conj(conj(x))` to `x` in `qadjoint` and the atom conjugation, and fold\n`conj` of a constant in `_rec`. `conj` of a symbolic variable is unchanged, and\ninvolution and Hermiticity are preserved. Full test suite passes (2830 tests).\n\n* add regression test\n\n* fix changelog",
          "timestamp": "2026-06-25T10:07:19+02:00",
          "tree_id": "56fd87c41c2028345e470423abb1dda711618fa4",
          "url": "https://github.com/qojulia/SecondQuantizedAlgebra.jl/commit/12976da14300a7535912e34533ae3b2cff3190fd"
        },
        "date": 1782375047281,
        "tool": "julia",
        "benches": [
          {
            "name": "Accumulation/Many-mode H/foldl M=16",
            "value": 11647,
            "unit": "ns",
            "extra": "gctime=0\nmemory=59472\nallocs=213\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/foldl M=24",
            "value": 26189,
            "unit": "ns",
            "extra": "gctime=0\nmemory=131408\nallocs=436\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/foldl M=8",
            "value": 2891.75,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12432\nallocs=63\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=16",
            "value": 1827.7,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6816\nallocs=24\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=24",
            "value": 2449.222222222222,
            "unit": "ns",
            "extra": "gctime=0\nmemory=7968\nallocs=32\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=8",
            "value": 736.5625,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=13\nparams={\"evals\":128,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Same-site/foldl",
            "value": 6872.75,
            "unit": "ns",
            "extra": "gctime=0\nmemory=30912\nallocs=138\nparams={\"evals\":4,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Same-site/sum",
            "value": 2385.777777777778,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4576\nallocs=29\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=1",
            "value": 1823.7,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4128\nallocs=39\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=2",
            "value": 16875,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26544\nallocs=276\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=3",
            "value": 62332,
            "unit": "ns",
            "extra": "gctime=0\nmemory=84096\nallocs=835\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=4",
            "value": 161068,
            "unit": "ns",
            "extra": "gctime=0\nmemory=187808\nallocs=1805\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=5",
            "value": 356158,
            "unit": "ns",
            "extra": "gctime=0\nmemory=365408\nallocs=3447\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, V]",
            "value": 181319,
            "unit": "ns",
            "extra": "gctime=0\nmemory=55248\nallocs=1283\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, [S, H0]]",
            "value": 329659,
            "unit": "ns",
            "extra": "gctime=0\nmemory=95872\nallocs=2167\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_Dicke, S_j]",
            "value": 44155,
            "unit": "ns",
            "extra": "gctime=0\nmemory=54304\nallocs=545\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_JC, σ_j]",
            "value": 58537,
            "unit": "ns",
            "extra": "gctime=0\nmemory=61632\nallocs=674\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/double-sum spin-spin",
            "value": 14442,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8800\nallocs=129\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/indexed JC H",
            "value": 1292.9,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3488\nallocs=19\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/double Σ_ij(J_ij*S_i*S_j)",
            "value": 8546,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12984\nallocs=104\nparams={\"evals\":3,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/single Σ_i(σ_i*σ_j)",
            "value": 9915,
            "unit": "ns",
            "extra": "gctime=0\nmemory=11560\nallocs=122\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=2",
            "value": 659.0736196319018,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1568\nallocs=8\nparams={\"evals\":163,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=3",
            "value": 925.3333333333334,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1904\nallocs=9\nparams={\"evals\":15,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=4",
            "value": 1279.9,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2352\nallocs=10\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=5",
            "value": 1707.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2848\nallocs=11\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Ground state/3-level rewrite",
            "value": 866.7547169811321,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2816\nallocs=16\nparams={\"evals\":53,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Multi-mode/2-mode 6-op chain",
            "value": 1055.6,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2080\nallocs=9\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H",
            "value": 1176.8,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3264\nallocs=17\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H²",
            "value": 2659.4444444444443,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5472\nallocs=27\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H",
            "value": 1244.8,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328\nallocs=17\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H²",
            "value": 4992.333333333333,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14432\nallocs=41\nparams={\"evals\":6,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H",
            "value": 1495.2,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=21\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H²",
            "value": 2552.6666666666665,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4992\nallocs=29\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
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
          "id": "f10f17875a660a0fbb99de7ec5eb5f75118e0915",
          "message": "perf: o(n) accumlation in average (#196)\n\n* perf: o(n) accumlation in average\n\n* also do prod",
          "timestamp": "2026-06-26T10:04:25+02:00",
          "tree_id": "660f54454b7e6faa2e0fe0faeb83b4987a811e8d",
          "url": "https://github.com/qojulia/SecondQuantizedAlgebra.jl/commit/f10f17875a660a0fbb99de7ec5eb5f75118e0915"
        },
        "date": 1782461641405,
        "tool": "julia",
        "benches": [
          {
            "name": "Accumulation/Many-mode H/foldl M=16",
            "value": 11206,
            "unit": "ns",
            "extra": "gctime=0\nmemory=59472\nallocs=213\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/foldl M=24",
            "value": 24691,
            "unit": "ns",
            "extra": "gctime=0\nmemory=131408\nallocs=436\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/foldl M=8",
            "value": 2875.25,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12432\nallocs=63\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=16",
            "value": 1719,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6816\nallocs=24\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=24",
            "value": 2329.222222222222,
            "unit": "ns",
            "extra": "gctime=0\nmemory=7968\nallocs=32\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=8",
            "value": 768.9888888888889,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=13\nparams={\"evals\":90,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Same-site/foldl",
            "value": 7412.25,
            "unit": "ns",
            "extra": "gctime=0\nmemory=30912\nallocs=138\nparams={\"evals\":4,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Same-site/sum",
            "value": 3026.375,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4576\nallocs=29\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=1",
            "value": 1660.2,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4128\nallocs=39\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=2",
            "value": 16968,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26544\nallocs=276\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=3",
            "value": 64432,
            "unit": "ns",
            "extra": "gctime=0\nmemory=83728\nallocs=831\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=4",
            "value": 161105,
            "unit": "ns",
            "extra": "gctime=0\nmemory=183248\nallocs=1790\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=5",
            "value": 356567,
            "unit": "ns",
            "extra": "gctime=0\nmemory=359952\nallocs=3424\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, V]",
            "value": 184804,
            "unit": "ns",
            "extra": "gctime=0\nmemory=55248\nallocs=1283\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, [S, H0]]",
            "value": 338523,
            "unit": "ns",
            "extra": "gctime=0\nmemory=95872\nallocs=2167\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_Dicke, S_j]",
            "value": 42640,
            "unit": "ns",
            "extra": "gctime=0\nmemory=54304\nallocs=545\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_JC, σ_j]",
            "value": 57882,
            "unit": "ns",
            "extra": "gctime=0\nmemory=61632\nallocs=674\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/double-sum spin-spin",
            "value": 16872,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8816\nallocs=130\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/indexed JC H",
            "value": 1174,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3488\nallocs=19\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/double Σ_ij(J_ij*S_i*S_j)",
            "value": 8433,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12984\nallocs=104\nparams={\"evals\":3,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/single Σ_i(σ_i*σ_j)",
            "value": 11015,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12288\nallocs=131\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=2",
            "value": 602.8977272727273,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1568\nallocs=8\nparams={\"evals\":176,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=3",
            "value": 856.4107142857143,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1904\nallocs=9\nparams={\"evals\":56,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=4",
            "value": 1132.7,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2352\nallocs=10\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=5",
            "value": 1509.6,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2848\nallocs=11\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Ground state/3-level rewrite",
            "value": 840.453125,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2816\nallocs=16\nparams={\"evals\":64,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Multi-mode/2-mode 6-op chain",
            "value": 925.2083333333334,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2080\nallocs=9\nparams={\"evals\":24,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H",
            "value": 1059.9,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3264\nallocs=17\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H²",
            "value": 2278.4444444444443,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5472\nallocs=27\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H",
            "value": 1087.2,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328\nallocs=17\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H²",
            "value": 4215.285714285715,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14432\nallocs=41\nparams={\"evals\":7,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H",
            "value": 1369.6,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=21\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H²",
            "value": 2220.222222222222,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4992\nallocs=29\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
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
          "id": "b78c7396f1f41cbf6e2de4174cd73965f4aac1a4",
          "message": "feat: proper operator substitution and time dependent parameter numerics (#198)\n\n* feat: proper operator substitution and time dependent parameter numerics\n\n* add tests\n\n* fix docs\n\n* more tests\n\n* add Coeff to API docs",
          "timestamp": "2026-06-28T19:52:34+02:00",
          "tree_id": "b7a21939e3399792063caec7f64080f9840df795",
          "url": "https://github.com/qojulia/SecondQuantizedAlgebra.jl/commit/b78c7396f1f41cbf6e2de4174cd73965f4aac1a4"
        },
        "date": 1782715066925,
        "tool": "julia",
        "benches": [
          {
            "name": "Accumulation/Many-mode H/foldl M=16",
            "value": 10921,
            "unit": "ns",
            "extra": "gctime=0\nmemory=59472\nallocs=213\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/foldl M=24",
            "value": 24356,
            "unit": "ns",
            "extra": "gctime=0\nmemory=131408\nallocs=436\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/foldl M=8",
            "value": 2789.6666666666665,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12432\nallocs=63\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=16",
            "value": 1683.1,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6816\nallocs=24\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=24",
            "value": 2269.8888888888887,
            "unit": "ns",
            "extra": "gctime=0\nmemory=7968\nallocs=32\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=8",
            "value": 700,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=13\nparams={\"evals\":145,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Same-site/foldl",
            "value": 6083.4,
            "unit": "ns",
            "extra": "gctime=0\nmemory=30912\nallocs=138\nparams={\"evals\":5,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Same-site/sum",
            "value": 2092.8888888888887,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4576\nallocs=29\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=1",
            "value": 1676.2,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4128\nallocs=39\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=2",
            "value": 17303,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26544\nallocs=276\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=3",
            "value": 64300,
            "unit": "ns",
            "extra": "gctime=0\nmemory=83728\nallocs=831\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=4",
            "value": 162434,
            "unit": "ns",
            "extra": "gctime=0\nmemory=183680\nallocs=1794\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=5",
            "value": 359613,
            "unit": "ns",
            "extra": "gctime=0\nmemory=360048\nallocs=3428\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, V]",
            "value": 201709,
            "unit": "ns",
            "extra": "gctime=0\nmemory=55248\nallocs=1283\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, [S, H0]]",
            "value": 368540,
            "unit": "ns",
            "extra": "gctime=0\nmemory=95872\nallocs=2167\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_Dicke, S_j]",
            "value": 43452,
            "unit": "ns",
            "extra": "gctime=0\nmemory=54304\nallocs=545\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_JC, σ_j]",
            "value": 63158,
            "unit": "ns",
            "extra": "gctime=0\nmemory=61632\nallocs=674\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/double-sum spin-spin",
            "value": 16079,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8752\nallocs=128\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/indexed JC H",
            "value": 1150.1,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3488\nallocs=19\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/double Σ_ij(J_ij*S_i*S_j)",
            "value": 8356,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12984\nallocs=104\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/single Σ_i(σ_i*σ_j)",
            "value": 11812,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12288\nallocs=131\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=2",
            "value": 637.2383720930233,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1568\nallocs=8\nparams={\"evals\":172,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=3",
            "value": 876.6538461538462,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1904\nallocs=9\nparams={\"evals\":52,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=4",
            "value": 1164.1,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2352\nallocs=10\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=5",
            "value": 1585,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2848\nallocs=11\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Ground state/3-level rewrite",
            "value": 796.1851851851852,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2816\nallocs=16\nparams={\"evals\":81,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Multi-mode/2-mode 6-op chain",
            "value": 971.8125,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2080\nallocs=9\nparams={\"evals\":16,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H",
            "value": 1062.9,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3264\nallocs=17\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H²",
            "value": 2370,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5472\nallocs=27\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H",
            "value": 1108.1,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328\nallocs=17\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H²",
            "value": 4333.857142857143,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14432\nallocs=41\nparams={\"evals\":7,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H",
            "value": 1359.6,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=21\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H²",
            "value": 2293.222222222222,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4992\nallocs=29\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
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
          "id": "77d5c36306e38f096191c54f9a0fa96f4a7f2b14",
          "message": "fix: ordering per space_index and fix printing bug (#206)",
          "timestamp": "2026-07-11T08:03:12+02:00",
          "tree_id": "1170ed1adf57663159d911aa738b7e37320fe2d8",
          "url": "https://github.com/qojulia/SecondQuantizedAlgebra.jl/commit/77d5c36306e38f096191c54f9a0fa96f4a7f2b14"
        },
        "date": 1783750268107,
        "tool": "julia",
        "benches": [
          {
            "name": "Accumulation/Many-mode H/foldl M=16",
            "value": 10600,
            "unit": "ns",
            "extra": "gctime=0\nmemory=59472\nallocs=213\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/foldl M=24",
            "value": 24066,
            "unit": "ns",
            "extra": "gctime=0\nmemory=131408\nallocs=436\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/foldl M=8",
            "value": 2763,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12432\nallocs=63\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=16",
            "value": 1724.3,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6816\nallocs=24\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=24",
            "value": 2276.4444444444443,
            "unit": "ns",
            "extra": "gctime=0\nmemory=7968\nallocs=32\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=8",
            "value": 712.55,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=13\nparams={\"evals\":140,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Same-site/foldl",
            "value": 5851,
            "unit": "ns",
            "extra": "gctime=0\nmemory=30912\nallocs=138\nparams={\"evals\":5,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Same-site/sum",
            "value": 2099.4444444444443,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4576\nallocs=29\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=1",
            "value": 1663.1,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4128\nallocs=39\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=2",
            "value": 17983,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26544\nallocs=276\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=3",
            "value": 67838,
            "unit": "ns",
            "extra": "gctime=0\nmemory=83728\nallocs=831\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=4",
            "value": 164799,
            "unit": "ns",
            "extra": "gctime=0\nmemory=183680\nallocs=1794\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=5",
            "value": 372339,
            "unit": "ns",
            "extra": "gctime=0\nmemory=360112\nallocs=3428\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, V]",
            "value": 203702,
            "unit": "ns",
            "extra": "gctime=0\nmemory=54736\nallocs=1263\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, [S, H0]]",
            "value": 373951,
            "unit": "ns",
            "extra": "gctime=0\nmemory=94816\nallocs=2134\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_Dicke, S_j]",
            "value": 43783,
            "unit": "ns",
            "extra": "gctime=0\nmemory=54304\nallocs=545\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_JC, σ_j]",
            "value": 62207,
            "unit": "ns",
            "extra": "gctime=0\nmemory=61376\nallocs=658\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/double-sum spin-spin",
            "value": 15950,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8688\nallocs=124\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/indexed JC H",
            "value": 1153.2,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3488\nallocs=19\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/double Σ_ij(J_ij*S_i*S_j)",
            "value": 8385.666666666666,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12904\nallocs=99\nparams={\"evals\":3,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/single Σ_i(σ_i*σ_j)",
            "value": 10620,
            "unit": "ns",
            "extra": "gctime=0\nmemory=11496\nallocs=118\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=2",
            "value": 602.4060606060606,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1568\nallocs=8\nparams={\"evals\":165,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=3",
            "value": 850.2272727272727,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1904\nallocs=9\nparams={\"evals\":22,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=4",
            "value": 1171.2,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2352\nallocs=10\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=5",
            "value": 1566.9,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2848\nallocs=11\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Ground state/3-level rewrite",
            "value": 805.2688172043011,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2816\nallocs=16\nparams={\"evals\":93,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Multi-mode/2-mode 6-op chain",
            "value": 953.8,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2080\nallocs=9\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H",
            "value": 1077,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3264\nallocs=17\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H²",
            "value": 2417.777777777778,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5472\nallocs=27\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H",
            "value": 1102.1,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328\nallocs=17\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H²",
            "value": 4467,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14432\nallocs=41\nparams={\"evals\":7,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H",
            "value": 1364.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=21\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H²",
            "value": 2299.777777777778,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4992\nallocs=29\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
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
          "id": "726f91eab66c5e1ae32660bca4cb4d12ec255f4b",
          "message": "fix: printing average between heterenous types (#207)",
          "timestamp": "2026-07-12T10:07:54+02:00",
          "tree_id": "a0e112fc6539e84abbad35f3a15a7d0e3ecffc00",
          "url": "https://github.com/qojulia/SecondQuantizedAlgebra.jl/commit/726f91eab66c5e1ae32660bca4cb4d12ec255f4b"
        },
        "date": 1783844128315,
        "tool": "julia",
        "benches": [
          {
            "name": "Accumulation/Many-mode H/foldl M=16",
            "value": 10890,
            "unit": "ns",
            "extra": "gctime=0\nmemory=59472\nallocs=213\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/foldl M=24",
            "value": 24326,
            "unit": "ns",
            "extra": "gctime=0\nmemory=131408\nallocs=436\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/foldl M=8",
            "value": 2778.5555555555557,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12432\nallocs=63\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=16",
            "value": 1701.2,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6816\nallocs=24\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=24",
            "value": 2314.3333333333335,
            "unit": "ns",
            "extra": "gctime=0\nmemory=7968\nallocs=32\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=8",
            "value": 713.9251700680272,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=13\nparams={\"evals\":147,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Same-site/foldl",
            "value": 6239.6,
            "unit": "ns",
            "extra": "gctime=0\nmemory=30912\nallocs=138\nparams={\"evals\":5,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Same-site/sum",
            "value": 2217.5555555555557,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4576\nallocs=29\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=1",
            "value": 1660.1,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4128\nallocs=39\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=2",
            "value": 17773,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26912\nallocs=280\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=3",
            "value": 65362,
            "unit": "ns",
            "extra": "gctime=0\nmemory=83360\nallocs=827\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=4",
            "value": 161393,
            "unit": "ns",
            "extra": "gctime=0\nmemory=181776\nallocs=1774\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=5",
            "value": 356668,
            "unit": "ns",
            "extra": "gctime=0\nmemory=355776\nallocs=3380\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, V]",
            "value": 207248,
            "unit": "ns",
            "extra": "gctime=0\nmemory=54736\nallocs=1263\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, [S, H0]]",
            "value": 375976,
            "unit": "ns",
            "extra": "gctime=0\nmemory=94816\nallocs=2134\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_Dicke, S_j]",
            "value": 43822,
            "unit": "ns",
            "extra": "gctime=0\nmemory=54304\nallocs=545\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_JC, σ_j]",
            "value": 64912,
            "unit": "ns",
            "extra": "gctime=0\nmemory=61376\nallocs=658\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/double-sum spin-spin",
            "value": 16661,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8688\nallocs=124\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/indexed JC H",
            "value": 1510.9,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3864\nallocs=25\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/double Σ_ij(J_ij*S_i*S_j)",
            "value": 8726.333333333334,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12904\nallocs=99\nparams={\"evals\":3,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/single Σ_i(σ_i*σ_j)",
            "value": 12082,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12224\nallocs=127\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=2",
            "value": 609.8241758241758,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1568\nallocs=8\nparams={\"evals\":182,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=3",
            "value": 886.0677966101695,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1904\nallocs=9\nparams={\"evals\":59,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=4",
            "value": 1178.2,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2352\nallocs=10\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=5",
            "value": 1573.9,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2848\nallocs=11\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Ground state/3-level rewrite",
            "value": 816.0430107526881,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2816\nallocs=16\nparams={\"evals\":93,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Multi-mode/2-mode 6-op chain",
            "value": 964.6111111111111,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2080\nallocs=9\nparams={\"evals\":18,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H",
            "value": 1076.1,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3264\nallocs=17\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H²",
            "value": 2456.777777777778,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5472\nallocs=27\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H",
            "value": 1119.1,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328\nallocs=17\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H²",
            "value": 4592.857142857143,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14432\nallocs=41\nparams={\"evals\":7,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H",
            "value": 1375.6,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=21\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H²",
            "value": 2289.777777777778,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4992\nallocs=29\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
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
          "id": "e8063ce2122658adb8336363f3b4f544a0a5f591",
          "message": "feat!: to_numeric returns a concrete (sparse) operator by default (#208)",
          "timestamp": "2026-07-12T16:44:45+02:00",
          "tree_id": "5a931454a29956e4dab34331e002aef77460a1e5",
          "url": "https://github.com/qojulia/SecondQuantizedAlgebra.jl/commit/e8063ce2122658adb8336363f3b4f544a0a5f591"
        },
        "date": 1783867936864,
        "tool": "julia",
        "benches": [
          {
            "name": "Accumulation/Many-mode H/foldl M=16",
            "value": 10760,
            "unit": "ns",
            "extra": "gctime=0\nmemory=59472\nallocs=213\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/foldl M=24",
            "value": 24176,
            "unit": "ns",
            "extra": "gctime=0\nmemory=131408\nallocs=436\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/foldl M=8",
            "value": 2725.1111111111113,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12432\nallocs=63\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=16",
            "value": 1696.2,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6816\nallocs=24\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=24",
            "value": 2308.777777777778,
            "unit": "ns",
            "extra": "gctime=0\nmemory=7968\nallocs=32\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=8",
            "value": 717.0364963503649,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=13\nparams={\"evals\":137,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Same-site/foldl",
            "value": 6061.4,
            "unit": "ns",
            "extra": "gctime=0\nmemory=30912\nallocs=138\nparams={\"evals\":5,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Same-site/sum",
            "value": 2168.4444444444443,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4576\nallocs=29\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=1",
            "value": 1726.3,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4128\nallocs=39\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=2",
            "value": 18545,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26912\nallocs=280\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=3",
            "value": 67886,
            "unit": "ns",
            "extra": "gctime=0\nmemory=83360\nallocs=827\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=4",
            "value": 166952,
            "unit": "ns",
            "extra": "gctime=0\nmemory=181472\nallocs=1770\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=5",
            "value": 365832,
            "unit": "ns",
            "extra": "gctime=0\nmemory=355200\nallocs=3372\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, V]",
            "value": 197969,
            "unit": "ns",
            "extra": "gctime=0\nmemory=54736\nallocs=1263\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, [S, H0]]",
            "value": 363508,
            "unit": "ns",
            "extra": "gctime=0\nmemory=94816\nallocs=2134\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_Dicke, S_j]",
            "value": 43722,
            "unit": "ns",
            "extra": "gctime=0\nmemory=54304\nallocs=545\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_JC, σ_j]",
            "value": 64821,
            "unit": "ns",
            "extra": "gctime=0\nmemory=61376\nallocs=658\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/double-sum spin-spin",
            "value": 15489,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8688\nallocs=124\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/indexed JC H",
            "value": 1562.9,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3864\nallocs=25\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/double Σ_ij(J_ij*S_i*S_j)",
            "value": 8532.333333333334,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12904\nallocs=99\nparams={\"evals\":3,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/single Σ_i(σ_i*σ_j)",
            "value": 12152,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12224\nallocs=127\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=2",
            "value": 598.5474860335196,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1568\nallocs=8\nparams={\"evals\":179,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=3",
            "value": 878.0967741935484,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1904\nallocs=9\nparams={\"evals\":62,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=4",
            "value": 1183.2,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2352\nallocs=10\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=5",
            "value": 1576,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2848\nallocs=11\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Ground state/3-level rewrite",
            "value": 826.25,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2816\nallocs=16\nparams={\"evals\":68,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Multi-mode/2-mode 6-op chain",
            "value": 975.1666666666666,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2080\nallocs=9\nparams={\"evals\":12,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H",
            "value": 1106.1,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3264\nallocs=17\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H²",
            "value": 2479,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5472\nallocs=27\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H",
            "value": 1156.2,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328\nallocs=17\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H²",
            "value": 4534.142857142857,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14432\nallocs=41\nparams={\"evals\":7,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H",
            "value": 1422.7,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=21\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H²",
            "value": 2360,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4992\nallocs=29\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
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
          "id": "0875b83e1be3f90a704fcd290e34e5f0f77e80db",
          "message": "Fix: diagonal coefficient leak for coefficient-only summation indices (#210)",
          "timestamp": "2026-07-15T13:47:53+02:00",
          "tree_id": "31e752af5c3776fb9898571f05cb13066a63143c",
          "url": "https://github.com/qojulia/SecondQuantizedAlgebra.jl/commit/0875b83e1be3f90a704fcd290e34e5f0f77e80db"
        },
        "date": 1784116553903,
        "tool": "julia",
        "benches": [
          {
            "name": "Accumulation/Many-mode H/foldl M=16",
            "value": 11241,
            "unit": "ns",
            "extra": "gctime=0\nmemory=59472\nallocs=213\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/foldl M=24",
            "value": 24977,
            "unit": "ns",
            "extra": "gctime=0\nmemory=131408\nallocs=436\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/foldl M=8",
            "value": 2949.25,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12432\nallocs=63\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=16",
            "value": 1773.4,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6816\nallocs=24\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=24",
            "value": 2372.222222222222,
            "unit": "ns",
            "extra": "gctime=0\nmemory=7968\nallocs=32\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=8",
            "value": 755.9913793103449,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=13\nparams={\"evals\":116,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Same-site/foldl",
            "value": 6893,
            "unit": "ns",
            "extra": "gctime=0\nmemory=30912\nallocs=138\nparams={\"evals\":4,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Same-site/sum",
            "value": 2311,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4576\nallocs=29\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=1",
            "value": 1716.2,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4128\nallocs=39\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=2",
            "value": 18204,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26912\nallocs=280\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=3",
            "value": 66155,
            "unit": "ns",
            "extra": "gctime=0\nmemory=83360\nallocs=827\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=4",
            "value": 163257,
            "unit": "ns",
            "extra": "gctime=0\nmemory=182144\nallocs=1778\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=5",
            "value": 361359,
            "unit": "ns",
            "extra": "gctime=0\nmemory=356144\nallocs=3384\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, V]",
            "value": 195597,
            "unit": "ns",
            "extra": "gctime=0\nmemory=54736\nallocs=1263\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, [S, H0]]",
            "value": 361450,
            "unit": "ns",
            "extra": "gctime=0\nmemory=94816\nallocs=2134\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_Dicke, S_j]",
            "value": 43662,
            "unit": "ns",
            "extra": "gctime=0\nmemory=54304\nallocs=545\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_JC, σ_j]",
            "value": 62958,
            "unit": "ns",
            "extra": "gctime=0\nmemory=61376\nallocs=658\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/double-sum spin-spin",
            "value": 15218,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8688\nallocs=124\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/indexed JC H",
            "value": 1518.9,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3864\nallocs=25\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/double Σ_ij(J_ij*S_i*S_j)",
            "value": 8282.333333333334,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12904\nallocs=99\nparams={\"evals\":3,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/single Σ_i(σ_i*σ_j)",
            "value": 10409,
            "unit": "ns",
            "extra": "gctime=0\nmemory=11496\nallocs=118\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=2",
            "value": 622.9640718562874,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1568\nallocs=8\nparams={\"evals\":167,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=3",
            "value": 905.4,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1904\nallocs=9\nparams={\"evals\":35,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=4",
            "value": 1203.2,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2352\nallocs=10\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=5",
            "value": 1581.9,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2848\nallocs=11\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Ground state/3-level rewrite",
            "value": 854.6190476190476,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2816\nallocs=16\nparams={\"evals\":63,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Multi-mode/2-mode 6-op chain",
            "value": 974.8,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2080\nallocs=9\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H",
            "value": 1113,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3264\nallocs=17\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H²",
            "value": 2509.1111111111113,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5472\nallocs=27\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H",
            "value": 1174.2,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328\nallocs=17\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H²",
            "value": 4686,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14432\nallocs=41\nparams={\"evals\":7,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H",
            "value": 1454.8,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=21\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H²",
            "value": 2395.6666666666665,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4992\nallocs=29\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
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
          "id": "e7e764890c6840fa576c167e2a69ef2563b8d839",
          "message": "Throw a clear error for String names (operators, Hilbert spaces, indices) (#212)",
          "timestamp": "2026-07-16T14:50:18+02:00",
          "tree_id": "13256f302766ee88d2a5409085d47d3287e0c4ce",
          "url": "https://github.com/qojulia/SecondQuantizedAlgebra.jl/commit/e7e764890c6840fa576c167e2a69ef2563b8d839"
        },
        "date": 1784206663537,
        "tool": "julia",
        "benches": [
          {
            "name": "Accumulation/Many-mode H/foldl M=16",
            "value": 10980,
            "unit": "ns",
            "extra": "gctime=0\nmemory=59472\nallocs=213\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/foldl M=24",
            "value": 24876,
            "unit": "ns",
            "extra": "gctime=0\nmemory=131408\nallocs=436\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/foldl M=8",
            "value": 2800.777777777778,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12432\nallocs=63\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=16",
            "value": 1697.1,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6816\nallocs=24\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=24",
            "value": 2305.4444444444443,
            "unit": "ns",
            "extra": "gctime=0\nmemory=7968\nallocs=32\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=8",
            "value": 729.7297297297297,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=13\nparams={\"evals\":148,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Same-site/foldl",
            "value": 6275.6,
            "unit": "ns",
            "extra": "gctime=0\nmemory=30912\nallocs=138\nparams={\"evals\":5,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Same-site/sum",
            "value": 2295.3333333333335,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4576\nallocs=29\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=1",
            "value": 1705.1,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4128\nallocs=39\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=2",
            "value": 17583,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26544\nallocs=276\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=3",
            "value": 66503,
            "unit": "ns",
            "extra": "gctime=0\nmemory=83728\nallocs=831\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=4",
            "value": 168372,
            "unit": "ns",
            "extra": "gctime=0\nmemory=187072\nallocs=1797\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=5",
            "value": 372361,
            "unit": "ns",
            "extra": "gctime=0\nmemory=363440\nallocs=3431\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, V]",
            "value": 198388,
            "unit": "ns",
            "extra": "gctime=0\nmemory=54736\nallocs=1263\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, [S, H0]]",
            "value": 367321,
            "unit": "ns",
            "extra": "gctime=0\nmemory=94816\nallocs=2134\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_Dicke, S_j]",
            "value": 44253,
            "unit": "ns",
            "extra": "gctime=0\nmemory=54304\nallocs=545\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_JC, σ_j]",
            "value": 63227,
            "unit": "ns",
            "extra": "gctime=0\nmemory=61376\nallocs=658\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/double-sum spin-spin",
            "value": 16551,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8688\nallocs=124\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/indexed JC H",
            "value": 1169.1,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3488\nallocs=19\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/double Σ_ij(J_ij*S_i*S_j)",
            "value": 8766,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12904\nallocs=99\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/single Σ_i(σ_i*σ_j)",
            "value": 10669,
            "unit": "ns",
            "extra": "gctime=0\nmemory=11496\nallocs=118\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=2",
            "value": 611.7403314917127,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1568\nallocs=8\nparams={\"evals\":181,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=3",
            "value": 892.1754385964912,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1904\nallocs=9\nparams={\"evals\":57,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=4",
            "value": 1183.2,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2352\nallocs=10\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=5",
            "value": 1565.9,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2848\nallocs=11\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Ground state/3-level rewrite",
            "value": 809.5487804878048,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2816\nallocs=16\nparams={\"evals\":82,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Multi-mode/2-mode 6-op chain",
            "value": 964.421052631579,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2080\nallocs=9\nparams={\"evals\":19,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H",
            "value": 1050,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3264\nallocs=17\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H²",
            "value": 2412.3333333333335,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5472\nallocs=27\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H",
            "value": 1103,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328\nallocs=17\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H²",
            "value": 4495.428571428572,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14432\nallocs=41\nparams={\"evals\":7,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H",
            "value": 1352.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=21\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H²",
            "value": 2276.4444444444443,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4992\nallocs=29\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
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
          "id": "6e26e93c6a2a6a88e208075833873f9b75a4cf26",
          "message": "fix: render index slot suffixes as a comma subscript (#213)",
          "timestamp": "2026-07-16T18:05:02+02:00",
          "tree_id": "efdaece266c89c9e51705ee6ca3326c5cd06a5a0",
          "url": "https://github.com/qojulia/SecondQuantizedAlgebra.jl/commit/6e26e93c6a2a6a88e208075833873f9b75a4cf26"
        },
        "date": 1784218268293,
        "tool": "julia",
        "benches": [
          {
            "name": "Accumulation/Many-mode H/foldl M=16",
            "value": 12359,
            "unit": "ns",
            "extra": "gctime=0\nmemory=59472\nallocs=213\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/foldl M=24",
            "value": 27461,
            "unit": "ns",
            "extra": "gctime=0\nmemory=131408\nallocs=436\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/foldl M=8",
            "value": 2963.125,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12432\nallocs=63\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=16",
            "value": 1908.4444444444443,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6816\nallocs=24\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=24",
            "value": 2582.6666666666665,
            "unit": "ns",
            "extra": "gctime=0\nmemory=7968\nallocs=32\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=8",
            "value": 768.7222222222222,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=13\nparams={\"evals\":108,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Same-site/foldl",
            "value": 6477.25,
            "unit": "ns",
            "extra": "gctime=0\nmemory=30912\nallocs=138\nparams={\"evals\":4,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Same-site/sum",
            "value": 2283.4444444444443,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4576\nallocs=29\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=1",
            "value": 1811.7,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4128\nallocs=39\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=2",
            "value": 17836,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26912\nallocs=280\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=3",
            "value": 65877,
            "unit": "ns",
            "extra": "gctime=0\nmemory=84464\nallocs=839\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=4",
            "value": 160948,
            "unit": "ns",
            "extra": "gctime=0\nmemory=184784\nallocs=1806\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=5",
            "value": 361504,
            "unit": "ns",
            "extra": "gctime=0\nmemory=361648\nallocs=3440\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, V]",
            "value": 178123,
            "unit": "ns",
            "extra": "gctime=0\nmemory=54736\nallocs=1263\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, [S, H0]]",
            "value": 324459,
            "unit": "ns",
            "extra": "gctime=0\nmemory=94816\nallocs=2134\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_Dicke, S_j]",
            "value": 44395,
            "unit": "ns",
            "extra": "gctime=0\nmemory=54304\nallocs=545\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_JC, σ_j]",
            "value": 58025,
            "unit": "ns",
            "extra": "gctime=0\nmemory=61376\nallocs=658\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/double-sum spin-spin",
            "value": 14582,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8688\nallocs=124\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/indexed JC H",
            "value": 1278.9,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3488\nallocs=19\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/double Σ_ij(J_ij*S_i*S_j)",
            "value": 7935,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12904\nallocs=99\nparams={\"evals\":3,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/single Σ_i(σ_i*σ_j)",
            "value": 9704,
            "unit": "ns",
            "extra": "gctime=0\nmemory=11496\nallocs=118\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=2",
            "value": 668.7375,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1568\nallocs=8\nparams={\"evals\":160,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=3",
            "value": 938.4705882352941,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1904\nallocs=9\nparams={\"evals\":17,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=4",
            "value": 1281.9,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2352\nallocs=10\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=5",
            "value": 1711.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2848\nallocs=11\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Ground state/3-level rewrite",
            "value": 862.7592592592592,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2816\nallocs=16\nparams={\"evals\":54,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Multi-mode/2-mode 6-op chain",
            "value": 1061.6,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2080\nallocs=9\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H",
            "value": 1197.7,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3264\nallocs=17\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H²",
            "value": 2697.75,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5472\nallocs=27\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H",
            "value": 1247.8,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328\nallocs=17\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H²",
            "value": 5147.666666666667,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14432\nallocs=41\nparams={\"evals\":6,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H",
            "value": 1533.2,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=21\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H²",
            "value": 2581.5555555555557,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4992\nallocs=29\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "chr.hotter@gmail.com",
            "name": "Christoph Hotter",
            "username": "ChristophHotter"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "764b0a2dd874fb110df6dbf700650e01ba653a91",
          "message": "CollectiveTransiton implementation (#211)\n\n* CollectiveTransiton implementation\n\n* format and typos\n\n* coll tansition test alloc\n\n* enum CollTransition 7\n\n* changelog\n\n* docs error fix\n\n* superradiant decay example\n\n* superraiance example\n\n* rename to tavis cummings\n\n* _commute_ops_at reduce cost\n\n* format",
          "timestamp": "2026-07-16T13:29:40-04:00",
          "tree_id": "b3389a71aabeef1d324e9201e1e70d09b5064b45",
          "url": "https://github.com/qojulia/SecondQuantizedAlgebra.jl/commit/764b0a2dd874fb110df6dbf700650e01ba653a91"
        },
        "date": 1784223480911,
        "tool": "julia",
        "benches": [
          {
            "name": "Accumulation/Many-mode H/foldl M=16",
            "value": 6999,
            "unit": "ns",
            "extra": "gctime=0\nmemory=59472\nallocs=213\nparams={\"evals\":3,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/foldl M=24",
            "value": 15130,
            "unit": "ns",
            "extra": "gctime=0\nmemory=131408\nallocs=436\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/foldl M=8",
            "value": 1979.4444444444443,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12432\nallocs=63\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=16",
            "value": 968.3,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6816\nallocs=24\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=24",
            "value": 1397.6,
            "unit": "ns",
            "extra": "gctime=0\nmemory=7968\nallocs=32\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=8",
            "value": 424.2525773195876,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=13\nparams={\"evals\":194,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Same-site/foldl",
            "value": 3845.5714285714284,
            "unit": "ns",
            "extra": "gctime=0\nmemory=30912\nallocs=138\nparams={\"evals\":7,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Same-site/sum",
            "value": 1347.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4576\nallocs=29\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=1",
            "value": 982.7,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4128\nallocs=39\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=2",
            "value": 10508,
            "unit": "ns",
            "extra": "gctime=0\nmemory=27104\nallocs=283\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=3",
            "value": 39295,
            "unit": "ns",
            "extra": "gctime=0\nmemory=84192\nallocs=840\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=4",
            "value": 97546,
            "unit": "ns",
            "extra": "gctime=0\nmemory=183456\nallocs=1802\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=5",
            "value": 214652,
            "unit": "ns",
            "extra": "gctime=0\nmemory=359008\nallocs=3436\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, V]",
            "value": 105690,
            "unit": "ns",
            "extra": "gctime=0\nmemory=55920\nallocs=1299\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, [S, H0]]",
            "value": 197319,
            "unit": "ns",
            "extra": "gctime=0\nmemory=97184\nallocs=2206\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_Dicke, S_j]",
            "value": 25150,
            "unit": "ns",
            "extra": "gctime=0\nmemory=54432\nallocs=547\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_JC, σ_j]",
            "value": 33352,
            "unit": "ns",
            "extra": "gctime=0\nmemory=61376\nallocs=658\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/double-sum spin-spin",
            "value": 8483,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8688\nallocs=124\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/indexed JC H",
            "value": 880.875,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3864\nallocs=25\nparams={\"evals\":24,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/double Σ_ij(J_ij*S_i*S_j)",
            "value": 4960.333333333333,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12904\nallocs=99\nparams={\"evals\":6,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/single Σ_i(σ_i*σ_j)",
            "value": 5981.2,
            "unit": "ns",
            "extra": "gctime=0\nmemory=11496\nallocs=118\nparams={\"evals\":5,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=2",
            "value": 368.00975609756097,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1568\nallocs=8\nparams={\"evals\":205,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=3",
            "value": 534.8,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1904\nallocs=9\nparams={\"evals\":190,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=4",
            "value": 748.4778761061947,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2352\nallocs=10\nparams={\"evals\":113,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=5",
            "value": 946.0714285714286,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2848\nallocs=11\nparams={\"evals\":14,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Ground state/3-level rewrite",
            "value": 532.9358288770053,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2816\nallocs=16\nparams={\"evals\":187,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Multi-mode/2-mode 6-op chain",
            "value": 605.6279069767442,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2080\nallocs=9\nparams={\"evals\":172,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H",
            "value": 674.5772357723578,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3264\nallocs=17\nparams={\"evals\":123,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H²",
            "value": 1428.4,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5472\nallocs=27\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H",
            "value": 684.5106382978723,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328\nallocs=17\nparams={\"evals\":141,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H²",
            "value": 2693.125,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14432\nallocs=41\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H",
            "value": 847.9056603773585,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=21\nparams={\"evals\":53,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H²",
            "value": 1398.4,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4992\nallocs=29\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
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
          "id": "7fd819bfd542ec73ffc01579ddffe4be1830b6aa",
          "message": "fix: revert version to 0.9.4 in Project.toml and update Changelog (#215)",
          "timestamp": "2026-07-16T20:13:03+02:00",
          "tree_id": "65142d85424ebcecf047a033f06f6e8cc65f21f3",
          "url": "https://github.com/qojulia/SecondQuantizedAlgebra.jl/commit/7fd819bfd542ec73ffc01579ddffe4be1830b6aa"
        },
        "date": 1784226193865,
        "tool": "julia",
        "benches": [
          {
            "name": "Accumulation/Many-mode H/foldl M=16",
            "value": 10790,
            "unit": "ns",
            "extra": "gctime=0\nmemory=59472\nallocs=213\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/foldl M=24",
            "value": 25428,
            "unit": "ns",
            "extra": "gctime=0\nmemory=138000\nallocs=439\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/foldl M=8",
            "value": 2763,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12432\nallocs=63\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=16",
            "value": 1705.2,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6816\nallocs=24\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=24",
            "value": 2308.777777777778,
            "unit": "ns",
            "extra": "gctime=0\nmemory=7968\nallocs=32\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=8",
            "value": 731.5179856115108,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=13\nparams={\"evals\":139,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Same-site/foldl",
            "value": 6275.6,
            "unit": "ns",
            "extra": "gctime=0\nmemory=30912\nallocs=138\nparams={\"evals\":5,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Same-site/sum",
            "value": 2163,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4576\nallocs=29\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=1",
            "value": 1660.2,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4128\nallocs=39\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=2",
            "value": 18395,
            "unit": "ns",
            "extra": "gctime=0\nmemory=27104\nallocs=283\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=3",
            "value": 68008,
            "unit": "ns",
            "extra": "gctime=0\nmemory=84928\nallocs=848\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=4",
            "value": 173205,
            "unit": "ns",
            "extra": "gctime=0\nmemory=186160\nallocs=1830\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=5",
            "value": 381763,
            "unit": "ns",
            "extra": "gctime=0\nmemory=364736\nallocs=3492\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, V]",
            "value": 201398,
            "unit": "ns",
            "extra": "gctime=0\nmemory=55920\nallocs=1299\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, [S, H0]]",
            "value": 376591,
            "unit": "ns",
            "extra": "gctime=0\nmemory=97184\nallocs=2206\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_Dicke, S_j]",
            "value": 43612,
            "unit": "ns",
            "extra": "gctime=0\nmemory=54432\nallocs=547\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_JC, σ_j]",
            "value": 62036,
            "unit": "ns",
            "extra": "gctime=0\nmemory=61376\nallocs=658\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/double-sum spin-spin",
            "value": 15839,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8688\nallocs=124\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/indexed JC H",
            "value": 1164.2,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3488\nallocs=19\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/double Σ_ij(J_ij*S_i*S_j)",
            "value": 8435.666666666666,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12904\nallocs=99\nparams={\"evals\":3,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/single Σ_i(σ_i*σ_j)",
            "value": 11361,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12224\nallocs=127\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=2",
            "value": 605.1657458563536,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1568\nallocs=8\nparams={\"evals\":181,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=3",
            "value": 884.3125,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1904\nallocs=9\nparams={\"evals\":64,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=4",
            "value": 1175.2,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2352\nallocs=10\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=5",
            "value": 1565.9,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2848\nallocs=11\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Ground state/3-level rewrite",
            "value": 809.4270833333334,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2816\nallocs=16\nparams={\"evals\":96,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Multi-mode/2-mode 6-op chain",
            "value": 983.9285714285714,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2080\nallocs=9\nparams={\"evals\":14,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H",
            "value": 1063.1,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3264\nallocs=17\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H²",
            "value": 2373.3333333333335,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5472\nallocs=27\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H",
            "value": 1103.1,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328\nallocs=17\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H²",
            "value": 4459.714285714285,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14432\nallocs=41\nparams={\"evals\":7,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H",
            "value": 1355.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=21\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H²",
            "value": 2284.3333333333335,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4992\nallocs=29\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
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
          "id": "37a8cf03598e7963ddbc4b530cc91706ca16b607",
          "message": "feat: public API for the indexed-sum node and Op/Index rebinding (#214)",
          "timestamp": "2026-07-16T20:13:29+02:00",
          "tree_id": "1f20b8595bf08f3446036fc529e2b31549bcdca7",
          "url": "https://github.com/qojulia/SecondQuantizedAlgebra.jl/commit/37a8cf03598e7963ddbc4b530cc91706ca16b607"
        },
        "date": 1784226243337,
        "tool": "julia",
        "benches": [
          {
            "name": "Accumulation/Many-mode H/foldl M=16",
            "value": 11658,
            "unit": "ns",
            "extra": "gctime=0\nmemory=59472\nallocs=213\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/foldl M=24",
            "value": 26270,
            "unit": "ns",
            "extra": "gctime=0\nmemory=131408\nallocs=436\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/foldl M=8",
            "value": 2937,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12432\nallocs=63\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=16",
            "value": 1878.4444444444443,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6816\nallocs=24\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=24",
            "value": 2391.3333333333335,
            "unit": "ns",
            "extra": "gctime=0\nmemory=7968\nallocs=32\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=8",
            "value": 757.0175438596491,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=13\nparams={\"evals\":114,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Same-site/foldl",
            "value": 6570,
            "unit": "ns",
            "extra": "gctime=0\nmemory=30912\nallocs=138\nparams={\"evals\":4,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Same-site/sum",
            "value": 2242.222222222222,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4576\nallocs=29\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=1",
            "value": 1801.7,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4128\nallocs=39\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=2",
            "value": 17496,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26736\nallocs=279\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=3",
            "value": 64547,
            "unit": "ns",
            "extra": "gctime=0\nmemory=84928\nallocs=848\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=4",
            "value": 161684,
            "unit": "ns",
            "extra": "gctime=0\nmemory=186864\nallocs=1838\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=5",
            "value": 354355,
            "unit": "ns",
            "extra": "gctime=0\nmemory=366080\nallocs=3512\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, V]",
            "value": 187483,
            "unit": "ns",
            "extra": "gctime=0\nmemory=55920\nallocs=1299\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, [S, H0]]",
            "value": 344670,
            "unit": "ns",
            "extra": "gctime=0\nmemory=97184\nallocs=2206\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_Dicke, S_j]",
            "value": 44597,
            "unit": "ns",
            "extra": "gctime=0\nmemory=54432\nallocs=547\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_JC, σ_j]",
            "value": 58488,
            "unit": "ns",
            "extra": "gctime=0\nmemory=61376\nallocs=658\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/double-sum spin-spin",
            "value": 14712,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8688\nallocs=124\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/indexed JC H",
            "value": 1265.9,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3488\nallocs=19\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/double Σ_ij(J_ij*S_i*S_j)",
            "value": 8039.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12904\nallocs=99\nparams={\"evals\":4,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/single Σ_i(σ_i*σ_j)",
            "value": 10556,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12224\nallocs=127\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=2",
            "value": 660.2405063291139,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1568\nallocs=8\nparams={\"evals\":158,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=3",
            "value": 912.9230769230769,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1904\nallocs=9\nparams={\"evals\":13,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=4",
            "value": 1269,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2352\nallocs=10\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=5",
            "value": 1679.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2848\nallocs=11\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Ground state/3-level rewrite",
            "value": 889.1355932203389,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2816\nallocs=16\nparams={\"evals\":59,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Multi-mode/2-mode 6-op chain",
            "value": 1038.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2080\nallocs=9\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H",
            "value": 1193.8,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3264\nallocs=17\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H²",
            "value": 2688.4444444444443,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5472\nallocs=27\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H",
            "value": 1228.9,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328\nallocs=17\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H²",
            "value": 5015.833333333333,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14432\nallocs=41\nparams={\"evals\":6,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H",
            "value": 1508.3,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=21\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H²",
            "value": 2550.5555555555557,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4992\nallocs=29\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
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
          "id": "3b4b786352fb429c9e086aa0793eafc57d01184f",
          "message": "feat: update symtype handling for Hermitian averages (#217)",
          "timestamp": "2026-07-17T14:41:43+02:00",
          "tree_id": "47cc5c9e22f898c57ae5e71f0e766a4185817d3a",
          "url": "https://github.com/qojulia/SecondQuantizedAlgebra.jl/commit/3b4b786352fb429c9e086aa0793eafc57d01184f"
        },
        "date": 1784292590983,
        "tool": "julia",
        "benches": [
          {
            "name": "Accumulation/Many-mode H/foldl M=16",
            "value": 10639,
            "unit": "ns",
            "extra": "gctime=0\nmemory=59472\nallocs=213\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/foldl M=24",
            "value": 24345,
            "unit": "ns",
            "extra": "gctime=0\nmemory=131408\nallocs=436\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/foldl M=8",
            "value": 2727.3333333333335,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12432\nallocs=63\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=16",
            "value": 1704.1,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6816\nallocs=24\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=24",
            "value": 2265.222222222222,
            "unit": "ns",
            "extra": "gctime=0\nmemory=7968\nallocs=32\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=8",
            "value": 708.1363636363636,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=13\nparams={\"evals\":132,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Same-site/foldl",
            "value": 6313.6,
            "unit": "ns",
            "extra": "gctime=0\nmemory=30912\nallocs=138\nparams={\"evals\":5,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Same-site/sum",
            "value": 2259.777777777778,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4576\nallocs=29\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=1",
            "value": 1683.2,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4128\nallocs=39\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=2",
            "value": 18324,
            "unit": "ns",
            "extra": "gctime=0\nmemory=27104\nallocs=283\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=3",
            "value": 67676,
            "unit": "ns",
            "extra": "gctime=0\nmemory=85296\nallocs=852\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=4",
            "value": 170367,
            "unit": "ns",
            "extra": "gctime=0\nmemory=187264\nallocs=1842\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=5",
            "value": 376862,
            "unit": "ns",
            "extra": "gctime=0\nmemory=366768\nallocs=3516\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, V]",
            "value": 205001,
            "unit": "ns",
            "extra": "gctime=0\nmemory=55920\nallocs=1299\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, [S, H0]]",
            "value": 383264,
            "unit": "ns",
            "extra": "gctime=0\nmemory=97184\nallocs=2206\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_Dicke, S_j]",
            "value": 42960,
            "unit": "ns",
            "extra": "gctime=0\nmemory=54432\nallocs=547\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_JC, σ_j]",
            "value": 62306,
            "unit": "ns",
            "extra": "gctime=0\nmemory=61376\nallocs=658\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/double-sum spin-spin",
            "value": 16832,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8688\nallocs=124\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/indexed JC H",
            "value": 1177.2,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3488\nallocs=19\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/double Σ_ij(J_ij*S_i*S_j)",
            "value": 8696.333333333334,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12904\nallocs=99\nparams={\"evals\":3,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/single Σ_i(σ_i*σ_j)",
            "value": 11411,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12224\nallocs=127\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=2",
            "value": 613.4715909090909,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1568\nallocs=8\nparams={\"evals\":176,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=3",
            "value": 882.5813953488372,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1904\nallocs=9\nparams={\"evals\":43,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=4",
            "value": 1172.1,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2352\nallocs=10\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=5",
            "value": 1570.9,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2848\nallocs=11\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Ground state/3-level rewrite",
            "value": 823.7402597402597,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2816\nallocs=16\nparams={\"evals\":77,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Multi-mode/2-mode 6-op chain",
            "value": 962.7,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2080\nallocs=9\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H",
            "value": 1091,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3264\nallocs=17\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H²",
            "value": 2466.8888888888887,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5472\nallocs=27\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H",
            "value": 1139.1,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328\nallocs=17\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H²",
            "value": 4554.285714285715,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14432\nallocs=41\nparams={\"evals\":7,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H",
            "value": 1384.6,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=21\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H²",
            "value": 2375.5555555555557,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4992\nallocs=29\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
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
          "id": "5a37827bed36cc25d4de790651e8bda049ec3b85",
          "message": "feat: backend agnostic numerical conversion (#200)",
          "timestamp": "2026-07-20T19:59:29+02:00",
          "tree_id": "fa52cc896b108cb001755a90b955a9afb90da4a1",
          "url": "https://github.com/qojulia/SecondQuantizedAlgebra.jl/commit/5a37827bed36cc25d4de790651e8bda049ec3b85"
        },
        "date": 1784570817441,
        "tool": "julia",
        "benches": [
          {
            "name": "Accumulation/Many-mode H/foldl M=16",
            "value": 10870,
            "unit": "ns",
            "extra": "gctime=0\nmemory=59472\nallocs=213\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/foldl M=24",
            "value": 24726,
            "unit": "ns",
            "extra": "gctime=0\nmemory=131408\nallocs=436\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/foldl M=8",
            "value": 2800.8888888888887,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12432\nallocs=63\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=16",
            "value": 1683.1,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6816\nallocs=24\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=24",
            "value": 2298.777777777778,
            "unit": "ns",
            "extra": "gctime=0\nmemory=7968\nallocs=32\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=8",
            "value": 713.2805755395683,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=13\nparams={\"evals\":139,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Same-site/foldl",
            "value": 6229.6,
            "unit": "ns",
            "extra": "gctime=0\nmemory=30912\nallocs=138\nparams={\"evals\":5,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Same-site/sum",
            "value": 2101.777777777778,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4576\nallocs=29\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=1",
            "value": 1678.1,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4128\nallocs=39\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=2",
            "value": 18254,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26736\nallocs=279\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=3",
            "value": 68689,
            "unit": "ns",
            "extra": "gctime=0\nmemory=84560\nallocs=844\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=4",
            "value": 170841,
            "unit": "ns",
            "extra": "gctime=0\nmemory=185728\nallocs=1826\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=5",
            "value": 377859,
            "unit": "ns",
            "extra": "gctime=0\nmemory=364112\nallocs=3492\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, V]",
            "value": 201047,
            "unit": "ns",
            "extra": "gctime=0\nmemory=55920\nallocs=1299\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, [S, H0]]",
            "value": 375935,
            "unit": "ns",
            "extra": "gctime=0\nmemory=97184\nallocs=2206\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_Dicke, S_j]",
            "value": 43342,
            "unit": "ns",
            "extra": "gctime=0\nmemory=54432\nallocs=547\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_JC, σ_j]",
            "value": 61826,
            "unit": "ns",
            "extra": "gctime=0\nmemory=61376\nallocs=658\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/double-sum spin-spin",
            "value": 16030,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8688\nallocs=124\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/indexed JC H",
            "value": 1135.1,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3488\nallocs=19\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/double Σ_ij(J_ij*S_i*S_j)",
            "value": 8596.333333333334,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12904\nallocs=99\nparams={\"evals\":3,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/single Σ_i(σ_i*σ_j)",
            "value": 11642,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12224\nallocs=127\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=2",
            "value": 602.224043715847,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1568\nallocs=8\nparams={\"evals\":183,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=3",
            "value": 869,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1904\nallocs=9\nparams={\"evals\":76,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=4",
            "value": 1154.2,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2352\nallocs=10\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=5",
            "value": 1539.9,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2848\nallocs=11\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Ground state/3-level rewrite",
            "value": 818.7183098591549,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2816\nallocs=16\nparams={\"evals\":71,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Multi-mode/2-mode 6-op chain",
            "value": 957.578947368421,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2080\nallocs=9\nparams={\"evals\":19,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H",
            "value": 1057,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3264\nallocs=17\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H²",
            "value": 2381.1111111111113,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5472\nallocs=27\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H",
            "value": 1083.1,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328\nallocs=17\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H²",
            "value": 4564.285714285715,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14432\nallocs=41\nparams={\"evals\":7,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H",
            "value": 1350.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=21\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H²",
            "value": 2286.5555555555557,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4992\nallocs=29\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
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
          "id": "d338716ee3e4ff1e4277b9ebdd404333c24a4f4d",
          "message": "benchmark against QuantumAlgebra (#163)",
          "timestamp": "2026-07-21T21:10:25+02:00",
          "tree_id": "c2bcb00fb757244ed04052b8bfff83455f0f18dd",
          "url": "https://github.com/qojulia/SecondQuantizedAlgebra.jl/commit/d338716ee3e4ff1e4277b9ebdd404333c24a4f4d"
        },
        "date": 1784661244009,
        "tool": "julia",
        "benches": [
          {
            "name": "Accumulation/Many-mode H/foldl M=16",
            "value": 10760,
            "unit": "ns",
            "extra": "gctime=0\nmemory=59472\nallocs=213\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/foldl M=24",
            "value": 24406,
            "unit": "ns",
            "extra": "gctime=0\nmemory=131408\nallocs=436\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/foldl M=8",
            "value": 2735.1111111111113,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12432\nallocs=63\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=16",
            "value": 1740.2,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6816\nallocs=24\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=24",
            "value": 2311,
            "unit": "ns",
            "extra": "gctime=0\nmemory=7968\nallocs=32\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=8",
            "value": 693.7214285714285,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=13\nparams={\"evals\":140,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Same-site/foldl",
            "value": 6289.8,
            "unit": "ns",
            "extra": "gctime=0\nmemory=30912\nallocs=138\nparams={\"evals\":5,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Same-site/sum",
            "value": 2144,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4576\nallocs=29\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=1",
            "value": 1712.2,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4128\nallocs=39\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=2",
            "value": 17934,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26736\nallocs=279\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=3",
            "value": 68719,
            "unit": "ns",
            "extra": "gctime=0\nmemory=84928\nallocs=848\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=4",
            "value": 170578,
            "unit": "ns",
            "extra": "gctime=0\nmemory=186432\nallocs=1834\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=5",
            "value": 384539,
            "unit": "ns",
            "extra": "gctime=0\nmemory=366112\nallocs=3508\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, V]",
            "value": 203120,
            "unit": "ns",
            "extra": "gctime=0\nmemory=55920\nallocs=1299\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, [S, H0]]",
            "value": 376524,
            "unit": "ns",
            "extra": "gctime=0\nmemory=97184\nallocs=2206\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_Dicke, S_j]",
            "value": 44453,
            "unit": "ns",
            "extra": "gctime=0\nmemory=54432\nallocs=547\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_JC, σ_j]",
            "value": 63929,
            "unit": "ns",
            "extra": "gctime=0\nmemory=61376\nallocs=658\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/double-sum spin-spin",
            "value": 15970,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8688\nallocs=124\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/indexed JC H",
            "value": 1159.2,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3488\nallocs=19\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/double Σ_ij(J_ij*S_i*S_j)",
            "value": 8422.666666666666,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12904\nallocs=99\nparams={\"evals\":3,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/single Σ_i(σ_i*σ_j)",
            "value": 10830,
            "unit": "ns",
            "extra": "gctime=0\nmemory=11496\nallocs=118\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=2",
            "value": 607.4269662921348,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1568\nallocs=8\nparams={\"evals\":178,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=3",
            "value": 884.6964285714286,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1904\nallocs=9\nparams={\"evals\":56,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=4",
            "value": 1189.2,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2352\nallocs=10\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=5",
            "value": 1571.9,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2848\nallocs=11\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Ground state/3-level rewrite",
            "value": 819.506329113924,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2816\nallocs=16\nparams={\"evals\":79,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Multi-mode/2-mode 6-op chain",
            "value": 970.4666666666667,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2080\nallocs=9\nparams={\"evals\":15,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H",
            "value": 1163.1,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3264\nallocs=17\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H²",
            "value": 2405.5555555555557,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5472\nallocs=27\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H",
            "value": 1112,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328\nallocs=17\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H²",
            "value": 4509.857142857143,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14432\nallocs=41\nparams={\"evals\":7,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H",
            "value": 1366.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=21\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H²",
            "value": 2318.777777777778,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4992\nallocs=29\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
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
          "id": "8496348c264086c3d8e9078f0338afd06e5dce53",
          "message": "Fold elementary functions of a literal zero (fixes e^0 in printed equations) (#221)",
          "timestamp": "2026-07-22T14:23:43+02:00",
          "tree_id": "c8453fec46ccfddf21ebd009986e3265999e089f",
          "url": "https://github.com/qojulia/SecondQuantizedAlgebra.jl/commit/8496348c264086c3d8e9078f0338afd06e5dce53"
        },
        "date": 1784723434234,
        "tool": "julia",
        "benches": [
          {
            "name": "Accumulation/Many-mode H/foldl M=16",
            "value": 11848,
            "unit": "ns",
            "extra": "gctime=0\nmemory=59472\nallocs=213\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/foldl M=24",
            "value": 26701,
            "unit": "ns",
            "extra": "gctime=0\nmemory=131408\nallocs=436\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/foldl M=8",
            "value": 2903.125,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12432\nallocs=63\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=16",
            "value": 1864.8,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6816\nallocs=24\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=24",
            "value": 2452.5555555555557,
            "unit": "ns",
            "extra": "gctime=0\nmemory=7968\nallocs=32\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=8",
            "value": 759.7540983606557,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=13\nparams={\"evals\":122,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Same-site/foldl",
            "value": 6367.6,
            "unit": "ns",
            "extra": "gctime=0\nmemory=30912\nallocs=138\nparams={\"evals\":5,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Same-site/sum",
            "value": 2222.222222222222,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4576\nallocs=29\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=1",
            "value": 1776.7,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4128\nallocs=39\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=2",
            "value": 17417,
            "unit": "ns",
            "extra": "gctime=0\nmemory=27104\nallocs=283\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=3",
            "value": 63987,
            "unit": "ns",
            "extra": "gctime=0\nmemory=84928\nallocs=848\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=4",
            "value": 159892,
            "unit": "ns",
            "extra": "gctime=0\nmemory=186160\nallocs=1830\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=5",
            "value": 355547,
            "unit": "ns",
            "extra": "gctime=0\nmemory=364816\nallocs=3492\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, V]",
            "value": 185771,
            "unit": "ns",
            "extra": "gctime=0\nmemory=55920\nallocs=1299\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, [S, H0]]",
            "value": 342798,
            "unit": "ns",
            "extra": "gctime=0\nmemory=97184\nallocs=2206\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_Dicke, S_j]",
            "value": 43836,
            "unit": "ns",
            "extra": "gctime=0\nmemory=54432\nallocs=547\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_JC, σ_j]",
            "value": 57757,
            "unit": "ns",
            "extra": "gctime=0\nmemory=61376\nallocs=658\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/double-sum spin-spin",
            "value": 14582,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8688\nallocs=124\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/indexed JC H",
            "value": 1265.9,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3488\nallocs=19\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/double Σ_ij(J_ij*S_i*S_j)",
            "value": 8142.333333333333,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12904\nallocs=99\nparams={\"evals\":3,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/single Σ_i(σ_i*σ_j)",
            "value": 9584,
            "unit": "ns",
            "extra": "gctime=0\nmemory=11496\nallocs=118\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=2",
            "value": 663.0125786163522,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1568\nallocs=8\nparams={\"evals\":159,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=3",
            "value": 934.9,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1904\nallocs=9\nparams={\"evals\":20,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=4",
            "value": 1263.9,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2352\nallocs=10\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=5",
            "value": 1697.6,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2848\nallocs=11\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Ground state/3-level rewrite",
            "value": 863,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2816\nallocs=16\nparams={\"evals\":47,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Multi-mode/2-mode 6-op chain",
            "value": 1058.6,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2080\nallocs=9\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H",
            "value": 1154.8,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3264\nallocs=17\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H²",
            "value": 2697.4444444444443,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5472\nallocs=27\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H",
            "value": 1203.8,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328\nallocs=17\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H²",
            "value": 4972.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14432\nallocs=41\nparams={\"evals\":6,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H",
            "value": 1478.3,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=21\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H²",
            "value": 2537.1111111111113,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4992\nallocs=29\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
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
          "id": "cb4e31fce5c0b15221f36c191c8c86cbcf4c1bdf",
          "message": "perf: reduce TTFX (type-stability fixes + tuned precompile workload) (#223)",
          "timestamp": "2026-07-22T18:40:15+02:00",
          "tree_id": "ae08716e3b9f24be68824485cabcf599ea53b99f",
          "url": "https://github.com/qojulia/SecondQuantizedAlgebra.jl/commit/cb4e31fce5c0b15221f36c191c8c86cbcf4c1bdf"
        },
        "date": 1784738642459,
        "tool": "julia",
        "benches": [
          {
            "name": "Accumulation/Many-mode H/foldl M=16",
            "value": 11557,
            "unit": "ns",
            "extra": "gctime=0\nmemory=59472\nallocs=213\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/foldl M=24",
            "value": 25939,
            "unit": "ns",
            "extra": "gctime=0\nmemory=131408\nallocs=436\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/foldl M=8",
            "value": 2972,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12432\nallocs=63\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=16",
            "value": 1896.9,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6816\nallocs=24\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=24",
            "value": 2417,
            "unit": "ns",
            "extra": "gctime=0\nmemory=7968\nallocs=32\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=8",
            "value": 772.4537037037037,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=13\nparams={\"evals\":108,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Same-site/foldl",
            "value": 6597.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=30912\nallocs=138\nparams={\"evals\":4,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Same-site/sum",
            "value": 2270.1111111111113,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4576\nallocs=29\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=1",
            "value": 1790.7,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4128\nallocs=39\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=2",
            "value": 17767,
            "unit": "ns",
            "extra": "gctime=0\nmemory=27104\nallocs=283\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=3",
            "value": 65448,
            "unit": "ns",
            "extra": "gctime=0\nmemory=84928\nallocs=848\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=4",
            "value": 161633,
            "unit": "ns",
            "extra": "gctime=0\nmemory=186160\nallocs=1830\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=5",
            "value": 361945,
            "unit": "ns",
            "extra": "gctime=0\nmemory=364736\nallocs=3492\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, V]",
            "value": 183687,
            "unit": "ns",
            "extra": "gctime=0\nmemory=55920\nallocs=1299\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, [S, H0]]",
            "value": 337300,
            "unit": "ns",
            "extra": "gctime=0\nmemory=97184\nallocs=2206\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_Dicke, S_j]",
            "value": 42844,
            "unit": "ns",
            "extra": "gctime=0\nmemory=54432\nallocs=547\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_JC, σ_j]",
            "value": 56936,
            "unit": "ns",
            "extra": "gctime=0\nmemory=61376\nallocs=658\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/double-sum spin-spin",
            "value": 14482,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8688\nallocs=124\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/indexed JC H",
            "value": 1614.4,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3864\nallocs=25\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/double Σ_ij(J_ij*S_i*S_j)",
            "value": 8217.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12904\nallocs=99\nparams={\"evals\":4,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/single Σ_i(σ_i*σ_j)",
            "value": 9825,
            "unit": "ns",
            "extra": "gctime=0\nmemory=11496\nallocs=118\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=2",
            "value": 658.1942857142857,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1568\nallocs=8\nparams={\"evals\":175,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=3",
            "value": 936.8484848484849,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1904\nallocs=9\nparams={\"evals\":33,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=4",
            "value": 1311,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2352\nallocs=10\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=5",
            "value": 1758.7,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2848\nallocs=11\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Ground state/3-level rewrite",
            "value": 904,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2816\nallocs=16\nparams={\"evals\":42,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Multi-mode/2-mode 6-op chain",
            "value": 1056.6,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2080\nallocs=9\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H",
            "value": 1196.8,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3264\nallocs=17\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H²",
            "value": 2691.777777777778,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5472\nallocs=27\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H",
            "value": 1250.9,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328\nallocs=17\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H²",
            "value": 5099.333333333333,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14432\nallocs=41\nparams={\"evals\":6,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H",
            "value": 1540.3,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=21\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H²",
            "value": 2618.4444444444443,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4992\nallocs=29\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
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
          "id": "f4f4e45e580a6fa06332496e292c71fa37f20917",
          "message": "feat: add support for p-aware coefficients in QuantumToolbox backend (#222)",
          "timestamp": "2026-07-22T18:41:12+02:00",
          "tree_id": "8d8f0d59db5228481d772105117a9f6dbbaaad3a",
          "url": "https://github.com/qojulia/SecondQuantizedAlgebra.jl/commit/f4f4e45e580a6fa06332496e292c71fa37f20917"
        },
        "date": 1784738686550,
        "tool": "julia",
        "benches": [
          {
            "name": "Accumulation/Many-mode H/foldl M=16",
            "value": 10881,
            "unit": "ns",
            "extra": "gctime=0\nmemory=59472\nallocs=213\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/foldl M=24",
            "value": 24366,
            "unit": "ns",
            "extra": "gctime=0\nmemory=131408\nallocs=436\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/foldl M=8",
            "value": 2745.1111111111113,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12432\nallocs=63\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=16",
            "value": 1669.1,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6816\nallocs=24\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=24",
            "value": 2293.222222222222,
            "unit": "ns",
            "extra": "gctime=0\nmemory=7968\nallocs=32\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Many-mode H/sum M=8",
            "value": 707.2058823529412,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=13\nparams={\"evals\":136,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Same-site/foldl",
            "value": 6103.6,
            "unit": "ns",
            "extra": "gctime=0\nmemory=30912\nallocs=138\nparams={\"evals\":5,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Accumulation/Same-site/sum",
            "value": 2118.3333333333335,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4576\nallocs=29\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=1",
            "value": 1650.1,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4128\nallocs=39\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=2",
            "value": 17944,
            "unit": "ns",
            "extra": "gctime=0\nmemory=27104\nallocs=283\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=3",
            "value": 66565,
            "unit": "ns",
            "extra": "gctime=0\nmemory=85296\nallocs=852\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=4",
            "value": 167373,
            "unit": "ns",
            "extra": "gctime=0\nmemory=188720\nallocs=1841\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Nested JC/depth=5",
            "value": 371055,
            "unit": "ns",
            "extra": "gctime=0\nmemory=367968\nallocs=3511\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, V]",
            "value": 202809,
            "unit": "ns",
            "extra": "gctime=0\nmemory=55920\nallocs=1299\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Commutator/Schrieffer-Wolff/[S, [S, H0]]",
            "value": 371937,
            "unit": "ns",
            "extra": "gctime=0\nmemory=97184\nallocs=2206\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_Dicke, S_j]",
            "value": 43150,
            "unit": "ns",
            "extra": "gctime=0\nmemory=54432\nallocs=547\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Diagonal collapse/[H_JC, σ_j]",
            "value": 61725,
            "unit": "ns",
            "extra": "gctime=0\nmemory=61376\nallocs=658\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/double-sum spin-spin",
            "value": 15910,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8688\nallocs=124\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Simplify/indexed JC H",
            "value": 1460.7,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3864\nallocs=25\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/double Σ_ij(J_ij*S_i*S_j)",
            "value": 8562.666666666666,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12904\nallocs=99\nparams={\"evals\":3,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Indexing/Sum construction/single Σ_i(σ_i*σ_j)",
            "value": 10490,
            "unit": "ns",
            "extra": "gctime=0\nmemory=11496\nallocs=118\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=2",
            "value": 601.7374301675977,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1568\nallocs=8\nparams={\"evals\":179,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=3",
            "value": 869.3076923076923,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1904\nallocs=9\nparams={\"evals\":52,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=4",
            "value": 1169.1,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2352\nallocs=10\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Fock (c*c')^n/n=5",
            "value": 1539.9,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2848\nallocs=11\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Ground state/3-level rewrite",
            "value": 807.31,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2816\nallocs=16\nparams={\"evals\":100,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Normal Order/Multi-mode/2-mode 6-op chain",
            "value": 964.6666666666666,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2080\nallocs=9\nparams={\"evals\":21,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H",
            "value": 1066,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3264\nallocs=17\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Jaynes-Cummings/H²",
            "value": 2386.6666666666665,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5472\nallocs=27\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H",
            "value": 1109,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3328\nallocs=17\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Two cavities/H²",
            "value": 4456.857142857143,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14432\nallocs=41\nparams={\"evals\":7,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H",
            "value": 1352.6,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=21\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          },
          {
            "name": "Simplify/Λ-system/H²",
            "value": 2284.3333333333335,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4992\nallocs=29\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":3,\"time_tolerance\":0.05}"
          }
        ]
      }
    ]
  }
}