# Benchmark history

Every suite leaf has exactly three slash-free keys: **topic / case / legend**.
For example, `Simplify/Jaynes-Cummings/H` and
`Superoperators/Numeric conversion/to_numeric(L), 4 modes` each produce a
series in a case chart. Group related measurements under the same case so
that the visualizer can compare their legends. Keep benchmark inputs and
measurement settings unchanged when reorganizing charts.

Both runners build `benchmark_suite()` from `suite.jl`. Its validator rejects
missing or extra levels, empty groups, and slashes inside keys before tuning.
`test_suite.jl` also checks the historical rename targets against the live
suite; CI runs it and the migration tests before recording measurements.

## Migrating existing results

`hierarchy.json` records the exact old names and their replacement paths.
Unlisted names with three components remain unchanged. Do not combine
historical benchmarks whose measured workloads differ, even if their labels
look similar.

The migration uses Python 3.9 or newer and has no third-party dependencies:

```sh
python3 benchmark/migrate_history.py path/to/data.js
python3 benchmark/migrate_history.py path/to/benchmark_data.json
python3 benchmark/migrate_history.py path/to/benchmark_cpu_timestamp.json
```

It accepts the published `window.BENCHMARK_DATA` assignment, the dashboard's
JSON download, and `BenchmarkTools.save` JSON (including full trials and
estimates). Only paths change: run metadata, measurements, version headers,
and existing group tags are retained. An in-place edit creates a `.bak`
backup; rerunning on canonical data leaves the file untouched. Unknown
layouts and collisions fail before writing. Supply a second path to write a
separate migrated copy instead.

For local historical results, migrate the files in `benchmark/data/` before
using `runbenchmarks_local.jl` or `compare_two.jl` to compare old and new runs.

The published history lives on **gh-pages**, separately from this source
branch. Apply the migration to the latest `benchmark/data.js` on that branch
and deploy it together with this change, before the first new benchmark run.
Re-run against fresh gh-pages data if benchmark CI has appended results in
the meantime. Do not replace a newer history with an older migrated snapshot.
The dashboard HTML, visualizer and CSS do not need changes. Keep the backup
outside the deployed files.

```sh
python3 -m unittest discover -s benchmark -p 'test_history.py' -v
julia --project=benchmark benchmark/test_suite.jl
```
