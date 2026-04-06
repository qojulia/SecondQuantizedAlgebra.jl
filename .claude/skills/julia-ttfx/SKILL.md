---
name: julia-ttfx
description: >
  Diagnose and fix Time-to-First-Execute (TTFX) issues in Julia packages by finding
  and resolving method invalidations and type instabilities. Uses SnoopCompile, JET, and
  Cthulhu (via Julia MCP when available) to trace invalidation trees, analyze inference
  triggers, identify root causes like type piracy and overly broad method signatures, and
  apply targeted fixes at the source. Use this skill whenever someone mentions TTFX, latency
  on first call, invalidations, compilation time, slow package loading, "time to first plot",
  recompilation, or inference time problems in Julia — even if they don't use the exact
  term "TTFX". Also trigger when someone asks why `using MyPackage` is slow or why
  precompiled code isn't being used.
---

# Julia TTFX Diagnosis and Fix

Help users find and fix the root causes of Time-to-First-Execute (TTFX) problems in Julia packages. TTFX is the delay on first function call in a session, caused by JIT compilation. The main culprit behind *unexpectedly bad* TTFX is **method invalidation** — loading a package causes previously compiled code to be discarded and recompiled.

The goal is to eliminate invalidations and type instabilities at their source, not to paper over them with precompile workloads.

## Package versions

This skill was written and tested against:
- **SnoopCompile** v3.2.5 (SnoopCompileCore v3.1.1)
- **JET** v0.11.3 (requires Julia 1.12+; use v0.9.x for Julia 1.11)
- **Cthulhu** v3.0.2

## Running Julia code

Execute all Julia code through whichever Julia MCP tool is available in the current session (e.g., `mcp__julia__julia_eval`, `mcp__julia-server__eval`, or similar). Check the available MCP tools and use the one that can evaluate Julia code.

All analysis in this skill is non-interactive. Do not use interactive TUI tools (`@descend`, `descend`, `ascend`, `descend_code_typed`, `descend_code_warntype`, `@descend_code_typed`, `@descend_code_warntype`) — they all enter an interactive terminal loop that MCP cannot drive. Use these non-interactive alternatives instead:
- `InteractiveUtils.code_warntype(f, (types,))` for type inspection
- JET's `@report_opt` / `report_opt` for type instability detection
- Cthulhu programmatic internals for deeper inference inspection (see Step 5)

## Workflow

### Step 1: Understand the problem

Determine from context:
- Which package is the target?
- Is the issue in the package itself or a dependency?
- What's the symptom? (slow `using`, slow first call to a specific function?)
- Did it get worse recently? (new dependency, Julia version change?)

Check the project's `Project.toml` and `src/` to understand the package structure and dependencies before running analysis.

### Step 2: Record invalidations with `@snoop_invalidations`

This is almost always the right starting point. It captures every invalidation that occurs when a package loads.

**Critical ordering**: SnoopCompileCore must be loaded *before* data collection. SnoopCompile is loaded *after*, because it's a large package that itself triggers invalidations and would contaminate measurements. Also, `@snoop_invalidations` and `@snoop_inference` only capture *new* inference/invalidations — they require a fresh Julia session to get accurate results. If the MCP session already has packages loaded, restart it first.

```julia
using SnoopCompileCore

invs = @snoop_invalidations begin
    using TargetPackage
end

# Now load the analysis package
using SnoopCompile

# Parse into invalidation trees
trees = invalidation_trees(invs)

# Filter trivial cases, sort by impact (most damaging first)
trees = filter(t -> countchildren(t) > 0, trees)
sort!(trees; by=countchildren, rev=true)

# Show the worst offenders
show(trees[1:min(10, length(trees))])
```

**Pretty-print a summary table** (requires PrettyTables.jl — run `using PrettyTables` first or install it):
```julia
report_invalidations(stdout, invs; n_rows=20)
```

**Print individual trees as text:**
```julia
using AbstractTrees
print_tree(trees[1])  # full text tree of the worst offender
```

**What to look for:**
- Trees with many children — these cause the most damage
- The *root* of each tree is the method definition that triggered the cascade
- Which package defines that root method — that's where the fix belongs
- Two classes of trees: `tree.backedges` (one method covers the call) and `tree.mt_backedges` (method table scan, often from poor inference)

**Filter to a specific module:**
```julia
filtermod(TargetPackage, trees)
```

**Find path to a specific MethodInstance:**
```julia
findcaller(some_method_instance, trees)
```

### Step 3: Analyze inference time with `@snoop_inference`

If invalidations alone don't explain the TTFX, profile where inference time is spent. This requires a fresh Julia session — restart the MCP session before running (the session from Step 2 already has packages loaded).

```julia
using SnoopCompileCore

tinf = @snoop_inference begin
    TargetPackage.main_function(args...)
end

using SnoopCompile
```

If you need both invalidation data and inference data, you can collect both in a single fresh session:
```julia
using SnoopCompileCore

invs = @snoop_invalidations begin
    using TargetPackage
end

tinf = @snoop_inference begin
    TargetPackage.main_function(args...)
end

using SnoopCompile
# Now analyze both invs and tinf
```

**Text-based flamegraph** (agent-readable alternative to the visual flamegraph):
```julia
fg = flamegraph(tinf)
using AbstractTrees
print_tree(fg)  # hierarchical text representation of inference time
```

**Flat per-method timings** (tabular view of where inference time is spent):
```julia
flattened = flatten(tinf)               # sorted by exclusive time
flatten(tinf; sortby=inclusive)          # sort by inclusive time
flatten(tinf; tmin=0.001)              # filter by minimum time (seconds)
```

**Extract inference triggers** (calls made by runtime dispatch where the callee hadn't been inferred):
```julia
itrigs = inference_triggers(tinf)

# Get automated fix suggestions for each trigger
suggest(itrigs[1])

# Filter out ignorable triggers
itrigsel = [itrig for itrig in itrigs if !isignorable(suggest(itrig))]

# Group by source method and filter to your package
mtrigs = accumulate_by_source(Method, itrigs)
modtrigs = filtermod(TargetPackage, mtrigs)
```

**Build a hierarchical view:**
```julia
itree = trigger_tree(itrigs)
using AbstractTrees
print_tree(itree)
```

**Connect invalidations to inference** (only works if `invs` and `tinf` were collected in the same session — use the combined collection pattern above):
```julia
precompile_blockers(invs, tinf)
```

**Find stale (invalidated) instances in the workload:**
```julia
staleinstances(tinf)
```

### Step 4: Use JET for type instability analysis

JET's optimization analysis detects runtime dispatches, captured variables, and type instabilities — all of which create more method instances (more invalidation targets).

```julia
using JET

# Analyze a specific function — use target_modules to filter noise from Base/stdlib
@report_opt target_modules=(TargetPackage,) some_function(example_args...)
```

**JET detects three categories:**
- **Runtime dispatch** — a call that can't be resolved at compile time
- **Captured variables** (`Core.Box`) — closures capturing mutable variables, blocking inference
- **Unresolvable recursive calls** — where the compiler gives up

**JET integrates directly with SnoopCompile inference triggers:**
```julia
using JET
itrig = itrigs[5]
report_callee(itrig)     # JET analysis of the callee
report_callees(itrig)    # JET analysis of all callees
report_caller(itrig)     # JET analysis of the caller
```

**Function-form API (useful for programmatic use):**
```julia
report_opt(some_function, (ArgType1, ArgType2); target_modules=(TargetPackage,))
```

### Step 5: Drill into specifics with code_warntype and Cthulhu internals

When you need to see exactly what the compiler infers for a particular call:

**Quick type inspection:**
```julia
using InteractiveUtils
code_warntype(target_function, (ArgType1, ArgType2))
```

**Deeper inspection via Cthulhu** (renders the full type-annotated view to a string, no terminal needed):
```julia
using Cthulhu
using Cthulhu: CthulhuState, view_function, CONFIG, set_config,
               find_method_instance, generate_code_instance, lookup

interp = Base.Compiler.NativeInterpreter()
provider = Cthulhu.AbstractProvider(interp)

# Get inference data
mi = find_method_instance(provider, target_function, Tuple{ArgType1, ArgType2})
ci = generate_code_instance(provider, mi)
ci.rettype  # inferred return type — if not concrete, there's instability

result = lookup(provider, ci, false)  # false = not optimized

# Render type-annotated IR to a string
config = set_config(CONFIG; view=:typed, debuginfo=:compact, iswarn=true, optimize=false)
state = CthulhuState(provider; config, ci, mi)
io = IOBuffer()
view_function(state)(io, provider, state, result)
output = String(take!(io))
println(output)
```

Available renderers (replace `view_function(state)` call):
- `Cthulhu.cthulhu_typed(io, provider, state, result)` — typed IR with instability highlighting
- `Cthulhu.cthulhu_source(io, provider, state, result)` — annotated source view
- `Cthulhu.cthulhu_llvm(io, provider, state, result)` — LLVM IR
- `Cthulhu.cthulhu_native(io, provider, state, result)` — native assembly

These are internal APIs used stably in Cthulhu's test suite. Fall back to `code_warntype` if they error on a newer Cthulhu version.

### Step 6: Identify root causes and apply fixes

Common root causes, ranked by typical impact:

#### Type piracy
A package defines methods on types and functions it doesn't own:
```julia
# BAD — invalidates all compiled code calling show on this type
Base.show(io::IO, x::SomeOtherPkg.TheirType) = ...

# FIX: move to the package that owns TheirType, or use a wrapper type
```

#### Overly broad signatures
```julia
# BAD — matches too many types
f(x::Any) = ...
f(x::AbstractArray) = ...

# FIX: narrow to what you actually need
f(x::MySpecificType) = ...
```

#### Extending Base functions unnecessarily
```julia
# RISKY — Base.convert is called everywhere
Base.convert(::Type{MyType}, x::SomeType) = ...

# SAFER — use a constructor instead
MyType(x::SomeType) = ...
```

#### Untyped containers
```julia
# BAD — creates Vector{Any}, breaks inference
items = []

# FIX: use typed containers
items = String[]
items = Union{Int, Char}[]
```

#### Abstract struct fields
```julia
# BAD — field stored as Any
struct Foo
    x
end

# FIX: concrete type or type parameter
struct Foo{T}
    x::T
end
```

#### Captured variables in closures (Core.Box)
```julia
# BAD — r gets boxed because it's reassigned
function abmult(r::Int)
    if r < 0
        r = -r
    end
    f = x -> x * r
    return f
end

# FIX: use let-block to rebind
function abmult(r::Int)
    if r < 0
        r = -r
    end
    f = let r = r
        x -> x * r
    end
    return f
end
```

#### Overspecialization (too many compiled variants)
```julia
# If a method gets compiled for many type combinations:
function f(@nospecialize(x::SomeType))
    # reduces number of compiled method instances
end

# Julia 1.10+: also blocks inference propagation
Base.@nospecializeinfer function f(@nospecialize(x::SomeType))
    # ...
end
```

### Step 7: Verify the fix

Record the counts from Step 2/3 before making changes, then re-run in a fresh session after applying fixes:

```julia
using SnoopCompileCore

invs = @snoop_invalidations begin
    using TargetPackage
end

tinf = @snoop_inference begin
    TargetPackage.main_function(args...)
end

using SnoopCompile
trees = invalidation_trees(invs)
trees = filter(t -> countchildren(t) > 0, trees)
println("Invalidation trees: $(length(trees)), total invalidated: $(sum(countchildren, trees))")

itrigs = inference_triggers(tinf)
println("Inference triggers: $(length(itrigs))")
```

Compare these numbers against the counts recorded before the fix to confirm improvement.

## When to stop

Diminishing returns kick in quickly. If the top invalidation trees each have fewer than ~10 children and `suggest()` marks remaining triggers as ignorable, the remaining TTFX is honest compilation cost rather than fixable invalidation damage.

## Quick reference

| Function | Purpose |
|---|---|
| `invalidation_trees(invs)` | Parse raw invalidations into trees |
| `filtermod(Module, x)` | Filter trees/triggers to a module |
| `countchildren(tree)` | Number of invalidated methods in a tree |
| `findcaller(mi, trees)` | Find path to a MethodInstance |
| `inference_triggers(tinf)` | Runtime dispatch entry points |
| `suggest(itrig)` | Automated fix suggestions |
| `isignorable(suggestion)` | Check if a trigger can be ignored |
| `accumulate_by_source(Method, itrigs)` | Group triggers by method |
| `trigger_tree(itrigs)` | Hierarchical trigger view |
| `precompile_blockers(invs, tinf)` | Link invalidations to precompile failures |
| `staleinstances(tinf)` | Invalidated instances in workload |
| `flatten(tinf)` | Per-method timing breakdown |
| `flamegraph(tinf)` | Inference time tree (use with `print_tree`) |
| `report_invalidations(stdout, invs)` | Pretty-print invalidation summary |
| `report_callee(itrig)` | JET analysis of trigger callee |
| `@report_opt target_modules=(M,) f(x)` | JET optimization analysis |
| `code_warntype(f, (types,))` | Non-interactive type inspection |
| `Cthulhu.find_method_instance(provider, f, tt)` | Get MethodInstance (internal API) |
| `Cthulhu.generate_code_instance(provider, mi)` | Get CodeInstance with inferred return type (internal API) |
| `Cthulhu.cthulhu_typed(io, provider, state, result)` | Render typed IR to IO (internal API) |
