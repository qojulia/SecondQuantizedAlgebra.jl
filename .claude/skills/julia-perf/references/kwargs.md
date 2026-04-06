# Keyword Arguments (kwargs) Performance

Keyword arguments introduce overhead that the compiler often cannot eliminate.
Avoid them in performance-sensitive code paths.

## Why kwargs are slow

1. **No specialization** — Julia specializes on positional argument types but keyword
   arguments are handled via `NamedTuple`s that may not be fully inferred, especially
   when splatted (`kwargs...`).

2. **Dynamic dispatch** — Splatted `kwargs...` prevents the compiler from knowing which
   keywords are present at compile time, blocking inlining and causing dynamic dispatch.

3. **Allocation** — Splatting keyword arguments (`func(; kwargs...)`) can allocate a
   `NamedTuple` or `pairs(::NamedTuple)` at runtime, adding GC pressure.

4. **Invalidation risk** — Methods with keyword arguments generate hidden "keyword sorter"
   methods internally. These can cause method invalidations when packages are loaded,
   hurting TTFX and precompilation effectiveness.

## When it matters

- **Hot inner loops** — avoid kwargs; use positional arguments.
- **User-facing API functions** — kwargs are fine; the ergonomic benefit outweighs the
  negligible cost when the function body does real work.

## Recommended pattern

Expose kwargs at the API boundary, forward to a positional-argument inner function:

```julia
# User-facing API with kwargs for convenience
function solve(problem; tol=1e-8, maxiter=1000)
    _solve(problem, tol, maxiter)  # forward to positional-arg inner function
end

# Performance-critical inner function with positional args only
function _solve(problem, tol, maxiter)
    # hot loop here — fully specializable and inlineable
end
```

## Anti-patterns

```julia
# BAD: splatted kwargs forwarding — blocks inference
function outer(x; kwargs...)
    inner(x; kwargs...)
end

# GOOD: explicit forwarding
function outer(x; tol=1e-8, maxiter=1000)
    inner(x, tol, maxiter)
end
```

```julia
# BAD: kwargs in a hot loop
for i in 1:1_000_000
    step(x; dt=0.01, method=:rk4)
end

# GOOD: positional args in hot loop
for i in 1:1_000_000
    _step(x, 0.01, :rk4)
end
```
