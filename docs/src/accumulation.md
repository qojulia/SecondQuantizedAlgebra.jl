```@meta
CurrentModule = SecondQuantizedAlgebra
```

# Building Large Sums Efficiently

Assembling a Hamiltonian, a Lindbladian, or a set of equations of motion usually means adding together many operator terms. SecondQuantizedAlgebra.jl accelerates that assembly: `sum` and `reduce(+, …)` over a vector of expressions accumulate in place, turning an `O(n^2)` chain of additions into an `O(n)` pass. For a many-mode or many-atom model this is the difference between an assembly step that scales gracefully and one that dominates your build time.

## Why repeated addition is slow

A [`QAdd`](@ref) stores its terms in a dictionary so that like terms collect automatically. Every binary `+` has *value semantics*: because its two operands are objects the caller still holds, the operation cannot mutate either one, so it copies the whole backing dictionary before inserting. That copy is the right thing to do for a single `a + b`, but it is wasteful when you fold a long list:

```math
\big(\dots((t_1 + t_2) + t_3) + \dots + t_n\big)
```

Each step copies the growing accumulator, so building an ``n``-term sum by repeated `+` copies ``1 + 2 + \dots + n`` terms in total, which is ``O(n^2)``.

The accumulator threads **one** shared dictionary through the entire reduction and materializes a single `QAdd` at the end. It reuses the exact same in-place insertion the binary `+` uses internally, so the result is byte-identical to the `+`-chain, but it pays the copy cost only once.

## A coupled-cavity chain

Consider a driven Bose-Hubbard / coupled-cavity array of ``M`` sites, a workhorse model in circuit QED and photonic lattices:

```math
H = \sum_{k=1}^{M} \left[ -\Delta\, a_k^\dagger a_k + \frac{U}{2}\, a_k^\dagger a_k^\dagger a_k a_k + \Omega\,(a_k + a_k^\dagger) \right] - J \sum_{k=1}^{M-1} \left( a_k^\dagger a_{k+1} + a_{k+1}^\dagger a_k \right).
```

Each site carries a detuning ``\Delta``, an on-site interaction ``U``, a coherent drive ``\Omega``, and nearest-neighbor hopping ``J``. Build it by summing over comprehensions:

```@example accum
using SecondQuantizedAlgebra

@variables Δ U Ω J

M = 6
h = reduce(⊗, (FockSpace(Symbol(:c, k)) for k in 1:M))
a = [Destroy(h, Symbol(:a, k), k) for k in 1:M]

onsite = sum([-Δ*a[k]'*a[k] + (U/2)*a[k]'*a[k]'*a[k]*a[k] + Ω*(a[k] + a[k]') for k in 1:M])
hopping = sum([-J*(a[k]'*a[k+1] + a[k+1]'*a[k]) for k in 1:(M-1)])

H = onsite + hopping
nothing # hide
```

```@example accum
println("H has ", length(H.arguments), " terms")
```

The `sum([… for k in …])` form is all you need: the comprehension produces a `Vector{QAdd}`, and `sum` over that vector takes the fast in-place path automatically. The result is identical to folding the same terms with `+`:

```@example accum
terms = [-Δ*a[k]'*a[k] + (U/2)*a[k]'*a[k]'*a[k]*a[k] + Ω*(a[k] + a[k]') for k in 1:M]
sum(terms) == foldl(+, terms)
```

The speedup grows with system size, because the wasted copying it removes is the part that scales quadratically. Indicative timings for the chain above (building the full `Vector{QAdd}` of terms, then reducing it):

| Sites ``M`` | terms | `foldl(+, …)` | `sum(…)` | speedup |
|:-----------:|:-----:|:-------------:|:--------:|:-------:|
| 8           | 46    | ~27 µs        | ~5 µs    | ~5×     |
| 16          | 94    | ~110 µs       | ~8 µs    | ~13×    |
| 24          | 142   | ~244 µs       | ~11 µs   | ~22×    |

## The fast-path rule

The accelerated methods are dispatched on an *array* of `QAdd`s, so the shape of how you call them matters:

```@example accum
H1 = sum([Ω*(a[k] + a[k]') for k in 1:M])    # fast: a Vector{QAdd} is built, then accumulated
H2 = reduce(+, [Ω*(a[k] + a[k]') for k in 1:M]) # fast: same path as sum
H1 == H2
```

!!! warning "Use brackets, not a bare generator or a `+=` loop"
    A bare generator cannot be intercepted by element type, and an explicit accumulation loop calls the copying binary `+` on every step. Both stay on the slow `O(n^2)` path:

    ```julia
    H = sum(Ω*(a[k] + a[k]') for k in 1:M)        # SLOW: bare generator → generic +-fold
    H = zero(QAdd); for k in 1:M; H += Ω*(a[k]+a[k]'); end  # SLOW: copies on every +=
    ```

    Wrap the generator in brackets (`sum([… for …])`) or collect the terms into a vector first. The extra `O(M)` intermediate vector is negligible next to the quadratic saving.

## Where it helps

The accumulator speeds up the *additive assembly* of an expression from many parts. The applications that benefit most are the ones whose term count grows with system size:

- **Many-mode and lattice Hamiltonians**: Bose-Hubbard chains, coupled-cavity arrays, optomechanical networks.
- **Atom and spin ensembles**: Tavis-Cummings and Dicke-type models, collective spin operators ``S_\pm = \sum_i \sigma_i^\pm``.
- **Lindbladians with many channels**: a master-equation generator that sums one dissipator per decay or dephasing channel.
- **Equations of motion / cumulant systems**: each derived ``\dot{\langle O\rangle}`` is a sum of many contributions, accumulated across a large operator set.
- **Truncated expansions**: Magnus, Dyson, Baker-Campbell-Hausdorff, and perturbation series, where the truncation order sets the number of terms.

For a two-mode toy model the difference is invisible; for tens or hundreds of sites it is the dominant assembly cost.

!!! note "Products are a separate cost"
    The accumulator only accelerates *addition*. The cost of the operator products themselves (high powers, dense normal-ordering) is intrinsic to canonicalization and is not affected. The product path is deliberately left as repeated `*`; see the [Developer Documentation](@ref) for the rationale.

## Under the hood

Internally the accumulation is exposed through the [MutableArithmetics.jl](https://github.com/jump-dev/MutableArithmetics.jl) interface on a small builder type. This is what makes `sum` and `reduce(+, …)` fast, and it also lets the [`@rewrite`](https://jump.dev/MutableArithmetics.jl/stable/) macro and manual `operate!!` loops compose with `QAdd` correctly:

```julia
import MutableArithmetics as MA

acc = sum([Ω*(a[k] + a[k]') for k in 1:M])
manual = MA.@rewrite sum(Ω*(a[k] + a[k]') for k in 1:M)
acc == manual   # true
```

!!! tip "For plain sums, prefer the bracketed `sum`"
    `@rewrite` threads the same in-place accumulator, but its macro-expansion and buffer setup add overhead that, for a straightforward sum, outweighs the saving. Reach for the bracketed `sum([… for …])` (or `reduce(+, [...])`) as the everyday fast path; `@rewrite` is for when you are already composing a larger in-place arithmetic expression.

The binary `+`, `-`, and scalar `*` keep their value semantics unchanged: only the reduction drivers that own a transient accumulator take the in-place path. The analytical sum [`Σ`](@ref) over a symbolic [`Index`](@ref) is orthogonal and already accumulates in place; see [Symbolic Sums and Indices](@ref).
