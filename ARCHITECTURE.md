# Architecture

## Overview

An expression flows through four layers.

```text
constructors ──► Op ──*,+──► QAdd  (canonical polynomial, Dict{QTerm, Coeff})
                                │
                 sin, cos, expim│                     taylor(expr, 0:n)
                                ▼                            │
                              QExpr (formal, ordered tree) ──┘ back to QAdd
                                │
                    average     ▼          to_numeric / numeric_average
             QAdd ──────► symbolic scalars           QAdd ──────► backend operator
                          (SymbolicUtils terms)                   (package extension)
```

Every operator is the single concrete struct `Op`, and its role (annihilation, transition, Pauli, …) is a `kind` tag. A polynomial expression is a `QAdd`: a dictionary from terms to coefficients, where a term is an ordered operator product together with the index inequalities that scope it. Products and sums are brought into canonical form the moment they are built. Coefficients are `Coeff` values: exact numbers, parameter polynomials, or raw SymbolicUtils expressions. Formal operator functions live in a separate tree, `QExpr`, which never enters the polynomial code. Averages turn a `QAdd` into SymbolicUtils terms for equation-of-motion packages such as QuantumCumulants.jl, and the numeric layer converts a `QAdd` to matrices through a backend defined in a package extension. Maps on operators, `SuperOperator`, `Lindbladian` and `SuperExpr`, form a parallel layer whose two sides are canonical operator products.

## Code map

| Directory | Contents |
| :-- | :-- |
| `src/operators/` | Hilbert spaces, the `Op` struct and its role constructors, and the operator hooks every algebra rule goes through |
| `src/expressions/` | coefficients (exact numbers, polynomials, `Coeff`), `QTerm`, `QAdd`, `QExpr`, `Index` and `Σ`, name interning |
| `src/algebra/` | the canonicalization pipeline, `simplify`/`normal_order`/`substitute`, sum accumulation, Weyl ordering, unitary transformations |
| `src/superoperators/` | `SuperOperator`, `Lindbladian`, `SuperExpr`, their shared index scope and frame transforms |
| `src/average.jl` | `average`, `undo_average` and the average and sum nodes |
| `src/numeric/` | the backend-neutral conversion core, superoperator conversion and the backend interface |
| `ext/` | the QuantumOpticsBase and QuantumToolbox backends |
| `src/printing/` | Unicode and LaTeX rendering |

## Invariants

### Canonical form

**Every product is canonical when it is built.** The result of `*` is already sorted, reduced and commuted, so `normal_order` and `simplify` are idempotent on it and like terms collect at once. A new code path that builds terms must end in the canonicalization pipeline.

**Rewrites come in three classes.** Local identities on a provably shared site (transition composition, Pauli products) always fire. Commutation swaps fire only when the site relation is known. Rewrites that multiply the number of terms, such as completeness, fire only on request. This is why the site comparison distinguishes `Equal` from `Undetermined`, and why it must keep doing so.

**The pipeline runs reduce, then commute, then reduce.** The first reduce composes transitions and Pauli operators, which compose rather than commute. Commutation then handles bosonic, spin and quadrature pairs, and the final reduce catches compositions that a commutator residual creates next to a neighbour.

**The ground-state projector $\sigma^{gg}$ stays atomic.** Rewriting it as $1 - \sum_{k \neq g} \sigma^{kk}$ multiplies a product with $k$ such factors into $n^k$ terms before anything can cancel, and a mean-field solver would then carry $n - 1$ dependent moments instead of one. `expand_completeness` is the explicit handle.

**Free indices stay undetermined.** Two operators on the same space with different free indices may or may not act on the same site, so no rule fires between them and they keep their written order. Only a `Σ` or `assume_distinct_index` resolves them.

**The canonical form does not depend on declaration order.** Operator and index names are interned to integer ids in insertion order, so ordering by id would make the canonical form depend on which name was declared first. All orderings go through the alphabetical rank of the name instead.

**A pass that forks must copy.** The pipeline streams terms through nested passes without materializing intermediate lists. A pass with one output may mutate the operator vector in place; a pass with several outputs copies it for every branch but one. A missed copy gives two terms that share one vector, a silent aliasing bug.

**Diagonal splitting substitutes before sorting.** When a product meets a summed index, the term splits into the off-diagonal part and the part where the two indices coincide. The coinciding part substitutes the index in the unsorted product: sorting first can put two same-site operators in the wrong physical order, so that $\sigma_j^{12}\sigma_i^{21}$ at $i = j$ would give $\sigma^{22}$ instead of $\sigma^{11}$.

**Every coincidence is counted once.** A summed index that may equal several free indices takes the first one it equals, so the diagonal for the $k$-th excludes the earlier ones. Several binders are split one at a time, each inside the parts of the previous one; splitting each against the original term would count a joint coincidence such as $i = k, j = l$ once per binder.

### Terms and indices

**A term is its operators together with its inequalities.** Like terms merge only when both the operator product and the inequality constraints match, so $\sum_i a_i$ and $\sum_{i \neq j} a_i$ remain separate entries. Zero coefficients are never stored.

**Every stored inequality refers to something the term can see.** A constraint must involve an index of one of the term's operators, of its coefficient, or of the enclosing sum. The `QAdd` constructor drops the others, since they would only prevent like terms from merging.

**The algebra never creates an `Index`.** Every index in an output was constructed by the user, so the user can match on it, substitute it and give it limits. Hence a product of two sums over the same bound index throws instead of renaming one of them, and $\sum_i \sigma_i + \sum_j \sigma_j$ stays two terms.

**Binary `+` copies; reductions do not.** Both operands of `+` are still held by the caller, so `+` copies the accumulated dictionary and a chain of $n$ additions costs $O(n^2)$. `sum` and `reduce(+, …)` over an array instead thread one accumulator through the fold and build the result once; the result must stay identical to the chain of `+`. Products are left as repeated `*`, because their cost is the distributive expansion itself, which an accumulator cannot save.

**Interned ids never leave the session.** The name and range tables in `intern.jl` are written only at construction, from one thread, and their ids follow insertion order. Ids are never serialized, and no SymbolicUtils symbol built from them may be cached at precompile time: it would be a stale copy of the symbol the new session hashconses.

### Expression layers

**`QExpr` never touches the polynomial path.** Polynomial work stays in canonical `QAdd` leaves, and a formal function is never expanded implicitly. `taylor` is the only way back to `QAdd`. Static `to_numeric` instead evaluates matrix functions on complete materialized arguments at the chosen finite basis, after all free parameters have been replaced by static values. A nonempty `time_parameter` is unsupported. `average` and `numeric_average` still refuse an unlowered `QExpr`. A change to `QExpr` must not cost the polynomial algebra anything.

### Superoperators

**Maps are not `QField`s.** Were `SuperOperator` or `Lindbladian` a `QField`, every operator method would accept them, and `a * L` or `average(L)` would return nonsense instead of an error. The three map types share the internal `Union` `AnyMap`, since `AbstractSuperOperator` belongs to QuantumInterface.

**Both sides share one scope.** In $\sum_i \sigma_i \rho \sigma_i^\dagger$ the same $i$ binds both sides, so a term stores one inequality list and one summation scope for the pair. The diagonal split runs on the joint unsorted term, and every binder the result no longer mentions is counted, repeatedly, since eliminating one can make another's count known. Without this, `transpose(Σ(dissipator(a_i'), i))(1)` would be $1 - N$ instead of zero.

**Ambiguous indices throw.** Composing maps that bind the same index, applying a map to an operator with that index free, or a side carrying its own `Σ` has no unique reading, and the algebra never mints an index to choose one.

**Lindbladian equality is structural.** Normalization moves jump constants into the Hamiltonian, numeric factors and unimodular phases into the rate, and drops the identity part of $H$, so equal spellings of one generator compare equal. A cross term is not expanded: `dissipator(a + b)` and its four-channel expansion are unequal as Lindbladians and equal as `SuperOperator`s, which keeps the user's channels through printing and frames. `σᵍᵍ` stays atomic here as well.

**GKSL extraction is exact by default.** Float maps carry roundoff that grows as the basis approaches linear dependence, so no fixed tolerance both keeps small genuine terms and accepts every valid float map; `rtol` is an opt-in, measured in map units so that rescaling the basis cannot hide a term. Trace preservation is checked first, so a trace defect is reported as one.

**Channel order is reproducible.** `channels(L)` sorts by `canonical_isless`, since `Dict` order depends on hashing and printing and positional access must not.

**The return type never depends on values.** A side with a formal function makes a `SuperExpr`, and a `SuperExpr` whose formal part is empty stays one until `SuperOperator(X)` is called.

### Coefficients

**`Complex{Num}` exists only at the boundary.** Coefficients are `Coeff` internally; `Complex{Num}` appears only where values leave the package, as in `get_prefactor` or `to_num`. In particular, `average` must not build a `Complex{Num}`: Symbolics turns it into an opaque `complex(re, im)` node that `simplify` cannot cancel and that generated code cannot call.

**Exactness is a type.** An exact value is a `GaussianRational`, and a `ComplexF64` is always inexact. Exact with exact stays exact, and anything combined with a float becomes a float, so no operation has to inspect a value.

**One tier per value.** Exact arithmetic runs on `Int` and is redone on `BigInt` when it overflows, so its steps must not mutate their inputs. A result returns to `Int` when it fits, so the same value never exists in both tiers.

**Radicals are canonical.** A radical is a product of prime square roots (or other roots) with exponents in $(0, 1)$, so $\sqrt{6}$, $\sqrt{2}\sqrt{3}$ and $\sqrt{24}/2$ are the same value with the same hash.

**Only plain products become compound atoms.** A function of an atom, or of a product of real atoms such as `cos(ω*t)`, is a polynomial atom, which keeps its arithmetic on the fast path. Scaled or summed arguments stay raw, since `sin(2θ)` must reach the trigonometric reduction and `exp(-r)` must cancel `1/exp(r)`.

**`GaussianRational` is not a `Number`.** Methods of a new `Number` subtype on `==`, `hash` and arithmetic invalidate SymbolicUtils code compiled for abstract numbers.

### Averages

**An average is real when its operator is Hermitian.** The average and sum nodes override `maketerm` to derive their type from the operator, since the default rebuild in `substitute` would silently turn a real average back into a complex one. The index scope of a sum node is an argument, not metadata, because SymbolicUtils ignores metadata in `isequal` and two sums with different scopes would otherwise cancel.

### Orderings

**Three orderings, three purposes.** The site comparison decides which neighbours may swap and is deliberately partial. The display order in `sorted_arguments` is total but may tie distinct operators. Only `order_key` is faithful to identity, `order_key(a) == order_key(b)` exactly when `isequal(a, b)`, and downstream packages use it through `canonical_isless` to choose canonical representatives. That comparison is deliberately not `Base.isless`, since `a + b < 2a` would read as a comparison of operator values. None of them is `hash`, which is not stable across Julia versions.

### Unitary transformations

**The affine action is the single representation.** A `UnitaryTransform` stores $X \mapsto L X + c$ as blocks, one per algebra: bosonic, phase space, spin and Pauli, finite-level. Each block is inverted with the formula of its own group, so `inv` never inverts a symbolic matrix, and the rewrite rules are compiled from the action. Caller-supplied matrices are checked for structure but never for unitarity, so no tolerance or projection is applied.

### Types and performance

**Struct fields are concrete, and `Op` and `Index` are bits types.** `Vector{Op}` is then an inline buffer that is hashed and compared on integers.

**The numeric core is backend-neutral.** `src/numeric/` never names a backend type, and the type of an assembled operator does not depend on the number of terms.

## Adding an operator kind

Every operator is an `Op`, and its role is the `kind` field of type `OpKind`. The algebra never inspects an operator directly. It talks to operators through six functions in `src/operators/operators.jl`, each a single method on `(::Op, ::Op)` that branches on `kind`:

| Function | Returns | Purpose |
| :-- | :-- | :-- |
| `site_compare(a, b, ne)` | `Less`, `Equal`, `Undetermined` or `Greater` | site relation of two operators, which drives the sort |
| `can_commute(a, b)` | `Bool` | whether a same-site pair may swap without a residual |
| `commute_pair(a, b)` | `(c1, r1, c2, r2, n)` | the residual $ab - ba = c_1 r_1 + c_2 r_2$, where each `r` is one operator for `n = 1` and the identity for `n = 0` |
| `may_reduce(a, b)` | `Bool` | whether `reduce_pair` can do anything, a cheap pre-check |
| `reduce_pair(a, b)` | `(kind, op, c)` | a local identity: the pair becomes one operator, a number, or stays |
| `ground_state_expand(op)` | `(is_gs, g, n, site)` | the completeness data of a ground-state projector |

Two residual slots exist because a collective commutator $[S^{ij}, S^{kl}] = \delta_{jk} S^{il} - \delta_{li} S^{kj}$ has two terms; the other kinds fill only the first.

A new kind touches these places:

1. **The kind.** An `OP_*` member of `OpKind` in `src/operators/op.jl`, an `is_*` predicate, and a `@public` entry.
2. **The constructor.** A function in the file of its Hilbert space that validates the arguments and returns an `Op`. Levels or axes go in the shared fields `l1`, `l2`, `g`, `nlev`; a kind that needs none leaves them zero. Write the constructor for the single space only; its `ProductSpace` methods are one-liners through `on_factor`, which checks the factor and moves the operator onto it.
3. **The algebra.** A branch in each function of the table that the kind needs. Kinds on one space that interact, like annihilation and creation, form a family in `site_family`.
4. **Adjoint.** A branch of `adjoint` in `src/operators/operators.jl`. `order_key` reads the shared fields and needs no change, as long as distinct operators of the kind differ in those fields.
5. **Numerics.** A branch of `numeric_operator` in every backend under `ext/`, or a deliberate `ArgumentError` where a backend has no representation.
6. **Printing.** The Unicode printer in `src/printing/printing.jl` and the LaTeX recipe in `src/printing/latexify_recipes.jl`.
7. **Tests.** Construction in `test/fundamentals/`, the algebra rules in `test/algebra/`, the conversion in `test/numeric/`, and the output in `test/presentation/`.

An extra payload field on `Op`, for instance a `Union{Nothing, …}`, is not an option: `Op` must stay a concrete bits type.
