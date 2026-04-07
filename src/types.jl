"""
    QField

Abstract supertype for all second-quantized operator expressions.

Subtypes:
- [`QSym`](@ref): leaf operators (e.g. [`Destroy`](@ref), [`Create`](@ref), [`Transition`](@ref))
- [`QTerm`](@ref): compound expressions (e.g. [`QAdd`](@ref))

Supports arithmetic (`+`, `-`, `*`, `^`, `/`), `adjoint`, and comparison via `==`/`isequal`.
All arithmetic eagerly applies the global [`OrderingConvention`](@ref) and returns [`QAdd`](@ref).
"""
abstract type QField end

"""
    QSym <: QField

Abstract type for fundamental (leaf) operators in the expression tree.

Every `QSym` carries four fields identifying its site:
- `name::Symbol` — display name
- `space_index::Int` — position in [`ProductSpace`](@ref) (1 for single spaces)
- `copy_index::Int` — copy label within a [`ClusterSpace`](@ref) (default 1)
- `index::Index` — symbolic summation index ([`Index`](@ref)), or `NO_INDEX`

Concrete subtypes: [`Destroy`](@ref), [`Create`](@ref), [`Transition`](@ref),
[`Pauli`](@ref), [`Spin`](@ref), [`Position`](@ref), [`Momentum`](@ref).
"""
abstract type QSym <: QField end

"""
    QTerm <: QField

Abstract type for compound noncommutative expressions built from [`QSym`](@ref) operators.

The only concrete subtype is [`QAdd`](@ref), a dict-based sum of ordered operator products.
"""
abstract type QTerm <: QField end

"""
    OrderingConvention

Abstract type for operator ordering conventions.

The global convention (set via [`set_ordering!`](@ref)) determines which commutation
relations are applied eagerly during `*`. Concrete subtypes:
- [`NormalOrder`](@ref) — creation operators left of annihilation operators (default)
- [`LazyOrder`](@ref) — no reordering; defer to [`simplify`](@ref) or [`normal_order`](@ref)
"""
abstract type OrderingConvention end

"""
    NormalOrder <: OrderingConvention

Normal ordering convention: creation operators are placed to the left of annihilation operators.

Applied commutation relations:
- **Fock**: ``[a, a^\\dagger] = 1``
- **Spin**: ``[S_j, S_k] = i\\epsilon_{jkl} S_l``
- **Phase space**: ``[p, x] = -i``

This is the default ordering convention. See also [`LazyOrder`](@ref), [`set_ordering!`](@ref).
"""
struct NormalOrder <: OrderingConvention end

"""
    LazyOrder <: OrderingConvention

Lazy ordering convention: no commutation rules are applied during multiplication (`*`).
Operator products are stored in the order they are written.

Use [`simplify`](@ref) to apply algebraic identities (Transition composition, Pauli products)
or [`normal_order`](@ref) to apply full normal-ordering commutation rules.

See also [`NormalOrder`](@ref), [`set_ordering!`](@ref).
"""
struct LazyOrder <: OrderingConvention end
