# Numeric backends and the open extension surface.
#
# `to_numeric`/`numeric_average` are backend-neutral in `src/`; every place that would
# name a concrete numeric type (a `Basis`, a lazy operator, an `expect`) goes through a
# generic hook whose methods live in a package extension (`ext/`). A backend is selected
# by a zero-field singleton so dispatch is static and the backend axis is free for
# inference. Loading `QuantumOpticsBase` or `QuantumToolbox` activates the matching
# extension, which adds the methods below.

"""
    NumericBackend

Abstract supertype for numeric-conversion backends. A backend is normally a concrete
singleton whose methods implement the exported numeric hooks. See [Adding a numeric
backend](@ref numeric-backend-interface) for the interface.

Bundled implementations are selected with [`QuantumOpticsBackend`](@ref) or
[`QuantumToolboxBackend`](@ref) after loading the corresponding package.
"""
abstract type NumericBackend end

"""
    QuantumOpticsBackend()

Select [QuantumOpticsBase.jl](https://github.com/qojulia/QuantumOpticsBase.jl). Active after
`using QuantumOpticsBase`.
"""
struct QuantumOpticsBackend <: NumericBackend end

"""
    QuantumToolboxBackend()

Select [QuantumToolbox.jl](https://github.com/qutip/QuantumToolbox.jl). Active after
`using QuantumToolbox`.
"""
struct QuantumToolboxBackend <: NumericBackend end

# --- Open extension hooks (bare generic functions; methods defined in ext/) ------------
#
# `be` is always a backend singleton. `basis` is opaque to the core: a `Basis` for
# QuantumOptics, integer dims (`Int`/`Vector{Int}`) for QuantumToolbox.

"""
    numeric_operator(be, op::Op, subbasis)

Map a single-site operator `op` to its numeric matrix on `subbasis` (the subsystem basis
for `op`'s slot). The main extension method is a closed `k === OP_*` value-branch on
`op.kind`; the separate extension point `numeric_operator(be, ::Val{kind}, op, subbasis)`
is the fallthrough for adding backend support for an existing role and throws by default.
"""
function numeric_operator end

"""
    numeric_basis(be, h::HilbertSpace, dims)

Build the full numeric basis for Hilbert space `h` with user dimensions `dims` (Fock
cutoff / spin number; ignored for spaces that carry their own level count). Returns a
`Basis` (QuantumOptics) or integer dims (QuantumToolbox).
"""
function numeric_basis end

"""
    numeric_subbasis(be, basis, slot::Int)

The subsystem basis for `slot` of a composite `basis` (the basis itself for a simple one).
"""
function numeric_subbasis end

"""
    numeric_embed(be, basis, slot::Int, leaf)

Place the single-site numeric operator `leaf` at `slot` of `basis`, returning an operator
on the full space (the leaf unchanged for a simple basis, a lazy/concrete tensor for a
composite one).
"""
function numeric_embed end

"""
    numeric_identity(be, basis)

The identity operator on `basis` (lazy on a composite basis).
"""
function numeric_identity end

"""
    numeric_num_subsystems(be, basis) -> Int

Number of tensor-product subsystems represented by `basis`. Backends use this hook to
validate indexed-site layouts before emitting any numeric operators.
"""
function numeric_num_subsystems end

"""
    numeric_assemble(be, basis, terms)

Assemble a static lazy operator from a backend-neutral `terms` list. Each entry is
`(c::ComplexF64, factors::Vector{<leaf>})` with `factors` already embedded at their slots.
The result is one concrete vector-backed lazy type for any term count.
"""
function numeric_assemble end

"""
    numeric_assemble_td(be, basis, td_terms)

Like [`numeric_assemble`](@ref) but each coefficient is `Union{ComplexF64, Function}`
(a constant or a `t -> ComplexF64`), producing the backend's native time-dependent
operator.
"""
function numeric_assemble_td end

"""
    numeric_materialize(be, op, op_type)

Materialize a static lazily assembled operator `op` according to `op_type`, applied exactly
once at the public boundary. Implementations must treat `nothing` as a request for their
ordinary eager operator and `identity` as a request for `op` unchanged. Other callables are
backend-defined. The eager representation need not be sparse, although both bundled
backends choose sparse storage.

This hook is not used by [`numeric_average`](@ref), which consumes the lazy assembly, or by
time-dependent conversion, which returns the backend's native time-dependent operator.
"""
function numeric_materialize end

"""
    numeric_expect(be, numop, state)

Expectation value `⟨state| numop |state⟩` (a `ComplexF64`), or a `Vector{ComplexF64}` when
`state` is a vector of states.
"""
function numeric_expect end

"""
    numeric_backend(state) -> NumericBackend

Return the numeric backend owning `state`. Third-party backends must implement this hook
for their state types to support [`numeric_average`](@ref) and state-based [`to_numeric`](@ref).
"""
function numeric_backend end

# The one-argument `numeric_basis(state)` is the public state-to-basis counterpart of the
# three-argument Hilbert-space builder declared above.
"""
    numeric_basis(state)

Return the backend basis/dimensions carried by `state`. Third-party backends must implement
this together with [`numeric_backend`](@ref) for state-based numeric conversion.
"""
numeric_basis(s::Union{StateVector, AbstractOperator}) = basis(s)

function _check_product_dims(h::ProductSpace, dims)
    nspaces = length(h.spaces)
    ndims = try
        length(dims)
    catch
        throw(
            ArgumentError(
                "dims for a ProductSpace must be an indexable collection with $nspaces entries",
            ),
        )
    end
    ndims == nspaces || throw(
        ArgumentError(
            "dims has $ndims entries, but the ProductSpace has $nspaces subspaces",
        ),
    )
    return dims
end

"""
    _default_backend()

The backend to use when none is given and no backend object pins one. Resolves to the sole
loaded backend extension, erroring helpfully when zero or both are loaded.
"""
function _default_backend()
    m = @__MODULE__
    qob = Base.get_extension(m, :SecondQuantizedAlgebraQuantumOpticsBaseExt) !== nothing
    qtb = Base.get_extension(m, :SecondQuantizedAlgebraQuantumToolboxExt) !== nothing
    qob && !qtb && return QuantumOpticsBackend()
    qtb && !qob && return QuantumToolboxBackend()
    if qob && qtb
        throw(
            ArgumentError(
                "both numeric backends are loaded; pass `backend = QuantumOpticsBackend()` " *
                    "or `backend = QuantumToolboxBackend()` explicitly.",
            ),
        )
    end
    throw(
        ArgumentError(
            "no numeric backend loaded; run `using QuantumOpticsBase` or `using QuantumToolbox` " *
                "to enable `to_numeric`/`numeric_average`.",
        ),
    )
end
