"""
    NumericBackend

Abstract supertype for numeric-conversion backends.

A backend is normally a concrete singleton whose methods implement the numeric hooks. See
[Adding a numeric backend](@ref) for the interface.

Bundled implementations are selected with [`QuantumOpticsBackend`](@ref) or
[`QuantumToolboxBackend`](@ref) after loading the corresponding package.
"""
abstract type NumericBackend end

"""
    QuantumOpticsBackend()

Numerical backend built on [QuantumOpticsBase.jl](https://github.com/qojulia/QuantumOpticsBase.jl).

Available after `using QuantumOpticsBase`. Operators become `Operator`s on QuantumOptics
bases, and time-dependent ones `TimeDependentSum`s, ready for the solvers of QuantumOptics.jl.

# Examples
```jldoctest
julia> using QuantumOpticsBase

julia> h = FockSpace(:cavity); a = Destroy(h, :a);

julia> to_numeric(a' * a, h, 4; backend = QuantumOpticsBackend()) ≈ number(FockBasis(4))
true
```

See also [`QuantumToolboxBackend`](@ref), [`to_numeric`](@ref).
"""
struct QuantumOpticsBackend <: NumericBackend end

"""
    QuantumToolboxBackend()

Numerical backend built on [QuantumToolbox.jl](https://github.com/qutip/QuantumToolbox.jl).

Available after `using QuantumToolbox`. Operators become `QuantumObject`s, and
time-dependent ones `QobjEvo`s, ready for the solvers of QuantumToolbox.jl.

# Examples
```jldoctest
julia> import QuantumToolbox

julia> h = FockSpace(:cavity); a = Destroy(h, :a);

julia> to_numeric(a' * a, h, 4; backend = QuantumToolboxBackend()) ≈ QuantumToolbox.num(5)
true
```

See also [`QuantumOpticsBackend`](@ref), [`to_numeric`](@ref).
"""
struct QuantumToolboxBackend <: NumericBackend end

"""
    numeric_operator(be, op::Op, subbasis)

Map a single-site operator `op` to its numeric matrix on `subbasis` (the subsystem basis
for `op`'s slot). Backends implement it as a closed `k === OP_*` value-branch on `op.kind`.
"""
function numeric_operator end

operator_leaf(be::NumericBackend, op::Op, subbasis) =
    is_custom(op) ? realized_leaf(be, op, subbasis) : numeric_operator(be, op, subbasis)

function realized_leaf(be::NumericBackend, op::Op, subbasis)
    data = kind_data(op)
    isempty(data.realization) &&
        throw(ArgumentError("operator kind $(data.name) has no realization, so it cannot be converted to a numeric operator"))
    leaf = 0.0im * numeric_operator(be, data.probe, subbasis)
    for (c, ops) in data.realization[op.l1]
        leaf = leaf + c * realized_product(be, ops, subbasis)
    end
    return leaf
end

function realized_product(be::NumericBackend, ops::Vector{Op}, subbasis)
    isempty(ops) && return numeric_identity(be, subbasis)
    m = numeric_operator(be, ops[1], subbasis)
    for k in 2:length(ops)
        m = m * numeric_operator(be, ops[k], subbasis)
    end
    return m
end

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

Time-dependent counterpart of [`numeric_assemble`](@ref).

Each coefficient is `Union{ComplexF64, Function}`: a constant, a `t -> ComplexF64`, or a
solver-parameter-aware `(p, t) -> ComplexF64`. The result is the backend's native
time-dependent operator.
"""
function numeric_assemble_td end

"""
    EagerOperator()

The default `op_type`: request the backend's ordinary eager operator.
"""
struct EagerOperator end

"""
    numeric_materialize(be, op, op_type)

Materialize a static lazily assembled operator `op` according to `op_type`, applied exactly
once at the public boundary. Implementations must treat [`EagerOperator`](@ref) as a request
for their ordinary eager operator and `identity` as a request for `op` unchanged. Other callables are
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

"""
    numeric_basis(state)

Return the backend basis/dimensions carried by `state`. Third-party backends must implement
this together with [`numeric_backend`](@ref) for state-based numeric conversion.
"""
numeric_basis(s::Union{StateVector, AbstractOperator}) = basis(s)

function check_product_dims(h::ProductSpace, dims)
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
    default_backend()

The backend to use when none is given and no backend object pins one. Resolves to the sole
loaded backend extension, erroring helpfully when zero or both are loaded.
"""
function default_backend()
    m = @__MODULE__
    qob = Base.get_extension(m, :SecondQuantizedAlgebraQuantumOpticsBaseExt) isa Module
    qtb = Base.get_extension(m, :SecondQuantizedAlgebraQuantumToolboxExt) isa Module
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
