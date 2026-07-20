# Public numeric API. Signatures use only QuantumInterface types (`Basis`, `StateVector`,
# `AbstractOperator`) plus the backend singletons, so `src/` never loads a numeric backend.
# The QuantumToolbox extension adds the parallel `QuantumObject` state methods.

const QuantumState = Union{StateVector, AbstractOperator}

"""
    to_numeric(op, h::HilbertSpace, dims; backend, kwargs...)
    to_numeric(op, basis [, d]; kwargs...)
    to_numeric(op, state [, d]; kwargs...)
    to_numeric(ops::AbstractVector, ...; kwargs...)

Convert a symbolic operator expression to a numeric operator on the chosen backend.

The first form is backend-neutral and is recommended for portable code. Load
`QuantumOpticsBase` or `QuantumToolbox` and pass `QuantumOpticsBackend()` or
`QuantumToolboxBackend()`; `backend` may be omitted when exactly one bundled backend is
loaded. Backend-native basis, state, and dimension forms are convenience overloads.

The representation contract is independent of the expression shape:

| Input | `op_type` | Result |
|:--|:--|:--|
| static | `nothing` (default) | ordinary eager backend operator |
| static | `identity` | backend's lazy assembly |
| static | another callable | backend-defined representation |
| time-dependent | `nothing` or `identity` | native time-dependent operator |

Both bundled backends use sparse storage for the default static result. Examples of explicit
conversions are `dense` for QuantumOptics and `QuantumToolbox.to_dense` for QuantumToolbox.
Other backends may choose a different eager representation. Time-dependent conversion
cannot apply a one-time materializer, so other `op_type` values are rejected.

Keyword arguments:

- `parameter=Dict()`: substitutions for scalar symbolic coefficients.
- `time_parameter=Dict()`: scalar values or functions of time. A non-empty mapping selects
  time-dependent conversion (`TimeDependentSum` for QuantumOptics, `QobjEvo` for
  QuantumToolbox).
- `operators=Dict()`: replacements for individual symbolic operators. Missing adjoints are
  inferred when `adjoint_ops=true` (the default). Positional `d` supplies exact replacements
  without this inference.
- `op_type=nothing`: representation selection as described above.

`dims` follows the symbolic space convention: a Fock value is a cutoff, so `5` creates a
six-dimensional space. In contrast, QuantumToolbox's native `to_numeric(op, dims)`
convenience form accepts raw matrix dimensions. Conversion throws `ArgumentError` when a
coefficient cannot be reduced to `ComplexF64`.

# Examples

```jldoctest
julia> using QuantumOpticsBase

julia> h = FockSpace(:f);

julia> @qnumbers a::Destroy(h);

julia> b = FockBasis(5); ψ = fockstate(b, 2);

julia> real(expect(to_numeric(a' * a, b), ψ)) ≈ 2
true
```

See also [`numeric_average`](@ref).
"""
function to_numeric end

# --- leaf / QAdd / scalar on an explicit basis (positional, no kwargs) -----------------
# The positional forms always materialize `sparse` (`op_type` is only on the keyword form).

to_numeric(op::Op, b::Basis, d::AbstractDict{<:QSym} = _NO_SUBS) =
    numeric_materialize(QuantumOpticsBackend(), _numeric_leaf(op, NumericContext(QuantumOpticsBackend(), b, d)), nothing)

to_numeric(q::QAdd, b::Basis, d::AbstractDict{<:QSym} = _NO_SUBS) =
    numeric_materialize(QuantumOpticsBackend(), _to_numeric_static(q, NumericContext(QuantumOpticsBackend(), b, d)), nothing)

to_numeric(x::Number, b::Basis, ::AbstractDict{<:QSym} = _NO_SUBS) =
    numeric_materialize(QuantumOpticsBackend(), _to_complex(x) * numeric_identity(QuantumOpticsBackend(), b), nothing)

# --- keyword form on an explicit basis -------------------------------------------------

function to_numeric(
        op::QField, b::Basis;
        parameter = Dict(), time_parameter = Dict(),
        operators = Dict{QSym, Any}(), adjoint_ops = true, op_type = nothing,
    )
    return _to_numeric_kw(QuantumOpticsBackend(), op, b; parameter, time_parameter, operators, adjoint_ops, op_type)
end

function to_numeric(
        x::Union{Number, SymbolicUtils.BasicSymbolic}, b::Basis;
        parameter = Dict(), time_parameter = Dict(),
        operators = Dict{QSym, Any}(), adjoint_ops = true, op_type = nothing,
    )
    return _to_numeric_kw(QuantumOpticsBackend(), x, b; parameter, time_parameter, operators, adjoint_ops, op_type)
end

# --- uniform HilbertSpace form (the only one that works for both backends) -------------

function to_numeric(
        op, h::HilbertSpace, dims;
        backend::NumericBackend = _default_backend(),
        parameter = Dict(), time_parameter = Dict(),
        operators = Dict{QSym, Any}(), adjoint_ops = true, op_type = nothing,
    )
    b = numeric_basis(backend, h, dims)
    return _to_numeric_kw(backend, op, b; parameter, time_parameter, operators, adjoint_ops, op_type)
end

# Vector of operators on the uniform HilbertSpace form (both backends).
to_numeric(ops::AbstractVector, h::HilbertSpace, dims; backend::NumericBackend = _default_backend(), kwargs...) =
    [to_numeric(op, h, dims; backend, kwargs...) for op in ops]

# --- state convenience (basis derived from the state) ----------------------------------

to_numeric(op, state::QuantumState; kwargs...) = to_numeric(op, numeric_basis(state); kwargs...)
to_numeric(op::QField, state::QuantumState, d::AbstractDict{<:QSym} = _NO_SUBS) =
    to_numeric(op, numeric_basis(state), d)
to_numeric(x::Number, state::QuantumState) = to_numeric(x, numeric_basis(state))

# --- vector of operators ---------------------------------------------------------------

to_numeric(ops::AbstractVector, b::Basis; kwargs...) = [to_numeric(op, b; kwargs...) for op in ops]

# =======================================================================================

"""
    numeric_average(op, state[, d::AbstractDict{<:QSym}]) -> ComplexF64
    numeric_average(op, states::AbstractVector[, d::AbstractDict{<:QSym}]) -> Vector

Expectation value ``\\langle \\psi | \\hat{O} | \\psi \\rangle`` of a symbolic operator
expression. Averaged `BasicSymbolic` expressions (from [`average`](@ref)) are unwrapped
automatically. `d` substitutes symbolic operators with custom numeric operators before
evaluation.

# Examples

```jldoctest
julia> using QuantumOpticsBase

julia> h = FockSpace(:f);

julia> @qnumbers a::Destroy(h);

julia> b = FockBasis(5); ψ = fockstate(b, 2);

julia> real(numeric_average(a' * a, ψ)) ≈ 2
true
```

See also [`to_numeric`](@ref), [`average`](@ref).
"""
function numeric_average end

numeric_average(op, state::QuantumState, d::AbstractDict{<:QSym} = _NO_SUBS) =
    _numeric_average(op, state, d)
function numeric_average(op, states::AbstractVector, d::AbstractDict{<:QSym} = _NO_SUBS)
    isempty(states) && throw(ArgumentError("numeric_average: states vector is empty"))
    return _numeric_average_vec(op, states, d)
end

# Leaf: assemble the lazy numeric operator on the state's basis (never materialized, so a
# `LazyKet` state works) and take the backend expectation.
_numeric_average(op::QField, state, d::AbstractDict{<:QSym}) =
    numeric_expect(numeric_backend(state), _to_numeric_lazy(op, NumericContext(numeric_backend(state), numeric_basis(state), d)), state)
_numeric_average(x::Number, ::Any, ::AbstractDict{<:QSym}) = _to_complex(x)
_numeric_average(x::Num, state, d::AbstractDict{<:QSym}) =
    _numeric_average(SymbolicUtils.unwrap(x), state, d)

# `arguments(::BasicSymbolic)` is statically `Any`; the `_to_complex` wraps reseal each
# recursive result to `ComplexF64`.
function _numeric_average(avg::SymbolicUtils.BasicSymbolic, state, d::AbstractDict{<:QSym})
    if SymbolicUtils.hasmetadata(avg, AverageOperator) ||
            (SymbolicUtils.iscall(avg) && SymbolicUtils.operation(avg) isa AvgFunc) ||
            is_indexed_sum(avg)
        return _to_complex(_numeric_average(undo_average(avg), state, d))
    end
    SymbolicUtils.isconst(avg) && return _to_complex(_numeric_average(avg.val, state, d))
    SymbolicUtils.iscall(avg) || throw(ArgumentError("numeric_average not implemented for $(typeof(avg))"))
    f, args = SymbolicUtils.operation(avg), SymbolicUtils.arguments(avg)
    f === (^) && return _to_complex(_numeric_average(args[1], state, d))^_to_complex(args[2])
    f === (+) && return _sum_avg(args, state, d)
    f === (*) && return _prod_avg(args, state, d)
    throw(ArgumentError("numeric_average not implemented for operation $f"))
end

function _sum_avg(args, state, d::AbstractDict{<:QSym})
    r = ComplexF64(0)
    for a in args
        r += _to_complex(_numeric_average(a, state, d))
    end
    return r
end
function _prod_avg(args, state, d::AbstractDict{<:QSym})
    r = ComplexF64(1)
    for a in args
        r *= _to_complex(_numeric_average(a, state, d))
    end
    return r
end

# Operator expressions assemble lazily once and expect over all states (the vector
# `numeric_expect` materialises internally); symbolic averages map.
_numeric_average_vec(op::QField, states, d::AbstractDict{<:QSym}) =
    numeric_expect(numeric_backend(first(states)), _to_numeric_lazy(op, NumericContext(numeric_backend(first(states)), numeric_basis(first(states)), d)), states)
_numeric_average_vec(op, states, d::AbstractDict{<:QSym}) =
    [_numeric_average(op, ψ, d) for ψ in states]

"""
    expect(op::QField, state[, d::AbstractDict{<:QSym}])

Alias for [`numeric_average`](@ref) on operator expressions. For symbolic scalar
expressions such as `average(op)`, call [`numeric_average`](@ref) directly.
"""
expect(op::QField, state, d::AbstractDict{<:QSym} = _NO_SUBS) = numeric_average(op, state, d)
