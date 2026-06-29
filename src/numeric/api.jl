# Public numeric API. Signatures use only QuantumInterface types (`Basis`, `StateVector`,
# `AbstractOperator`) plus the backend singletons, so `src/` never loads a numeric backend.
# The QuantumToolbox extension adds the parallel `QuantumObject` state methods.

const QuantumState = Union{StateVector, AbstractOperator}

"""
    to_numeric(op, basis [, d::AbstractDict{<:QSym}])
    to_numeric(op, state [, d::AbstractDict{<:QSym}])
    to_numeric(op, h::HilbertSpace, dims; backend, parameter, time_parameter, operators, adjoint_ops)
    to_numeric(op, basis; parameter=Dict(), time_parameter=Dict(), operators=Dict(), adjoint_ops=true)
    to_numeric(ops::AbstractVector, basis; kwargs...)

Convert a symbolic operator expression to a numeric operator on the chosen backend.

A backend is selected by loading `QuantumOpticsBase` (`QuantumOpticsBackend`) or
`QuantumToolbox` (`QuantumToolboxBackend`); the `HilbertSpace` form takes an explicit
`backend = …` (default: the single loaded backend). `d`/`operators` substitute individual
`QSym`s with custom numeric operators (missing adjoints are added when `adjoint_ops=true`).
The keyword form first substitutes scalar `parameter`s. If `time_parameter` is non-empty,
values may be numbers or functions of time and the result is the backend's **native**
time-dependent operator (`TimeDependentSum` / `QobjEvo`), callable at a time.

A `QAdd` (any sum or product) returns a **lazy** operator (`LazySum` / `QobjEvo` over a
vector-backed sum); a single operator returns the bare leaf. Materialise with
`dense`/`sparse` (QuantumOptics) or `concretize`/`sparse` (QuantumToolbox). Throws
`ArgumentError` if any prefactor cannot be reduced to a concrete `ComplexF64`.

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

to_numeric(op::Op, b::Basis, d::AbstractDict{<:QSym} = _NO_SUBS) =
    _numeric_leaf(op, NumericContext(QuantumOpticsBackend(), b, d))

to_numeric(q::QAdd, b::Basis, d::AbstractDict{<:QSym} = _NO_SUBS) =
    _to_numeric_static(q, NumericContext(QuantumOpticsBackend(), b, d))

to_numeric(x::Number, b::Basis, ::AbstractDict{<:QSym} = _NO_SUBS) =
    _to_complex(x) * numeric_identity(QuantumOpticsBackend(), b)

# --- keyword form on an explicit basis -------------------------------------------------

function to_numeric(
        op::QField, b::Basis;
        parameter = Dict(), time_parameter = Dict(),
        operators = Dict{QSym, Any}(), adjoint_ops = true,
    )
    return _to_numeric_kw(QuantumOpticsBackend(), op, b; parameter, time_parameter, operators, adjoint_ops)
end

function to_numeric(
        x::Union{Number, SymbolicUtils.BasicSymbolic}, b::Basis;
        parameter = Dict(), time_parameter = Dict(),
        operators = Dict{QSym, Any}(), adjoint_ops = true,
    )
    return _to_numeric_kw(QuantumOpticsBackend(), x, b; parameter, time_parameter, operators, adjoint_ops)
end

# --- uniform HilbertSpace form (the only one that works for both backends) -------------

function to_numeric(
        op, h::HilbertSpace, dims;
        backend::NumericBackend = _default_backend(),
        parameter = Dict(), time_parameter = Dict(),
        operators = Dict{QSym, Any}(), adjoint_ops = true,
    )
    b = numeric_basis(backend, h, dims)
    return _to_numeric_kw(backend, op, b; parameter, time_parameter, operators, adjoint_ops)
end

# --- state convenience (basis derived from the state) ----------------------------------

to_numeric(op, state::QuantumState; kwargs...) = to_numeric(op, _basis_of(state); kwargs...)
to_numeric(op::QField, state::QuantumState, d::AbstractDict{<:QSym} = _NO_SUBS) =
    to_numeric(op, _basis_of(state), d)
to_numeric(x::Number, state::QuantumState) = to_numeric(x, _basis_of(state))

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

# Leaf: build the numeric operator on the state's basis and take the backend expectation.
_numeric_average(op::QField, state, d::AbstractDict{<:QSym}) =
    numeric_expect(_backend_of(state), to_numeric(op, state, d), state)
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

# Operator expressions materialise once and expect over all states; symbolic averages map.
_numeric_average_vec(op::QField, states, d::AbstractDict{<:QSym}) =
    numeric_expect(_backend_of(first(states)), to_numeric(op, _basis_of(first(states)), d), states)
_numeric_average_vec(op, states, d::AbstractDict{<:QSym}) =
    [_numeric_average(op, ψ, d) for ψ in states]

"""
    expect(op::QField, state[, d::AbstractDict{<:QSym}])

Alias for [`numeric_average`](@ref) on operator expressions. For symbolic scalar
expressions such as `average(op)`, call [`numeric_average`](@ref) directly.
"""
expect(op::QField, state, d::AbstractDict{<:QSym} = _NO_SUBS) = numeric_average(op, state, d)
