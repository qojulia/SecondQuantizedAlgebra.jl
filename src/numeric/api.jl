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
| static | `EagerOperator()` (default) | ordinary eager backend operator |
| static | `identity` | backend's lazy assembly |
| static | another callable | backend-defined representation |
| time-dependent | `EagerOperator()` or `identity` | native time-dependent operator |

Both bundled backends use sparse storage for the default static result. Examples of explicit
conversions are `dense` for QuantumOptics and `QuantumToolbox.to_dense` for QuantumToolbox.
Other backends may choose a different eager representation. Time-dependent conversion
cannot apply a one-time materializer, so other `op_type` values are rejected.

Keyword arguments:

- `parameter`: substitutions for scalar symbolic coefficients (default: none).
- `time_parameter`: scalar values or functions of time (default: none). A non-empty mapping selects
  time-dependent conversion (`TimeDependentSum` for QuantumOptics, `QobjEvo` for
  QuantumToolbox). A value may also be a two-argument function `(p, t) -> value` of the solver
  parameter vector `p` and time `t`; the resulting `QobjEvo` is then differentiable with
  respect to `p` (for `sesolve`/`mesolve` optimal control). This `(p, t)` form is
  QuantumToolbox-only; QuantumOptics has no solver parameter and raises an `ArgumentError`. If
  any value function reads `p`, one-argument entries `t -> f(t)` and scalars are lifted to
  `(p, t)` form automatically.
- `operators`: replacements for individual symbolic operators (default: none). Missing adjoints are
  inferred when `adjoint_ops=true` (the default). Positional `d` supplies exact replacements
  without this inference.
- `op_type=EagerOperator()`: representation selection as described above.

`dims` follows the symbolic space convention: a Fock value is a cutoff, so `5` creates a
six-dimensional space. In contrast, QuantumToolbox's native `to_numeric(op, dims)`
convenience form accepts raw matrix dimensions. For QuantumToolbox composite dimensions, pass a
tuple when inference matters; a runtime `Vector`
remains supported as a convenience, but its tensor arity is not part of the Julia type.
Conversion throws `ArgumentError` when a coefficient cannot be reduced to `ComplexF64`.

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

to_numeric(op::Op, b::Basis, d::AbstractDict{<:QSym} = NO_SUBS) =
    numeric_materialize(QuantumOpticsBackend(), numeric_leaf(op, NumericContext(QuantumOpticsBackend(), b, d)), EagerOperator())
to_numeric(op::Op, b::CompositeBasis, d::AbstractDict{<:QSym} = NO_SUBS) =
    static_operator(single_qadd(CNUM_ONE, Op[op]), NumericContext(QuantumOpticsBackend(), b, d), EagerOperator())

to_numeric(q::QAdd, b::Basis, d::AbstractDict{<:QSym} = NO_SUBS) =
    static_operator(q, NumericContext(QuantumOpticsBackend(), b, d), EagerOperator())

to_numeric(x::Number, b::Basis, ::AbstractDict{<:QSym} = NO_SUBS) =
    numeric_materialize(QuantumOpticsBackend(), to_complex(x) * numeric_identity(QuantumOpticsBackend(), b), EagerOperator())

function to_numeric(
        op::QField, b::Basis;
        parameter = NO_PARAMETER, time_parameter = NO_TIME_PARAMETER,
        operators = NO_SUBS, adjoint_ops = true, op_type = EagerOperator(),
    )
    return to_numeric_kw(QuantumOpticsBackend(), op, b; parameter, time_parameter, operators, adjoint_ops, op_type)
end

function to_numeric(
        x::Union{Number, SymbolicUtils.BasicSymbolic}, b::Basis;
        parameter = NO_PARAMETER, time_parameter = NO_TIME_PARAMETER,
        operators = NO_SUBS, adjoint_ops = true, op_type = EagerOperator(),
    )
    return to_numeric_kw(QuantumOpticsBackend(), x, b; parameter, time_parameter, operators, adjoint_ops, op_type)
end

function to_numeric(
        op, h::H, dims;
        backend::B = default_backend(),
        parameter = NO_PARAMETER, time_parameter = NO_TIME_PARAMETER,
        operators = NO_SUBS, adjoint_ops = true, op_type = EagerOperator(),
    ) where {H <: HilbertSpace, B <: NumericBackend}
    b = numeric_basis(backend, h, dims)
    return to_numeric_kw(backend, op, b; parameter, time_parameter, operators, adjoint_ops, op_type)
end

to_numeric(ops::AbstractVector{T}, h::H, dims; backend::B = default_backend(), kwargs...) where {T, H <: HilbertSpace, B <: NumericBackend} =
    [to_numeric(op, h, dims; backend, kwargs...) for op in ops]

to_numeric(op, state::QuantumState; kwargs...) = to_numeric(op, numeric_basis(state); kwargs...)
to_numeric(op::QField, state::QuantumState, d::AbstractDict{<:QSym} = NO_SUBS) =
    to_numeric(op, numeric_basis(state), d)
to_numeric(x::Number, state::QuantumState) = to_numeric(x, numeric_basis(state))

to_numeric(ops::AbstractVector, b::Basis; kwargs...) = [to_numeric(op, b; kwargs...) for op in ops]

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

numeric_average(op, state::QuantumState, d::AbstractDict{<:QSym} = NO_SUBS) =
    numeric_average_impl(op, state, d)
function numeric_average(op, states::AbstractVector, d::AbstractDict{<:QSym} = NO_SUBS)
    isempty(states) && throw(ArgumentError("numeric_average: states vector is empty"))
    return numeric_average_vec(op, states, d)
end

numeric_average_impl(op::QField, state, d::AbstractDict{<:QSym}) =
    numeric_expect(numeric_backend(state), to_numeric_lazy(op, NumericContext(numeric_backend(state), numeric_basis(state), d)), state)
numeric_average_impl(x::Number, ::Any, ::AbstractDict{<:QSym}) = to_complex(x)
numeric_average_impl(x::Num, state, d::AbstractDict{<:QSym}) =
    numeric_average_impl(SymbolicUtils.unwrap(x), state, d)

function numeric_average_impl(avg::RawExpression, state, d::AbstractDict{<:QSym})
    if SymbolicUtils.hasmetadata(avg, AverageOperator) ||
            (SymbolicUtils.iscall(avg) && SymbolicUtils.operation(avg) isa AvgFunc) ||
            is_indexed_sum(avg)
        return to_complex(numeric_average_impl(undo_average(avg), state, d))
    end
    if SymbolicUtils.isconst(avg)
        v = avg.val
        v isa Op && return to_complex(numeric_average_impl(v, state, d))
        v isa QAdd && return to_complex(numeric_average_impl(v, state, d))
        return to_complex(avg)
    end
    SymbolicUtils.iscall(avg) || throw(ArgumentError("numeric_average not implemented for $(typeof(avg))"))
    f, args = SymbolicUtils.operation(avg), SymbolicUtils.arguments(avg)
    f === (^) && return to_complex(numeric_average_impl(args[1], state, d))^to_complex(args[2])
    f === (+) && return sum_avg(args, state, d)
    f === (*) && return prod_avg(args, state, d)
    f === expim && return exp(im * to_complex(numeric_average_impl(only(args), state, d)))
    throw(ArgumentError("numeric_average not implemented for operation $f"))
end

function sum_avg(args, state, d::AbstractDict{<:QSym})
    r = ComplexF64(0)
    for a in args
        r += to_complex(numeric_average_impl(a, state, d))
    end
    return r
end
function prod_avg(args, state, d::AbstractDict{<:QSym})
    r = ComplexF64(1)
    for a in args
        r *= to_complex(numeric_average_impl(a, state, d))
    end
    return r
end

numeric_average_vec(op::QField, states, d::AbstractDict{<:QSym}) =
    numeric_expect(numeric_backend(first(states)), to_numeric_lazy(op, NumericContext(numeric_backend(first(states)), numeric_basis(first(states)), d)), states)
numeric_average_vec(op, states, d::AbstractDict{<:QSym}) =
    [numeric_average(op, ψ, d) for ψ in states]

"""
    expect(op::QField, state[, d::AbstractDict{<:QSym}])

Alias for [`numeric_average`](@ref) on operator expressions. For symbolic scalar
expressions such as `average(op)`, call [`numeric_average`](@ref) directly.
"""
expect(op::QField, state, d::AbstractDict{<:QSym} = NO_SUBS) = numeric_average(op, state, d)
