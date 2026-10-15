const QuantumState = Union{StateVector, AbstractOperator}

"""
    to_numeric(op, h::HilbertSpace, dims; backend, parameter, time_parameter, operators, adjoint_ops, op_type)
    to_numeric(op, basis; kwargs...)
    to_numeric(op, state; kwargs...)
    to_numeric(op, basis_or_state, operators)
    to_numeric(ops::AbstractVector, args...; kwargs...)

Matrix of an operator expression on a truncated basis, for a numerical backend.

The first form builds the basis from the symbolic space `h` and works with every backend.
`dims` holds one entry per factor: the photon cutoff for a mode, so `5` gives six states,
and the spin length for a spin. The backend is `QuantumOpticsBackend()` or
`QuantumToolboxBackend()`, and may be left out when only one of the two packages is loaded.
The other forms take a basis of the backend, or a state whose basis is used; QuantumToolbox
also accepts its own dimension tuple, which counts states rather than photons. The third
argument `operators` replaces operators by given matrices, as the keyword does with
`adjoint_ops = false`. A vector of expressions gives a vector of matrices.

Keywords:

- `parameter`: values of the symbolic parameters. A coefficient that does not reduce to a
  number throws an `ArgumentError`.
- `time_parameter`: values that depend on time, as functions `t -> value` or numbers; the
  time variable itself is given as `t => identity`. A non-empty `time_parameter` makes the
  result time dependent, a `TimeDependentSum` for QuantumOptics and a `QobjEvo` for
  QuantumToolbox. With QuantumToolbox, a function `(p, t) -> value` of the solver parameters
  `p` makes the `QobjEvo` differentiable in `p`, for optimal control, and may be mixed with
  the other entries.
- `operators`: matrices on the full basis to use for given operators in place of the
  standard ones.
- `adjoint_ops = true`: also replace the adjoint of each operator in `operators` by the
  adjoint matrix.
- `op_type`: the default gives a sparse matrix. `identity` keeps the backend's lazy
  representation, which saves memory on large product spaces, and another function, such as
  `dense` or `QuantumToolbox.to_dense`, is applied to the result. A time-dependent result is
  always the backend's time-dependent operator, so only the default and `identity` are
  accepted.

Formal operator functions throw an `ArgumentError`; expand them with [`taylor`](@ref) first.

# Examples
```jldoctest
julia> using QuantumOpticsBase

julia> h = FockSpace(:cavity) ⊗ NLevelSpace(:atom, 2);

julia> @qnumbers a::Destroy(h, 1) σ::Transition(h, 1, 2, 2);

julia> @variables g;

julia> Hn = to_numeric(g * (a' * σ + σ' * a), h, (5, 2); backend = QuantumOpticsBackend(), parameter = Dict(g => 0.1));

julia> size(Hn)
(12, 12)
```

See also [`numeric_average`](@ref), [`QuantumOpticsBackend`](@ref), [`QuantumToolboxBackend`](@ref).
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
    numeric_average(x, state[, operators])
    numeric_average(x, states::AbstractVector[, operators])

Expectation value of an operator expression or of an expression of averages in a numerical state.

`x` is either an operator expression or a result of [`average`](@ref), with every parameter
already replaced by a number. The state is a ket or a density matrix of the backend, and its
basis fixes the truncation; a vector of states gives a vector of values. `operators` maps
operators to the matrices to use in place of the standard ones.

# Examples
```jldoctest
julia> using QuantumOpticsBase

julia> a = Destroy(FockSpace(:cavity), :a);

julia> ψ = coherentstate(FockBasis(20), 0.5);

julia> numeric_average(a' * a, ψ) ≈ 0.25
true

julia> numeric_average(average(a) + average(a'), ψ) ≈ 1
true
```

See also [`expect`](@ref expect(::QField, ::Any)), [`to_numeric`](@ref), [`average`](@ref).
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
    avg === Symbolics.IM && return ComplexF64(im)
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
    expect(op, state[, operators])

Expectation value of an operator expression, the same as [`numeric_average`](@ref).

The name matches `expect` of QuantumOptics.jl and QuantumToolbox.jl. Unlike
`numeric_average`, it does not accept expressions of averages.

# Examples
```jldoctest
julia> using QuantumOpticsBase

julia> a = Destroy(FockSpace(:cavity), :a);

julia> expect(a' * a, fockstate(FockBasis(5), 2)) ≈ 2
true
```

See also [`numeric_average`](@ref).
"""
expect(op::QField, state, d::AbstractDict{<:QSym} = NO_SUBS) = numeric_average(op, state, d)
