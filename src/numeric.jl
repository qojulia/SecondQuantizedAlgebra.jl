const QuantumState = Union{QuantumOpticsBase.StateVector, QuantumOpticsBase.AbstractOperator}

"""
    to_numeric(op, basis; kwargs...)
    to_numeric(op, state; kwargs...)
    to_numeric(op, basis, d::AbstractDict{<:QSym}; kwargs...)

Convert a symbolic operator expression to a numeric QuantumOpticsBase operator.

`basis` may be a QuantumOpticsBase `Basis` (`FockBasis`, `NLevelBasis`, `SpinBasis`,
`CompositeBasis`) or a quantum state (from which the basis is extracted). The optional
substitution dictionary `d` maps symbolic operators to custom numeric operators before
conversion. Returns a QuantumOpticsBase operator (dense, sparse, or `LazyTensor` for
composite bases).

# Examples

```jldoctest
julia> using QuantumOpticsBase

julia> h = FockSpace(:f);

julia> @qnumbers a::Destroy(h);

julia> b = FockBasis(5); ψ = fockstate(b, 2);

julia> op_num = to_numeric(a' * a, ψ);

julia> real(QuantumOpticsBase.expect(op_num, ψ)) ≈ 2
true
```

See also [`numeric_average`](@ref).
"""
function to_numeric(op::Destroy, b::QuantumOpticsBase.FockBasis; kwargs...)
    return QuantumOpticsBase.destroy(b)
end
function to_numeric(op::Create, b::QuantumOpticsBase.FockBasis; kwargs...)
    return QuantumOpticsBase.create(b)
end

function to_numeric(op::Transition, b::QuantumOpticsBase.NLevelBasis; kwargs...)
    return QuantumOpticsBase.transition(b, op.i, op.j)
end

function to_numeric(op::Pauli, b::QuantumOpticsBase.SpinBasis; kwargs...)
    b.spinnumber == 1 // 2 || throw(ArgumentError("Pauli operators require SpinBasis(1//2), got SpinBasis($(b.spinnumber))"))
    if op.axis == 1
        return QuantumOpticsBase.sigmax(b)
    elseif op.axis == 2
        return QuantumOpticsBase.sigmay(b)
    else
        return QuantumOpticsBase.sigmaz(b)
    end
end

function to_numeric(op::Spin, b::QuantumOpticsBase.SpinBasis; kwargs...)
    if op.axis == 1
        return 0.5 * QuantumOpticsBase.sigmax(b)
    elseif op.axis == 2
        return 0.5 * QuantumOpticsBase.sigmay(b)
    else
        return 0.5 * QuantumOpticsBase.sigmaz(b)
    end
end

function to_numeric(op::Position, b::QuantumOpticsBase.FockBasis; kwargs...)
    return (QuantumOpticsBase.destroy(b) + QuantumOpticsBase.create(b)) / sqrt(2)
end

function to_numeric(op::Momentum, b::QuantumOpticsBase.FockBasis; kwargs...)
    return im * (QuantumOpticsBase.create(b) - QuantumOpticsBase.destroy(b)) / sqrt(2)
end

function to_numeric(op::QSym, b::QuantumOpticsBase.CompositeBasis; kwargs...)
    idx = op.space_index
    op_num = to_numeric(op, b.bases[idx]; kwargs...)
    return QuantumOpticsBase.LazyTensor(b, [idx], (op_num,))
end

function _numeric_product(ops::Vector{QSym}, b::QuantumOpticsBase.Basis; kwargs...)
    acc = to_numeric(first(ops), b; kwargs...)
    for i in 2:length(ops)
        acc *= to_numeric(ops[i], b; kwargs...)
    end
    return acc
end

function _term_to_numeric(c::CNum, ops::Vector{QSym}, b::QuantumOpticsBase.Basis; kwargs...)
    if isempty(ops)
        return _to_number(c) * _lazy_one(b)
    end
    return _to_number(c) * _numeric_product(ops, b; kwargs...)
end

function to_numeric(s::QAdd, b::QuantumOpticsBase.Basis; kwargs...)
    return sum(_term_to_numeric(c, term.ops, b; kwargs...) for (term, c) in s.arguments)
end

function to_numeric(x::Number, b::QuantumOpticsBase.Basis; kwargs...)
    return _to_number(x) * _lazy_one(b)
end

_to_number(x::Number) = x
function _to_number(x::Num)
    v = SymbolicUtils.unwrap(x)
    n = Symbolics.value(v)
    return n isa Number ? n : x
end
function _to_number(x::CNum)
    r = _to_number(real(x))
    i = _to_number(imag(x))
    iszero(i) && return r
    return complex(r, i)
end
function _to_number(x::SymbolicUtils.BasicSymbolic)
    SymbolicUtils.isconst(x) && return _to_number(x.val)
    return _to_number(Num(x))
end

function to_numeric(op::QField, state::QuantumState; kwargs...)
    return to_numeric(op, QuantumOpticsBase.basis(state); kwargs...)
end
function to_numeric(x::Number, state::QuantumState; kwargs...)
    return to_numeric(x, QuantumOpticsBase.basis(state); kwargs...)
end

_lazy_one(b::QuantumOpticsBase.Basis) = one(b)
function _lazy_one(b::QuantumOpticsBase.CompositeBasis)
    return QuantumOpticsBase.LazyTensor(
        b, collect(1:length(b.bases)), Tuple(one(bi) for bi in b.bases)
    )
end

function to_numeric(op::QSym, b::QuantumOpticsBase.Basis, d::AbstractDict{<:QSym}; kwargs...)
    haskey(d, op) && return d[op]
    return to_numeric(op, b; kwargs...)
end
function to_numeric(op::QField, state::QuantumState, d::AbstractDict{<:QSym}; kwargs...)
    return to_numeric(op, QuantumOpticsBase.basis(state), d; kwargs...)
end

function _numeric_product(
        ops::Vector{QSym},
        b::QuantumOpticsBase.Basis,
        d::AbstractDict{<:QSym};
        kwargs...,
    )
    acc = to_numeric(first(ops), b, d; kwargs...)
    for i in 2:length(ops)
        acc *= to_numeric(ops[i], b, d; kwargs...)
    end
    return acc
end

function _term_to_numeric(
        c::CNum,
        ops::Vector{QSym},
        b::QuantumOpticsBase.Basis,
        d::AbstractDict{<:QSym};
        kwargs...,
    )
    if isempty(ops)
        return _to_number(c) * _lazy_one(b)
    end
    return _to_number(c) * _numeric_product(ops, b, d; kwargs...)
end

function to_numeric(s::QAdd, b::QuantumOpticsBase.Basis, d::AbstractDict{<:QSym}; kwargs...)
    return sum(_term_to_numeric(c, term.ops, b, d; kwargs...) for (term, c) in s.arguments)
end
function to_numeric(x::Number, b::QuantumOpticsBase.Basis, d::AbstractDict{<:QSym}; kwargs...)
    return _to_number(x) * _lazy_one(b)
end

"""
    numeric_average(op, state; kwargs...)
    numeric_average(op, states::AbstractVector; kwargs...)
    numeric_average(op, state, d::AbstractDict{<:QSym}; kwargs...)
    numeric_average(op, states::AbstractVector, d::AbstractDict{<:QSym}; kwargs...)

Compute the expectation value ``\\langle \\psi | \\hat{O} | \\psi \\rangle`` of a
symbolic operator expression. For a vector of states, returns the corresponding
vector of expectations. Averaged `BasicSymbolic` expressions (from [`average`](@ref))
are unwrapped automatically. The optional substitution dictionary `d` maps symbolic
operators to custom numeric operators before evaluation.

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
function numeric_average(op, state::QuantumState; kwargs...)
    return _numeric_average(op, state; kwargs...)
end
function numeric_average(op, states::AbstractVector; kwargs...)
    isempty(states) && throw(ArgumentError("numeric_average: states vector is empty"))
    return _numeric_average_vec(op, states; kwargs...)
end
function numeric_average(
        op,
        state::QuantumState,
        d::AbstractDict{<:QSym};
        kwargs...,
    )
    return _numeric_average(op, state, d; kwargs...)
end
function numeric_average(
        op,
        states::AbstractVector,
        d::AbstractDict{<:QSym};
        kwargs...,
    )
    isempty(states) && throw(ArgumentError("numeric_average: states vector is empty"))
    return _numeric_average_vec(op, states, d; kwargs...)
end

function _numeric_average(op::QField, state::QuantumState; kwargs...)
    op_num = to_numeric(op, state; kwargs...)
    return QuantumOpticsBase.expect(op_num, state)
end
_numeric_average(x::Number, ::QuantumState; kwargs...) = x
function _numeric_average(x::Num, state::QuantumState; kwargs...)
    return _numeric_average(SymbolicUtils.unwrap(x), state; kwargs...)
end
function _numeric_average(avg::SymbolicUtils.BasicSymbolic, state::QuantumState; kwargs...)
    if SymbolicUtils.iscall(avg) && SymbolicUtils.operation(avg) isa AvgFunc
        return _numeric_average(undo_average(avg), state; kwargs...)
    end
    if SymbolicUtils.isconst(avg)
        return _numeric_average(avg.val, state; kwargs...)
    end
    if SymbolicUtils.iscall(avg)
        f = SymbolicUtils.operation(avg)
        args = SymbolicUtils.arguments(avg)
        if f === (^)
            return _numeric_average(args[1], state; kwargs...)^_to_number(args[2])
        end
        return f((_numeric_average(a, state; kwargs...) for a in args)...)
    end
    error("numeric_average not implemented for $(typeof(avg))")
end

function _numeric_average(op::QField, state::QuantumState, d::AbstractDict{<:QSym}; kwargs...)
    op_num = to_numeric(op, state, d; kwargs...)
    return QuantumOpticsBase.expect(op_num, state)
end
_numeric_average(x::Number, ::QuantumState, ::AbstractDict{<:QSym}; kwargs...) = x
function _numeric_average(
        x::Num,
        state::QuantumState,
        d::AbstractDict{<:QSym};
        kwargs...,
    )
    return _numeric_average(SymbolicUtils.unwrap(x), state, d; kwargs...)
end
function _numeric_average(
        avg::SymbolicUtils.BasicSymbolic,
        state::QuantumState,
        d::AbstractDict{<:QSym};
        kwargs...,
    )
    if SymbolicUtils.iscall(avg) && SymbolicUtils.operation(avg) isa AvgFunc
        return _numeric_average(undo_average(avg), state, d; kwargs...)
    end
    if SymbolicUtils.isconst(avg)
        return _numeric_average(avg.val, state, d; kwargs...)
    end
    if SymbolicUtils.iscall(avg)
        f = SymbolicUtils.operation(avg)
        args = SymbolicUtils.arguments(avg)
        if f === (^)
            return _numeric_average(args[1], state, d; kwargs...)^_to_number(args[2])
        end
        return f((_numeric_average(a, state, d; kwargs...) for a in args)...)
    end
    error("numeric_average not implemented for $(typeof(avg))")
end

function _numeric_average_vec(
        op::QField,
        states::AbstractVector;
        kwargs...,
    )
    op_num = QuantumOpticsBase.sparse(to_numeric(op, first(states); kwargs...))
    return QuantumOpticsBase.expect(op_num, states)
end
function _numeric_average_vec(op, states::AbstractVector; kwargs...)
    return [_numeric_average(op, ψ; kwargs...) for ψ in states]
end

function _numeric_average_vec(
        op::QField,
        states::AbstractVector,
        d::AbstractDict{<:QSym};
        kwargs...,
    )
    op_num = QuantumOpticsBase.sparse(to_numeric(op, first(states), d; kwargs...))
    return QuantumOpticsBase.expect(op_num, states)
end
function _numeric_average_vec(op, states::AbstractVector, d::AbstractDict{<:QSym}; kwargs...)
    return [_numeric_average(op, ψ, d; kwargs...) for ψ in states]
end

"""
    expect(op, state; kwargs...)
    expect(op, state, d::AbstractDict{<:QSym}; kwargs...)

Compute the expectation value ``\\langle \\hat{O} \\rangle``. Alias for
[`numeric_average`](@ref) on symbolic operator expressions. For symbolic scalar
expressions such as `average(op)`, call [`numeric_average`](@ref) directly.
"""
expect(op::QField, state; kwargs...) = numeric_average(op, state; kwargs...)
expect(op::QField, state, d::AbstractDict{<:QSym}; kwargs...) = numeric_average(op, state, d; kwargs...)
# TODO: the four methods below are type piracy — `expect` is from QuantumOpticsBase
# and `Num`/`BasicSymbolic` are from Symbolics/SymbolicUtils. Aqua suppresses this via
# `treat_as_own = [expect]` in test/aqua_test.jl. Resolve by either dropping these
# overloads (callers can use `numeric_average` directly) or wrapping symbolic operands
# in a SQA-owned type.
# expect(avg::SymbolicUtils.BasicSymbolic, state; kwargs...) = numeric_average(avg, state; kwargs...)
# expect(avg::SymbolicUtils.BasicSymbolic, state, d::Dict; kwargs...) = numeric_average(avg, state, d; kwargs...)
# expect(x::Num, state; kwargs...) = numeric_average(x, state; kwargs...)
# expect(x::Num, state, d::Dict; kwargs...) = numeric_average(x, state, d; kwargs...)
