const QuantumState = Union{QuantumOpticsBase.StateVector, QuantumOpticsBase.AbstractOperator}

"""
    to_numeric(op, basis; kwargs...)
    to_numeric(op, state; kwargs...)
    to_numeric(op, basis, d::Dict; kwargs...)

Convert a symbolic operator expression to its numerical matrix representation
using QuantumOpticsBase.

# Arguments
- `op` — a [`QSym`](@ref), [`QAdd`](@ref), or `Number`
- `basis` — a QuantumOpticsBase `Basis` (e.g. `FockBasis`, `NLevelBasis`, `SpinBasis`,
  `CompositeBasis`), or a quantum state (from which the basis is extracted)
- `d::Dict` — optional operator substitution dictionary (keys: `QSym`, values: numeric operators)

Returns a QuantumOpticsBase operator (dense, sparse, or `LazyTensor` for composite bases).

# Examples
```julia
using QuantumOpticsBase
h = FockSpace(:f)
@qnumbers a::Destroy(h)
b = FockBasis(10)
op_num = to_numeric(a' * a, b)    # number operator as 11×11 matrix
```

See also [`numeric_average`](@ref).
"""
function to_numeric(op::Destroy, b::QuantumOpticsBase.FockBasis; kwargs...)
    return QuantumOpticsBase.destroy(b)
end
function to_numeric(op::Create, b::QuantumOpticsBase.FockBasis; kwargs...)
    return QuantumOpticsBase.create(b)
end

# Transition
function to_numeric(op::Transition, b::QuantumOpticsBase.NLevelBasis; kwargs...)
    return QuantumOpticsBase.transition(b, op.i, op.j)
end

# Pauli — requires spin-1/2 basis
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

# Spin
function to_numeric(op::Spin, b::QuantumOpticsBase.SpinBasis; kwargs...)
    if op.axis == 1
        return 0.5 * QuantumOpticsBase.sigmax(b)
    elseif op.axis == 2
        return 0.5 * QuantumOpticsBase.sigmay(b)
    else
        return 0.5 * QuantumOpticsBase.sigmaz(b)
    end
end

# Position: X = (a + a†) / √2
function to_numeric(op::Position, b::QuantumOpticsBase.FockBasis; kwargs...)
    return (QuantumOpticsBase.destroy(b) + QuantumOpticsBase.create(b)) / sqrt(2)
end

# Momentum: P = i(a† - a) / √2
function to_numeric(op::Momentum, b::QuantumOpticsBase.FockBasis; kwargs...)
    return im * (QuantumOpticsBase.create(b) - QuantumOpticsBase.destroy(b)) / sqrt(2)
end

# Composite basis — embed single operator
function to_numeric(op::QSym, b::QuantumOpticsBase.CompositeBasis; kwargs...)
    idx = op.space_index
    op_num = to_numeric(op, b.bases[idx]; kwargs...)
    return QuantumOpticsBase.LazyTensor(b, [idx], (op_num,))
end

# Internal: convert a single term (ops, prefactor) to numeric
function _term_to_numeric(c::CNum, ops::Vector{QSym}, b::QuantumOpticsBase.Basis; kwargs...)
    if isempty(ops)
        return _to_number(c) * _lazy_one(b)
    end
    ops_num = [to_numeric(op, b; kwargs...) for op in ops]
    return _to_number(c) * prod(ops_num)
end

# QAdd
function to_numeric(s::QAdd, b::QuantumOpticsBase.Basis; kwargs...)
    return sum(_term_to_numeric(c, ops, b; kwargs...) for (ops, c) in s.arguments)
end

# Number
function to_numeric(x::Number, b::QuantumOpticsBase.Basis; kwargs...)
    return _to_number(x) * _lazy_one(b)
end

# Convert CNum/Num to plain Julia number for numeric evaluation.
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

# State-based dispatch
function to_numeric(op::QField, state::QuantumState; kwargs...)
    return to_numeric(op, QuantumOpticsBase.basis(state); kwargs...)
end
function to_numeric(x::Number, state::QuantumState; kwargs...)
    return to_numeric(x, QuantumOpticsBase.basis(state); kwargs...)
end

# Lazy identity helper
_lazy_one(b::QuantumOpticsBase.Basis) = one(b)
function _lazy_one(b::QuantumOpticsBase.CompositeBasis)
    return QuantumOpticsBase.LazyTensor(
        b, collect(1:length(b.bases)), Tuple(one(bi) for bi in b.bases)
    )
end

"""
    numeric_average(op, state; kwargs...)
    numeric_average(op, state, d::Dict; kwargs...)

Compute the expectation value ``\\langle \\psi | \\hat{O} | \\psi \\rangle`` of a symbolic
operator expression `op` with a QuantumOpticsBase quantum state.

Converts `op` to numeric form via [`to_numeric`](@ref), then calls
`QuantumOpticsBase.expect`. Also handles averaged `BasicSymbolic` expressions
by first calling [`undo_average`](@ref).

# Examples
```julia
using QuantumOpticsBase
h = FockSpace(:f)
@qnumbers a::Destroy(h)
b = FockBasis(10)
ψ = coherentstate(b, 2.0)
numeric_average(a' * a, ψ)    # ≈ 4.0
```

See also [`to_numeric`](@ref), [`average`](@ref).
"""
function numeric_average(op::QField, state::QuantumState; kwargs...)
    op_num = to_numeric(op, state; kwargs...)
    return QuantumOpticsBase.expect(op_num, state)
end
function numeric_average(op::QField, states::AbstractVector; kwargs...)
    isempty(states) && throw(ArgumentError("numeric_average: states vector is empty"))
    op_num = QuantumOpticsBase.sparse(to_numeric(op, first(states); kwargs...))
    return QuantumOpticsBase.expect(op_num, states)
end
function numeric_average(x::Number, state::QuantumState; kwargs...)
    return x
end
function numeric_average(x::Num, state::QuantumState; kwargs...)
    return numeric_average(SymbolicUtils.unwrap(x), state; kwargs...)
end

# Average expressions: unwrap and compute expectation value
function numeric_average(avg::SymbolicUtils.BasicSymbolic, state::QuantumState; kwargs...)
    if SymbolicUtils.iscall(avg) && SymbolicUtils.operation(avg) isa AvgFunc
        op = undo_average(avg)
        return numeric_average(op, state; kwargs...)
    end
    if SymbolicUtils.isconst(avg)
        return numeric_average(avg.val, state; kwargs...)
    end
    if SymbolicUtils.iscall(avg)
        f = SymbolicUtils.operation(avg)
        args = SymbolicUtils.arguments(avg)
        if f === (^)
            return numeric_average(args[1], state; kwargs...)^_to_number(args[2])
        end
        return f(numeric_average.(args, Ref(state); kwargs...)...)
    end
    error("numeric_average not implemented for $(typeof(avg))")
end

# --- to_numeric / numeric_average with Dict parameter substitution ---

function to_numeric(op::QSym, b::QuantumOpticsBase.Basis, d::Dict; kwargs...)
    haskey(d, op) && return d[op]
    return to_numeric(op, b; kwargs...)
end
function to_numeric(op::QField, state::QuantumState, d::Dict; kwargs...)
    return to_numeric(op, QuantumOpticsBase.basis(state), d; kwargs...)
end

function _term_to_numeric(c::CNum, ops::Vector{QSym}, b::QuantumOpticsBase.Basis, d::Dict; kwargs...)
    if isempty(ops)
        return _to_number(c) * _lazy_one(b)
    end
    ops_num = [to_numeric(op, b, d; kwargs...) for op in ops]
    return _to_number(c) * prod(ops_num)
end

function to_numeric(s::QAdd, b::QuantumOpticsBase.Basis, d::Dict; kwargs...)
    return sum(_term_to_numeric(c, ops, b, d; kwargs...) for (ops, c) in s.arguments)
end
function to_numeric(x::Number, b::QuantumOpticsBase.Basis, d::Dict; kwargs...)
    return _to_number(x) * _lazy_one(b)
end

function numeric_average(op::QField, state::QuantumState, d::Dict; kwargs...)
    op_num = to_numeric(op, state, d; kwargs...)
    return QuantumOpticsBase.expect(op_num, state)
end
function numeric_average(op::QField, states::AbstractVector, d::Dict; kwargs...)
    isempty(states) && throw(ArgumentError("numeric_average: states vector is empty"))
    op_num = to_numeric(op, first(states), d; kwargs...)
    return QuantumOpticsBase.expect(op_num, states)
end
function numeric_average(x::Number, state::QuantumState, d::Dict; kwargs...)
    return x
end
function numeric_average(x::Num, state::QuantumState, d::Dict; kwargs...)
    return numeric_average(SymbolicUtils.unwrap(x), state, d; kwargs...)
end
function numeric_average(avg::SymbolicUtils.BasicSymbolic, state::QuantumState, d::Dict; kwargs...)
    if SymbolicUtils.iscall(avg) && SymbolicUtils.operation(avg) isa AvgFunc
        op = undo_average(avg)
        return numeric_average(op, state, d; kwargs...)
    end
    if SymbolicUtils.isconst(avg)
        return numeric_average(avg.val, state, d; kwargs...)
    end
    if SymbolicUtils.iscall(avg)
        f = SymbolicUtils.operation(avg)
        args = SymbolicUtils.arguments(avg)
        if f === (^)
            return numeric_average(args[1], state, d; kwargs...)^_to_number(args[2])
        end
        return f(numeric_average.(args, Ref(state), Ref(d); kwargs...)...)
    end
    error("numeric_average not implemented for $(typeof(avg))")
end

"""
    expect(op, state; kwargs...)
    expect(op, state, d::Dict; kwargs...)

Compute the expectation value ``\\langle \\hat{O} \\rangle``. Alias for
[`numeric_average`](@ref); imported from `QuantumOpticsBase` so the same name
works for both numeric and symbolic operands.
"""
expect(op::QField, state; kwargs...) = numeric_average(op, state; kwargs...)
expect(op::QField, state, d::Dict; kwargs...) = numeric_average(op, state, d; kwargs...)
expect(avg::SymbolicUtils.BasicSymbolic, state; kwargs...) = numeric_average(avg, state; kwargs...)
expect(avg::SymbolicUtils.BasicSymbolic, state, d::Dict; kwargs...) = numeric_average(avg, state, d; kwargs...)
expect(x::Num, state; kwargs...) = numeric_average(x, state; kwargs...)
expect(x::Num, state, d::Dict; kwargs...) = numeric_average(x, state, d; kwargs...)
