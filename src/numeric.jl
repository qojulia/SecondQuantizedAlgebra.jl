"""
    to_numeric(op, basis; kwargs...)

Convert a symbolic operator to its numerical matrix form on the given basis.
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
function to_numeric(op::QField, state; kwargs...)
    return to_numeric(op, QuantumOpticsBase.basis(state); kwargs...)
end
function to_numeric(x::Number, state; kwargs...)
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

Compute the expectation value of a symbolic operator with a quantum state.
"""
function numeric_average(op::QField, state; kwargs...)
    op_num = to_numeric(op, state; kwargs...)
    return QuantumOpticsBase.expect(op_num, state)
end
function numeric_average(x::Number, state; kwargs...)
    return x
end

# Average expressions: unwrap and compute expectation value
function numeric_average(avg::SymbolicUtils.BasicSymbolic, state; kwargs...)
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

"""
    to_numeric(op, basis, d::Dict; kwargs...)

Convert a symbolic operator to numeric form, substituting operators found in `d`.
"""
function to_numeric(op::QSym, b::QuantumOpticsBase.Basis, d::Dict; kwargs...)
    haskey(d, op) && return d[op]
    return to_numeric(op, b; kwargs...)
end
function to_numeric(op::QField, state, d::Dict; kwargs...)
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

"""
    numeric_average(op, state, d::Dict; kwargs...)

Compute expectation value, substituting operators from `d` before numeric conversion.
"""
function numeric_average(op::QField, state, d::Dict; kwargs...)
    op_num = to_numeric(op, state, d; kwargs...)
    return QuantumOpticsBase.expect(op_num, state)
end
function numeric_average(x::Number, state, d::Dict; kwargs...)
    return x
end
function numeric_average(avg::SymbolicUtils.BasicSymbolic, state, d::Dict; kwargs...)
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

# --- to_numeric / numeric_average with ranges ---

"""
    _ranges_position(space_index, copy_index, ranges) -> Int

Compute the position in a composite basis from `ranges`.
Position = sum(ranges[1:space_index-1]) + copy_index.
"""
function _ranges_position(space_index::Int, copy_index::Int, ranges::Vector{Int})
    offset = sum(ranges[k] for k in 1:(space_index - 1); init = 0)
    return offset + copy_index
end

"""
    to_numeric(op, basis, ranges::Vector{Int})

Convert a symbolic operator to numeric form using `ranges` to map indexed operators
to positions in a composite basis. Position = sum(ranges[1:space_index-1]) + copy_index.
"""
function to_numeric(op::QSym, b::QuantumOpticsBase.CompositeBasis, ranges::Vector{Int})
    pos = _ranges_position(op.space_index, op.copy_index, ranges)
    op_num = to_numeric(op, b.bases[pos])
    return QuantumOpticsBase.LazyTensor(b, [pos], (op_num,))
end

function _term_to_numeric(c::CNum, ops::Vector{QSym}, b::QuantumOpticsBase.CompositeBasis, ranges::Vector{Int})
    if isempty(ops)
        return _to_number(c) * _lazy_one(b)
    end
    ops_num = [to_numeric(op, b, ranges) for op in ops]
    return _to_number(c) * prod(ops_num)
end

function to_numeric(s::QAdd, b::QuantumOpticsBase.CompositeBasis, ranges::Vector{Int})
    return sum(_term_to_numeric(c, ops, b, ranges) for (ops, c) in s.arguments)
end

function to_numeric(op::QField, state, ranges::Vector{Int})
    return to_numeric(op, QuantumOpticsBase.basis(state), ranges)
end

"""
    numeric_average(op, state, ranges::Vector{Int})

Compute expectation value of an indexed operator using `ranges` for basis mapping.
"""
function numeric_average(op::QField, state, ranges::Vector{Int})
    op_num = to_numeric(op, state, ranges)
    return QuantumOpticsBase.expect(op_num, state)
end
function numeric_average(avg::SymbolicUtils.BasicSymbolic, state, ranges::Vector{Int})
    if SymbolicUtils.iscall(avg) && SymbolicUtils.operation(avg) isa AvgFunc
        op = undo_average(avg)
        return numeric_average(op, state, ranges)
    end
    if SymbolicUtils.isconst(avg)
        return numeric_average(avg.val, state, ranges)
    end
    if SymbolicUtils.iscall(avg)
        f = SymbolicUtils.operation(avg)
        args = SymbolicUtils.arguments(avg)
        if f === (^)
            return numeric_average(args[1], state, ranges)^_to_number(args[2])
        end
        return f(map(a -> numeric_average(a, state, ranges), args)...)
    end
    error("numeric_average not implemented for $(typeof(avg))")
end
function numeric_average(x::Number, state, ranges::Vector{Int})
    return x
end
