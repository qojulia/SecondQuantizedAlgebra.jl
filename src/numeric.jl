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

# QMul
function to_numeric(m::QMul, b::QuantumOpticsBase.Basis; kwargs...)
    if isempty(m.args_nc)
        return _to_number(m.arg_c) * _lazy_one(b)
    end
    ops_num = [to_numeric(op, b; kwargs...) for op in m.args_nc]
    return _to_number(m.arg_c) * prod(ops_num)
end

# QAdd
function to_numeric(s::QAdd, b::QuantumOpticsBase.Basis; kwargs...)
    terms_num = [to_numeric(t, b; kwargs...) for t in s.arguments]
    return sum(terms_num)
end

# Number
function to_numeric(x::Number, b::QuantumOpticsBase.Basis; kwargs...)
    return _to_number(x) * _lazy_one(b)
end

# Convert CNum/Num to plain Julia number for numeric evaluation.
# If the value is a literal number, extract it. If symbolic, return as-is.
_to_number(x::Number) = x
function _to_number(x::Num)
    v = Symbolics.unwrap(x)
    n = Symbolics.value(v)
    return n isa Number ? n : x  # return Num as-is if still symbolic
end
function _to_number(x::CNum)
    r = _to_number(real(x))
    i = _to_number(imag(x))
    iszero(i) && return r
    return complex(r, i)
end

# State-based dispatch
function to_numeric(op, state; kwargs...)
    return to_numeric(op, QuantumOpticsBase.basis(state); kwargs...)
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
