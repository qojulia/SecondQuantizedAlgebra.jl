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

# Composite basis — embed single operator
function to_numeric(op::QSym, b::QuantumOpticsBase.CompositeBasis; kwargs...)
    idx = op.space_index
    op_num = to_numeric(op, b.bases[idx]; kwargs...)
    return QuantumOpticsBase.LazyTensor(b, [idx], (op_num,))
end

# QMul
function to_numeric(m::QMul, b::QuantumOpticsBase.Basis; kwargs...)
    if isempty(m.args_nc)
        return m.arg_c * _lazy_one(b)
    end
    ops_num = [to_numeric(op, b; kwargs...) for op in m.args_nc]
    return m.arg_c * prod(ops_num)
end

# QAdd
function to_numeric(s::QAdd, b::QuantumOpticsBase.Basis; kwargs...)
    terms_num = [to_numeric(t, b; kwargs...) for t in s.arguments]
    return sum(terms_num)
end

# Number
function to_numeric(x::Number, b::QuantumOpticsBase.Basis; kwargs...)
    return x * _lazy_one(b)
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
