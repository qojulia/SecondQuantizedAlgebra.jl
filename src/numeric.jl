const QuantumState = Union{QuantumOpticsBase.StateVector, QuantumOpticsBase.AbstractOperator}

# `Union{}` value keeps the haskey-fallback branch concrete-typed.
const _NO_SUBS = Dict{QSym, Union{}}()

"""
    to_numeric(op, basis [, d::AbstractDict{<:QSym}])
    to_numeric(op, state [, d::AbstractDict{<:QSym}])

Convert a symbolic operator expression to a numeric QuantumOpticsBase operator.
`d` substitutes individual `QSym`s with custom numeric operators. Throws
`ArgumentError` if any prefactor cannot be reduced to a concrete `ComplexF64`.

# Examples

```jldoctest
julia> using QuantumOpticsBase

julia> h = FockSpace(:f);

julia> @qnumbers a::Destroy(h);

julia> b = FockBasis(5); ψ = fockstate(b, 2);

julia> real(QuantumOpticsBase.expect(to_numeric(a' * a, ψ), ψ)) ≈ 2
true
```

See also [`numeric_average`](@ref).
"""
function to_numeric end

# Per-basis isa-chain (≤4 alternatives stays in Julia's union-split budget).
function to_numeric(op::QSym, b::QuantumOpticsBase.FockBasis)
    op isa Destroy  && return QuantumOpticsBase.destroy(b)
    op isa Create   && return QuantumOpticsBase.create(b)
    op isa Position && return (QuantumOpticsBase.destroy(b) + QuantumOpticsBase.create(b)) / sqrt(2)
    op isa Momentum && return im * (QuantumOpticsBase.create(b) - QuantumOpticsBase.destroy(b)) / sqrt(2)
    throw(ArgumentError("$(typeof(op)) does not act on a FockBasis"))
end
function to_numeric(op::QSym, b::QuantumOpticsBase.NLevelBasis)
    op isa Transition && return QuantumOpticsBase.transition(b, op.i, op.j)
    throw(ArgumentError("$(typeof(op)) does not act on an NLevelBasis"))
end
function to_numeric(op::QSym, b::QuantumOpticsBase.SpinBasis)
    if op isa Pauli
        b.spinnumber == 1 // 2 || throw(ArgumentError("Pauli operators require SpinBasis(1//2), got SpinBasis($(b.spinnumber))"))
        op.axis == 1 && return QuantumOpticsBase.sigmax(b)
        op.axis == 2 && return QuantumOpticsBase.sigmay(b)
        return QuantumOpticsBase.sigmaz(b)
    end
    if op isa Spin
        op.axis == 1 && return 0.5 * QuantumOpticsBase.sigmax(b)
        op.axis == 2 && return 0.5 * QuantumOpticsBase.sigmay(b)
        return 0.5 * QuantumOpticsBase.sigmaz(b)
    end
    throw(ArgumentError("$(typeof(op)) does not act on a SpinBasis"))
end
function to_numeric(op::QSym, b::QuantumOpticsBase.CompositeBasis)
    return QuantumOpticsBase.LazyTensor(b, [op.space_index], (to_numeric(op, b.bases[op.space_index]),))
end

function to_numeric(op::QSym, b::QuantumOpticsBase.Basis, d::AbstractDict{<:QSym})
    haskey(d, op) && return d[op]
    return to_numeric(op, b)
end

# Scalar-term and operator-product branches return different operator types
# (Diagonal identity vs SparseMatrixCSC), so the conditional stays in-loop.
function to_numeric(s::QAdd, b::QuantumOpticsBase.Basis, d::AbstractDict{<:QSym} = _NO_SUBS)
    iter = iterate(s.arguments)
    iter === nothing && return _to_complex(_CNUM_ZERO) * _lazy_one(b)
    (term, c), st = iter
    result = isempty(term.ops) ? _to_complex(c) * _lazy_one(b) : _to_complex(c) * _product(term.ops, b, d)
    while true
        next = iterate(s.arguments, st)
        next === nothing && return result
        (term, c), st = next
        result += isempty(term.ops) ? _to_complex(c) * _lazy_one(b) : _to_complex(c) * _product(term.ops, b, d)
    end
    return
end

function _product(ops::Vector{QSym}, b::QuantumOpticsBase.Basis, d::AbstractDict{<:QSym} = _NO_SUBS)
    acc = to_numeric(first(ops), b, d)
    for i in 2:length(ops)
        acc *= to_numeric(ops[i], b, d)
    end
    return acc
end

to_numeric(x::Number, b::QuantumOpticsBase.Basis, ::AbstractDict{<:QSym} = _NO_SUBS) = _to_complex(x) * _lazy_one(b)
to_numeric(op::QField, state::QuantumState, d::AbstractDict{<:QSym} = _NO_SUBS) = to_numeric(op, QuantumOpticsBase.basis(state), d)
to_numeric(x::Number, state::QuantumState) = to_numeric(x, QuantumOpticsBase.basis(state))

_lazy_one(b::QuantumOpticsBase.Basis) = one(b)
_lazy_one(b::QuantumOpticsBase.CompositeBasis) = QuantumOpticsBase.LazyTensor(b, collect(1:length(b.bases)), Tuple(one(bi) for bi in b.bases))

# One method (union-split budget) routing every input through `convert ∘ Complex`,
# the only pattern that infers to ComplexF64 from `Any` after `Symbolics.value`.
function _to_complex(x)
    x isa ComplexF64   && return x
    x isa Complex{Num} && return convert(ComplexF64, Complex(Symbolics.value(real(x)), Symbolics.value(imag(x))))
    x isa Complex      && return convert(ComplexF64, x)
    x isa Num          && return convert(ComplexF64, Complex(Symbolics.value(x), false))
    if x isa SymbolicUtils.BasicSymbolic
        SymbolicUtils.isconst(x) || throw(ArgumentError("cannot reduce symbolic expression $x to a concrete number"))
        return convert(ComplexF64, Complex(x.val, false))
    end
    return convert(ComplexF64, Complex(x, false))
end

"""
    numeric_average(op, state[, d::AbstractDict{<:QSym}]) -> ComplexF64
    numeric_average(op, states::AbstractVector[, d::AbstractDict{<:QSym}]) -> Vector

Expectation value ``\\langle \\psi | \\hat{O} | \\psi \\rangle`` of a symbolic
operator expression. Averaged `BasicSymbolic` expressions (from
[`average`](@ref)) are unwrapped automatically. `d` substitutes symbolic
operators with custom numeric operators before evaluation.

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

numeric_average(op, state::QuantumState, d::AbstractDict{<:QSym} = _NO_SUBS) = _numeric_average(op, state, d)
function numeric_average(op, states::AbstractVector, d::AbstractDict{<:QSym} = _NO_SUBS)
    isempty(states) && throw(ArgumentError("numeric_average: states vector is empty"))
    return _numeric_average_vec(op, states, d)
end

_numeric_average(op::QField, state::QuantumState, d::AbstractDict{<:QSym}) = ComplexF64(QuantumOpticsBase.expect(to_numeric(op, state, d), state))
_numeric_average(x::Number, ::QuantumState, ::AbstractDict{<:QSym}) = _to_complex(x)
_numeric_average(x::Num, state::QuantumState, d::AbstractDict{<:QSym}) = _numeric_average(SymbolicUtils.unwrap(x), state, d)

# `arguments(::BasicSymbolic)` is statically `Any`; the `_to_complex` wraps
# reseal each recursive result to `ComplexF64`.
function _numeric_average(avg::SymbolicUtils.BasicSymbolic, state::QuantumState, d::AbstractDict{<:QSym})
    if SymbolicUtils.iscall(avg) && SymbolicUtils.operation(avg) isa AvgFunc
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

function _sum_avg(args, state::QuantumState, d::AbstractDict{<:QSym})
    r = ComplexF64(0)
    for a in args
        r += _to_complex(_numeric_average(a, state, d))
    end
    return r
end
function _prod_avg(args, state::QuantumState, d::AbstractDict{<:QSym})
    r = ComplexF64(1)
    for a in args
        r *= _to_complex(_numeric_average(a, state, d))
    end
    return r
end

_numeric_average_vec(op::QField, states::AbstractVector, d::AbstractDict{<:QSym}) = QuantumOpticsBase.expect(QuantumOpticsBase.sparse(to_numeric(op, first(states), d)), states)
_numeric_average_vec(op, states::AbstractVector, d::AbstractDict{<:QSym}) = [_numeric_average(op, ψ, d) for ψ in states]

"""
    expect(op::QField, state[, d::AbstractDict{<:QSym}])

Alias for [`numeric_average`](@ref) on operator expressions. For symbolic
scalar expressions such as `average(op)`, call [`numeric_average`](@ref)
directly.
"""
expect(op::QField, state, d::AbstractDict{<:QSym} = _NO_SUBS) = numeric_average(op, state, d)
