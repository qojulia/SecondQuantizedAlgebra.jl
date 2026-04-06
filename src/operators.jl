"""
    fundamental_operators(h::HilbertSpace; names=nothing) -> Vector{QSym}

Return all fundamental operators for a Hilbert space.
For example, a `FockSpace` has one fundamental operator (`Destroy`),
while an `NLevelSpace` with `n` levels has `n(n+1)/2 - 1` transitions
(excluding the ground state projector).
"""
function fundamental_operators(h::FockSpace, si::Int = 1; names = nothing)
    name = names === nothing ? :a : names[si]
    return QSym[Destroy(name, si)]
end

function fundamental_operators(h::NLevelSpace, si::Int = 1; names = nothing)
    name = names === nothing ? :σ : names[si]
    ops = Transition[]
    for i in 1:(h.n)
        for j in i:(h.n)
            (i == j) && i == h.ground_state && continue
            push!(ops, Transition(name, i, j, si))
        end
    end
    return ops
end

function fundamental_operators(h::PauliSpace, si::Int = 1; names = nothing)
    name = names === nothing ? :σ : names[si]
    return QSym[Pauli(name, 1, si), Pauli(name, 2, si), Pauli(name, 3, si)]
end

function fundamental_operators(h::SpinSpace, si::Int = 1; names = nothing)
    name = names === nothing ? :S : names[si]
    return QSym[Spin(name, 1, si), Spin(name, 2, si), Spin(name, 3, si)]
end

function fundamental_operators(h::PhaseSpace, si::Int = 1; names = nothing)
    name_pair = names === nothing ? (:x, :p) : names[si]
    return QSym[Position(name_pair[1], si), Momentum(name_pair[2], si)]
end

function fundamental_operators(h::ProductSpace; names = nothing)
    ops = QSym[]
    for (i, space) in enumerate(h.spaces)
        space_ops = fundamental_operators(_unwrap_space(space), i; names = names)
        append!(ops, space_ops)
    end
    return ops
end

"""
    unique_ops(ops) -> Vector

Return only unique operators from `ops`, considering that `op` and `op'`
represent the same degree of freedom.
"""
function unique_ops(ops::AbstractVector)
    ops_ = copy(ops)
    unique_ops!(ops_)
    return ops_
end

"""
    unique_ops!(ops) -> Vector

In-place version of [`unique_ops`](@ref).
"""
function unique_ops!(ops::AbstractVector)
    hashes = map(hash, ops)
    hashes_adj = map(x -> hash(adjoint(x)), ops)
    seen = UInt[]
    i = 1
    while i <= length(ops)
        if hashes[i] in seen || hashes_adj[i] in seen
            deleteat!(ops, i)
            deleteat!(hashes, i)
            deleteat!(hashes_adj, i)
        else
            push!(seen, hashes[i])
            hashes[i] == hashes_adj[i] || push!(seen, hashes_adj[i])
            i += 1
        end
    end
    return ops
end

"""
    find_operators(h::HilbertSpace, order::Int; names=nothing) -> Vector

Find all operator products that fully define a system up to the given `order`.
Returns unique operators and their products (up to `order` factors),
excluding zero terms and duplicates under adjoint.
"""
function find_operators(h::HilbertSpace, order::Int; names = nothing)
    # Auto-generate names when ProductSpace has duplicate space types
    if names === nothing && h isa ProductSpace
        space_types = typeof.(_unwrap_space.(h.spaces))
        if length(unique(space_types)) != length(space_types)
            alph = 'a':'z'
            names = Symbol.(alph[1:length(h.spaces)])
        end
    end
    fund_ops = fundamental_operators(h; names = names)
    fund_ops = unique(vcat(fund_ops, adjoint.(fund_ops)))

    all_ops = QField[]
    for i in 1:order
        for c in with_replacement_combinations(fund_ops, i)
            c_ = prod(reverse(c))
            iszero(c_) || push!(all_ops, c_)
        end
    end

    filter!(x -> !(x isa QAdd), all_ops)
    unique_ops!(all_ops)
    return all_ops
end


# --- Conjugation helpers ---

"""
    _conj(v)

Recursively conjugate a symbolic expression. For `BasicSymbolic` trees,
applies `conj` to leaves and reconstructs via `maketerm`.
"""
function _conj(v::SymbolicUtils.BasicSymbolic)
    if SymbolicUtils.iscall(v)
        f = SymbolicUtils.operation(v)
        args = map(_conj, SymbolicUtils.arguments(v))
        return TermInterface.maketerm(typeof(v), f, args, TermInterface.metadata(v))
    else
        return conj(v)
    end
end
_conj(x::Number) = conj(x)
_conj(x::QField) = adjoint(x)

"""
    _inconj(v)

Conjugate an expression containing averages, preserving operator ordering inside
averages (uses `adjoint` on the inner operator rather than naive `conj`).
"""
function _inconj(v::SymbolicUtils.BasicSymbolic)
    if SymbolicUtils.iscall(v)
        f = SymbolicUtils.operation(v)
        if f isa AvgFunc
            arg = SymbolicUtils.arguments(v)[1]
            inner = SymbolicUtils.isconst(arg) ? arg.val : arg
            return _average(adjoint(inner))
        elseif f === conj
            # conj(avg(...)) → recurse into the argument
            return _inconj(SymbolicUtils.arguments(v)[1])
        else
            args = map(_inconj, SymbolicUtils.arguments(v))
            return TermInterface.maketerm(typeof(v), f, args, TermInterface.metadata(v))
        end
    else
        return conj(v)
    end
end
_inconj(x::Number) = conj(x)
_inconj(x::QField) = adjoint(x)

"""
    _adjoint(x)

Adjoint that works on both operators (uses `adjoint`) and symbolic numbers
(uses `_conj`).
"""
_adjoint(op::QField) = adjoint(op)
_adjoint(s::SymbolicUtils.BasicSymbolic) = _conj(s)
_adjoint(x::Number) = adjoint(x)
