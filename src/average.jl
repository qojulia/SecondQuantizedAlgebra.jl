"""
    AvgSym

Marker type for the `symtype` of averaged operator expressions.
Average nodes are `BasicSymbolic{SymReal}` `Term` nodes with `symtype(x) === AvgSym`.
"""
struct AvgSym <: Number end

"""Metadata key for summation indices on averaged expressions."""
struct SumIndices end

"""Metadata key for non-equal index pairs on averaged expressions."""
struct SumNonEqual end

"""
    AvgFunc

Singleton callable struct used as the `operation` of average `Term` nodes.
Using a custom struct (instead of a SymbolicUtils `Sym`) lets us define
`SymbolicUtils.show_call(::IO, ::AvgFunc, ...)` for `⟨op⟩` display
without type piracy.
"""
struct AvgFunc end

"""Singleton instance of [`AvgFunc`](@ref)."""
const sym_average = AvgFunc()

Base.nameof(::AvgFunc) = :avg
Base.show(io::IO, ::AvgFunc) = print(io, "avg")

# ⟨op⟩ display for average Term nodes
function SymbolicUtils.show_call(
        io::IO, ::AvgFunc, x::SymbolicUtils.BasicSymbolic; kw...
    )
    args = SymbolicUtils.arguments(x)
    print(io, "⟨")
    for (i, arg) in enumerate(args)
        i > 1 && print(io, ", ")
        print(io, arg)
    end
    return print(io, "⟩")
end

"""
    is_average(x) -> Bool

Check whether `x` is a symbolic average `Term` node.
"""
is_average(::QField) = false
is_average(::Number) = false
function is_average(x::SymbolicUtils.BasicSymbolic)
    return SymbolicUtils.iscall(x) && SymbolicUtils.symtype(x) === AvgSym
end
is_average(x) = false

# --- Internal constructor ---

function _average(op::QField)
    return SymbolicUtils.Term{SymbolicUtils.SymReal}(sym_average, Any[op]; type = AvgSym)
end

# --- Public API ---

"""
    average(expr)

Compute the symbolic average of a quantum operator expression.
Returns a `BasicSymbolic{SymReal}` scalar with `symtype = AvgSym` that participates
in Symbolics arithmetic.

Linearity: distributes over `QAdd` and pulls c-number prefactors out of `QMul`.
"""
function average end

average(op::QSym) = _average(op)
average(x::Number) = x
average(x::SymbolicUtils.BasicSymbolic) = x

function average(op::QMul)
    isempty(op.args_nc) && return op.arg_c
    # Normalize: single-operator QMul → bare QSym so that _average always wraps
    # the same type regardless of whether the caller wrote average(σ) or average(1*σ).
    # For multi-operator QMul, strip the prefactor to get a unit-prefactor QMul.
    inner = length(op.args_nc) == 1 ? only(op.args_nc) : QMul(op.args_nc)
    # Fast path: unit prefactor, avoid allocating a new QMul
    isone(op.arg_c) && return _average(inner)
    # Pull out prefactor: _average wraps the pure operator part.
    avg = _average(inner)
    # CNum = Complex{Num}. Multiplying Complex{Num} * BasicSymbolic may not
    # dispatch correctly, so split into real/imag parts which are Num values
    # and multiply through Symbolics arithmetic.
    r, i = real(op.arg_c), imag(op.arg_c)
    iszero(i) && return r * avg
    iszero(r) && return im * i * avg
    return (r + im * i) * avg
end

function average(op::QAdd)
    result = mapreduce(t -> average(t), +, terms(op))
    if !isempty(op.indices)
        result = SymbolicUtils.setmetadata(result, SumIndices, op.indices)
        result = SymbolicUtils.setmetadata(result, SumNonEqual, op.non_equal)
    end
    return result
end

# --- Internal helpers ---

"""Restore summation metadata from a SymbolicUtils node onto the result QField expression."""
function _restore_sum_metadata(result, x::SymbolicUtils.BasicSymbolic)
    if SymbolicUtils.hasmetadata(x, SumIndices)
        indices = SymbolicUtils.getmetadata(x, SumIndices)
        non_equal = SymbolicUtils.getmetadata(x, SumNonEqual)
        if result isa QAdd
            return QAdd(result.arguments, indices, non_equal)
        elseif result isa QMul
            return QAdd(QMul[result], indices, non_equal)
        elseif result isa QSym
            return QAdd(QMul[_to_qmul(result)], indices, non_equal)
        end
    end
    return result
end

# --- undo_average ---

"""
    undo_average(expr)

Recursively strip `avg(...)` from a symbolic expression, recovering operator expressions.
Summation metadata is restored to the resulting `QAdd`.

Note: this function is inherently non-inferrable because it walks SymbolicUtils expression
trees via `operation(x)` (type `Any`) and `f(args...)` (dynamic dispatch with splatting).
This is acceptable — `undo_average` is not a hot path.

For Add/Mul nodes, `f(args...)` calls `+`/`*` on QField values which creates QAdd/QMul
via SQA's algebra. For non-standard operations (Pow, Div, etc.), `f(args...)` reconstructs
through Julia's generic dispatch — this is correct but may not preserve SymbolicUtils node
structure exactly. In practice, QC only produces Add/Mul trees from averaging.
"""
function undo_average(x::SymbolicUtils.BasicSymbolic)
    if SymbolicUtils.iscall(x)
        f = SymbolicUtils.operation(x)
        if f isa AvgFunc
            arg = SymbolicUtils.arguments(x)[1]
            # SymbolicUtils wraps non-symbolic args as Const{SymReal}(val);
            # extract the original QField. Const has no public extraction API,
            # so we access .val directly (SymbolicUtils internal).
            result = SymbolicUtils.isconst(arg) ? arg.val : arg
            return _restore_sum_metadata(result, x)
        elseif f === (+) || f === (*)
            # Additive/multiplicative nodes may contain averages as children.
            # Recurse into args and reconstruct with QField algebra.
            args = map(undo_average, SymbolicUtils.arguments(x))
            result = f(args...)
            return _restore_sum_metadata(result, x)
        else
            # Any other call (complex, ^, conj, etc.) is a pure c-number node.
            # Convert to CNum to stay in the QField type domain.
            return _to_cnum(x)
        end
    else
        # Non-call BasicSymbolic (bare symbol or Const) — convert to CNum
        return _to_cnum(x)
    end
end

undo_average(x::Number) = x
function undo_average(x::Num)
    inner = undo_average(Symbolics.unwrap(x))
    # Num only wraps BasicSymbolic; if undo_average recovered a QField/Number, return as-is
    return inner isa SymbolicUtils.BasicSymbolic ? Num(inner) : inner
end
undo_average(x::QField) = x
undo_average(x) = x  # catch-all for types not covered above (e.g. Symbol)

# Symbolics.Equation fields are BasicSymbolic{SymReal}, so we can only reconstruct an
# Equation when the un-averaged results are still BasicSymbolic. Return a Pair otherwise.
function undo_average(eq::Symbolics.Equation)
    lhs = undo_average(eq.lhs)
    rhs = undo_average(eq.rhs)
    if lhs isa SymbolicUtils.BasicSymbolic && rhs isa SymbolicUtils.BasicSymbolic
        return Symbolics.Equation(lhs, rhs)
    end
    return lhs => rhs
end

# --- Metadata helpers ---

"""
    has_sum_metadata(x) -> Bool

Check whether a symbolic expression carries summation index metadata.
"""
has_sum_metadata(::Any) = false
function has_sum_metadata(x::SymbolicUtils.BasicSymbolic)
    return SymbolicUtils.hasmetadata(x, SumIndices)
end

"""
    get_sum_indices(x::BasicSymbolic) -> Vector{Index}

Retrieve summation indices from a symbolic expression's metadata.
"""
function get_sum_indices(x::SymbolicUtils.BasicSymbolic)
    return SymbolicUtils.getmetadata(x, SumIndices)
end

"""
    get_sum_non_equal(x::BasicSymbolic) -> Vector{Tuple{Index,Index}}

Retrieve non-equal index pairs from a symbolic expression's metadata.
"""
function get_sum_non_equal(x::SymbolicUtils.BasicSymbolic)
    return SymbolicUtils.getmetadata(x, SumNonEqual)
end

# --- get_indices for BasicSymbolic (averages, sums of averages) ---

function get_indices(x::SymbolicUtils.BasicSymbolic)
    if SymbolicUtils.isconst(x)
        return get_indices(x.val)
    end
    if SymbolicUtils.iscall(x)
        f = SymbolicUtils.operation(x)
        if f isa AvgFunc
            arg = SymbolicUtils.arguments(x)[1]
            inner = SymbolicUtils.isconst(arg) ? arg.val : arg
            return get_indices(inner)
        end
        inds = Index[]
        for arg in SymbolicUtils.arguments(x)
            for idx in get_indices(arg)
                idx ∉ inds && push!(inds, idx)
            end
        end
        return inds
    end
    return Index[]
end

# --- acts_on ---

"""
    acts_on(expr) -> Vector{Int}

Return sorted unique Hilbert space indices that an operator or averaged expression acts on.
Not intended for hot-path use — allocates a fresh `Vector{Int}` on each call.
"""
acts_on(op::QSym) = Int[op.space_index]

function acts_on(op::QMul)
    aon = Int[x.space_index for x in op.args_nc]
    unique!(sort!(aon))
    return aon
end

function acts_on(op::QAdd)
    aon = Int[]
    for t in terms(op)
        append!(aon, acts_on(t))
    end
    unique!(sort!(aon))
    return aon
end

function acts_on(s::SymbolicUtils.BasicSymbolic)
    # Const nodes wrap non-symbolic values (e.g., QField inside avg args).
    # Const has no public extraction API, so we access .val directly (SymbolicUtils internal).
    if SymbolicUtils.isconst(s)
        return acts_on(s.val)
    end
    if SymbolicUtils.iscall(s)
        f = SymbolicUtils.operation(s)
        if f isa AvgFunc
            return acts_on(SymbolicUtils.arguments(s)[1])
        else
            aon = Int[]
            for arg in SymbolicUtils.arguments(s)
                append!(aon, acts_on(arg))
            end
            unique!(sort!(aon))
            return aon
        end
    else
        return Int[]
    end
end

acts_on(::Number) = Int[]
