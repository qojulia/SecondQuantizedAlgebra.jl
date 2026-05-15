"""
    IndexedVariable(name::Symbol, i::Index) -> Num

Create a symbolic single-indexed variable `name(i)` as a `Num`.
"""
function IndexedVariable(name::Symbol, i::Index)
    f = SymbolicUtils.Sym{SymbolicUtils.SymReal}(
        name;
        type = SymbolicUtils.FnType{Tuple{Int}, Real, Nothing},
        shape = UnitRange{Int}[]
    )
    return Num(f(SymbolicUtils.unwrap(i.sym)))
end

struct NotIdentical end

"""
    DoubleIndexedVariable(name::Symbol, i::Index, j::Index; identical::Bool=true) -> Num
"""
function DoubleIndexedVariable(
        name::Symbol, i::Index, j::Index;
        identical::Bool = true
    )
    if !identical && i == j
        return Num(0)
    end
    f = SymbolicUtils.Sym{SymbolicUtils.SymReal}(
        name;
        type = SymbolicUtils.FnType{Tuple{Int, Int}, Real, Nothing},
        shape = UnitRange{Int}[]
    )
    node = f(SymbolicUtils.unwrap(i.sym), SymbolicUtils.unwrap(j.sym))
    if !identical
        node = SymbolicUtils.setmetadata(node, NotIdentical, true)
    end
    return Num(node)
end

change_index(x::Number, ::Index, ::Index) = x
function change_index(x::Num, from::Index, to::Index)
    raw = SymbolicUtils.unwrap(x)
    isym = SymbolicUtils.unwrap(from.sym)
    vars = Symbolics.get_variables(raw)
    any(v -> isequal(v, isym), vars) || return x
    result = Symbolics.substitute(raw, Dict(isym => SymbolicUtils.unwrap(to.sym)))
    return _check_not_identical(result)
end
function change_index(x::CNum, from::Index, to::Index)
    return Complex(change_index(real(x), from, to), change_index(imag(x), from, to))
end

function change_index(op::Destroy, from::Index, to::Index)
    idx = op.index == from ? to : op.index
    return Destroy(op.name, op.space_index, idx)
end
function change_index(op::Create, from::Index, to::Index)
    idx = op.index == from ? to : op.index
    return Create(op.name, op.space_index, idx)
end
function change_index(op::Transition, from::Index, to::Index)
    idx = op.index == from ? to : op.index
    return Transition(op.name, op.i, op.j, op.space_index, idx, op.ground_state, op.n_levels)
end
function change_index(op::Pauli, from::Index, to::Index)
    idx = op.index == from ? to : op.index
    return Pauli(op.name, op.axis, op.space_index, idx)
end
function change_index(op::Spin, from::Index, to::Index)
    idx = op.index == from ? to : op.index
    return Spin(op.name, op.axis, op.space_index, idx)
end
function change_index(op::Position, from::Index, to::Index)
    idx = op.index == from ? to : op.index
    return Position(op.name, op.space_index, idx)
end
function change_index(op::Momentum, from::Index, to::Index)
    idx = op.index == from ? to : op.index
    return Momentum(op.name, op.space_index, idx)
end

function change_index(s::QAdd, from::Index, to::Index)
    new_d = QTermDict()
    for (term, c) in s.arguments
        new_c = change_index(c, from, to)
        new_ops = QSym[change_index(op, from, to) for op in _flatten_chains(term.chains)]
        new_bound = _substitute_bound(term.bound, from, to)
        new_ne = _substitute_ne(term.ne, from, to)
        _accumulate_normalized!(new_d, new_c, new_ops, new_bound, new_ne, get_ordering())
    end
    new_indices = [idx == from ? to : idx for idx in s.indices]
    return _qadd(new_d, new_indices)
end

"""
    get_indices(expr) -> Vector{Index}

Collect all symbolic summation indices appearing in `expr`.
"""
get_indices(::Number) = Index[]
get_indices(x::Num) = get_indices(SymbolicUtils.unwrap(x))
function get_indices(op::QSym)
    return has_index(op.index) ? Index[op.index] : Index[]
end
function get_indices(s::QAdd)
    inds = copy(s.indices)
    for term in keys(s.arguments)
        for idx in term.bound
            idx ∉ inds && push!(inds, idx)
        end
        for op in term.ops
            has_index(op.index) && op.index ∉ inds && push!(inds, op.index)
        end
    end
    return inds
end

function create_index_arrays(indices::Vector{Index}, ranges::Vector{<:AbstractRange})
    length(indices) == 1 && return collect(only(ranges))
    return vec(collect(Iterators.product(ranges...)))
end

function _check_not_identical(result::SymbolicUtils.BasicSymbolic)
    if SymbolicUtils.iscall(result) &&
            SymbolicUtils.hasmetadata(result, NotIdentical) &&
            length(SymbolicUtils.arguments(result)) == 2
        a1, a2 = SymbolicUtils.arguments(result)
        isequal(a1, a2) && return Num(0)
    end
    return Num(result)
end
_check_not_identical(result::Number) = Num(result)

# ============================================================================
#  Index dependence queries
# ============================================================================

function _depends_on_index_term(c::CNum, ops::Vector{QSym}, idx::Index)
    for op in ops
        op.index == idx && return true
    end
    isym = SymbolicUtils.unwrap(idx.sym)
    for part in (real(c), imag(c))
        vars = Symbolics.get_variables(part)
        any(v -> isequal(v, isym), vars) && return true
    end
    return false
end

function _any_depends_on_index(s::QAdd, idx::Index)
    for (term, c) in s.arguments
        _depends_on_index_term(c, term.ops, idx) && return true
    end
    return false
end

# ============================================================================
#  Σ — symbolic sums
# ============================================================================

"""
    Σ(expr, i::Index, non_equal::Vector{Index} = Index[])

Build the symbolic sum `sum_i expr` over index `i`.
"""
function Σ(expr::QAdd, i::Index, non_equal::Vector{Index} = Index[])
    if !_any_depends_on_index(expr, i)
        return expr * i.range
    end

    dependent = QTermDict()
    independent = QTermDict()
    for (term, c) in expr.arguments
        if _depends_on_index_term(c, term.ops, i)
            new_bound = _merge_bound_idx(term.bound, i)
            new_ne = term.ne
            for j in non_equal
                new_ne = _merge_ne_pair(new_ne, i, j)
            end
            _accumulate_normalized!(dependent, c, _flatten_chains(term.chains), new_bound, new_ne, get_ordering())
        else
            scaled_c = _mul_cnum(c, _to_cnum(i.range))
            _addto_key!(independent, _copy_key(term), scaled_c)
        end
    end

    result = _qadd(independent)
    isempty(dependent) || (result = result + _qadd(dependent))
    return result
end

function Σ(expr::QSym, i::Index, non_equal::Vector{Index} = Index[])
    return Σ(_single_qadd(_CNUM_ONE, QSym[expr]), i, non_equal)
end

function Σ(expr::Number, i::Index, non_equal::Vector{Index} = Index[])
    return Σ(_single_qadd(_to_cnum(expr), QSym[]), i, non_equal)
end

function Σ(expr::Union{QField, Number}, i::Index, j::Index, rest::Index...)
    inner = Σ(expr, i)
    return Σ(inner, j, rest...)
end

const ∑ = Σ
