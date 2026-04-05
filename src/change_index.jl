"""
    change_index(expr, from::Index, to::Index)

Substitute index `from` with `to` throughout an expression tree.
Operator indices are replaced directly. Prefactors are substituted
via `Symbolics.substitute` using the `sym` fields of the indices.
"""
change_index(x::Number, ::Index, ::Index) = x
function change_index(x::Num, from::Index, to::Index)
    return Num(
        Symbolics.substitute(
            Symbolics.unwrap(x),
            Dict(Symbolics.unwrap(from.sym) => Symbolics.unwrap(to.sym))
        )
    )
end
function change_index(x::CNum, from::Index, to::Index)
    return Complex(change_index(real(x), from, to), change_index(imag(x), from, to))
end

function change_index(op::Destroy, from::Index, to::Index)
    idx = op.index == from ? to : op.index
    return Destroy(op.name, op.space_index, op.copy_index, idx)
end
function change_index(op::Create, from::Index, to::Index)
    idx = op.index == from ? to : op.index
    return Create(op.name, op.space_index, op.copy_index, idx)
end
function change_index(op::Transition, from::Index, to::Index)
    idx = op.index == from ? to : op.index
    return Transition(op.name, op.i, op.j, op.space_index, op.copy_index, idx)
end
function change_index(op::Pauli, from::Index, to::Index)
    idx = op.index == from ? to : op.index
    return Pauli(op.name, op.axis, op.space_index, op.copy_index, idx)
end
function change_index(op::Spin, from::Index, to::Index)
    idx = op.index == from ? to : op.index
    return Spin(op.name, op.axis, op.space_index, op.copy_index, idx)
end
function change_index(op::Position, from::Index, to::Index)
    idx = op.index == from ? to : op.index
    return Position(op.name, op.space_index, op.copy_index, idx)
end
function change_index(op::Momentum, from::Index, to::Index)
    idx = op.index == from ? to : op.index
    return Momentum(op.name, op.space_index, op.copy_index, idx)
end

function change_index(m::QMul, from::Index, to::Index)
    new_c = change_index(m.arg_c, from, to)
    new_ops = QSym[change_index(op, from, to) for op in m.args_nc]
    return QMul(new_c, new_ops)
end

function change_index(s::QAdd, from::Index, to::Index)
    new_terms = QMul[change_index(t, from, to) for t in s.arguments]
    new_indices = [idx == from ? to : idx for idx in s.indices]
    new_ne = Tuple{Index, Index}[
        (a == from ? to : a, b == from ? to : b)
            for (a, b) in s.non_equal
    ]
    return QAdd(new_terms, new_indices, new_ne)
end

"""
    get_indices(expr) -> Vector{Index}

Collect all non-NO_INDEX indices in an expression.
"""
get_indices(::Number) = Index[]
get_indices(::Num) = Index[]
function get_indices(op::QSym)
    return has_index(op.index) ? Index[op.index] : Index[]
end
function get_indices(m::QMul)
    inds = Index[]
    for op in m.args_nc
        has_index(op.index) && op.index ∉ inds && push!(inds, op.index)
    end
    return inds
end
function get_indices(s::QAdd)
    inds = copy(s.indices)
    for t in s.arguments
        for idx in get_indices(t)
            idx ∉ inds && push!(inds, idx)
        end
    end
    return inds
end
