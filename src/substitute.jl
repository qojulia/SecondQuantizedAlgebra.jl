"""
    substitute(expr, d::Dict)

Substitute symbolic parameters and operators in `expr` according to `d`.

Keys in `d` can be:
- Symbolic variables (`Num` / `BasicSymbolic`) — substituted in c-number prefactors
- Operators (`QSym`) — replaced in the operator sequence; if the replacement is
  a number or symbolic variable, it folds into the prefactor.
"""
substitute(x::Number, ::Dict) = x

function substitute(op::QSym, d::Dict)
    haskey(d, op) && return d[op]
    return op
end

# Internal: always returns QMul (type-stable for QAdd aggregation)
function _substitute_qmul(m::QMul, sym_dict::Dict, op_dict::Dict{QSym, Any})
    # Substitute in prefactor
    new_c = if isempty(sym_dict)
        m.arg_c
    else
        r = Symbolics.substitute(Symbolics.unwrap(real(m.arg_c)), sym_dict)
        i = Symbolics.substitute(Symbolics.unwrap(imag(m.arg_c)), sym_dict)
        Complex(Num(r), Num(i))
    end

    # Substitute operators — replacements that are numbers fold into prefactor
    new_ops = QSym[]
    extra_c = _CNUM_ONE
    for op in m.args_nc
        if haskey(op_dict, op)
            v = op_dict[op]
            if v isa QSym
                push!(new_ops, v)
            elseif v isa QMul
                extra_c *= v.arg_c
                append!(new_ops, v.args_nc)
            else
                extra_c *= _to_cnum(v)
            end
        else
            push!(new_ops, op)
        end
    end

    return QMul(new_c * extra_c, new_ops)
end

# Split a user Dict into symbolic and operator substitution dicts
function _split_sub_dict(d::Dict)
    sym_dict = Dict{Any, Any}()
    op_dict = Dict{QSym, Any}()
    for (k, v) in d
        if k isa QSym
            op_dict[k] = v
        else
            sym_dict[Symbolics.unwrap(k)] = Symbolics.unwrap(v isa Num ? v : Num(v))
        end
    end
    return sym_dict, op_dict
end

function substitute(m::QMul, d::Dict)
    sym_dict, op_dict = _split_sub_dict(d)
    result = _substitute_qmul(m, sym_dict, op_dict)
    # Simplify output: unwrap trivial cases
    iszero(result) && return 0
    if isempty(result.args_nc)
        c = result.arg_c
        return iszero(imag(c)) ? real(c) : c
    end
    if length(result.args_nc) == 1 && isequal(result.arg_c, _CNUM_ONE)
        return result.args_nc[1]
    end
    return result
end

function substitute(s::QAdd, d::Dict)
    sym_dict, op_dict = _split_sub_dict(d)
    new_terms = QMul[_substitute_qmul(QMul(c, ops), sym_dict, op_dict) for (ops, c) in s.arguments]
    return QAdd(new_terms)
end
