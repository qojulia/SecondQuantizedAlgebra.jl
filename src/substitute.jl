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

function substitute(m::QMul, d::Dict)
    # Split dict into operator substitutions and symbolic substitutions
    sym_dict = Dict{Any, Any}()
    op_dict = Dict{QSym, Any}()
    for (k, v) in d
        if k isa QSym
            op_dict[k] = v
        else
            sym_dict[Symbolics.unwrap(k)] = Symbolics.unwrap(v isa Num ? v : Num(v))
        end
    end

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
            elseif v isa Number || v isa Num
                extra_c *= _to_cnum(v)
            else
                extra_c *= _to_cnum(v)
            end
        else
            push!(new_ops, op)
        end
    end

    result_c = new_c * extra_c
    _iszero_cnum(result_c) && return 0
    if isempty(new_ops)
        return iszero(imag(result_c)) ? real(result_c) : result_c
    end
    if length(new_ops) == 1 && isequal(result_c, _CNUM_ONE)
        return new_ops[1]
    end
    return QMul(result_c, new_ops)
end

function substitute(s::QAdd, d::Dict)
    result_terms = Any[]
    for (ops, c) in s.arguments
        t = substitute(QMul(c, ops), d)
        push!(result_terms, t)
    end
    # Rebuild: sum all substituted terms
    return sum(result_terms)
end
