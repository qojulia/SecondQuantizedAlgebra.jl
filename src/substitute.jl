"""
    substitute(expr, d::Dict)

Substitute symbolic parameters and/or operators in `expr` according to dictionary `d`.

Supports two kinds of substitution rules, which can be mixed in a single `Dict`:
- **Symbolic variables** (`Num` / `BasicSymbolic` keys) — substituted in c-number prefactors
  via `Symbolics.substitute`.
- **Operators** ([`QSym`](@ref) keys) — replaced in the operator sequence. If the replacement
  value is a number or symbolic variable, it is folded into the prefactor. If it is a
  single-term [`QAdd`](@ref), its operators and prefactor are spliced in.

# Examples
```julia
h = FockSpace(:f)
@qnumbers a::Destroy(h)
@variables ω
expr = ω * a' * a
substitute(expr, Dict(ω => 2.0))        # 2.0 * a†a
substitute(expr, Dict(a => 0.5 * a))    # 0.5ω * a†a
```
"""
substitute(x::Number, ::Dict) = x

function substitute(op::QSym, d::Dict)
    haskey(d, op) && return d[op]
    return op
end

# Internal: substitute within a single term (prefactor + operator sequence)
function _substitute_term(
        c::CNum, ops::Vector{QSym},
        sym_dict::Dict{SymbolicUtils.BasicSymbolic, SymbolicUtils.BasicSymbolic},
        op_dict::Dict{QSym, Any}
    )
    new_c = if isempty(sym_dict)
        c
    else
        r = Symbolics.substitute(SymbolicUtils.unwrap(real(c)), sym_dict)
        i = Symbolics.substitute(SymbolicUtils.unwrap(imag(c)), sym_dict)
        Complex(Num(r), Num(i))
    end

    new_ops = QSym[]
    extra_c = _CNUM_ONE
    for op in ops
        if haskey(op_dict, op)
            v = op_dict[op]
            if v isa QSym
                push!(new_ops, v)
            elseif v isa QAdd && length(v.arguments) == 1
                (vops, vc) = only(v.arguments)
                extra_c *= vc
                append!(new_ops, vops)
            else
                extra_c *= _to_cnum(v)
            end
        else
            push!(new_ops, op)
        end
    end

    return (new_c * extra_c)::CNum, new_ops
end

# Split a user Dict into symbolic and operator substitution dicts
function _split_sub_dict(d::Dict)
    sym_dict = Dict{SymbolicUtils.BasicSymbolic, SymbolicUtils.BasicSymbolic}()
    op_dict = Dict{QSym, Any}()
    for (k, v) in d
        if k isa QSym
            op_dict[k] = v
        else
            sym_dict[SymbolicUtils.unwrap(k)] = SymbolicUtils.unwrap(v isa Num ? v : Num(v))
        end
    end
    return sym_dict, op_dict
end

function substitute(s::QAdd, d::Dict)
    sym_dict, op_dict = _split_sub_dict(d)
    new_d = QTermDict()
    for (ops, c) in s.arguments
        new_c, new_ops = _substitute_term(c, ops, sym_dict, op_dict)
        _iszero_cnum(new_c) && continue
        _addto!(new_d, new_ops, new_c)
    end
    return QAdd(new_d, s.indices, s.non_equal)
end
