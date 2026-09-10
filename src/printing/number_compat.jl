# Generic fallback for symbolic numeric wrappers. `Complex{Num}` keeps the specialized
# renderer in `printing.jl`; other `Number` wrappers can be printed directly without SQA
# depending on their concrete wrapper type.
show_display(io::IO, c::Number) = print(io, c)

function needs_pf_parens(c::Number)
    u = SymbolicUtils.unwrap(c)
    u isa SymbolicUtils.BasicSymbolic || return false
    SymbolicUtils.iscall(u) || return false
    op = SymbolicUtils.operation(u)
    return op === (+) || op === (/)
end
