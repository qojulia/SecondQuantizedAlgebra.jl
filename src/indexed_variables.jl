"""
    IndexedVariable(name::Symbol, i::Index) -> Num

Create a symbolic indexed variable `name(i)`.
Returns a `Num` from Symbolics.jl.
"""
function IndexedVariable(name::Symbol, i::Index)
    f = SymbolicUtils.Sym{SymbolicUtils.SymReal}(name;
        type = SymbolicUtils.FnType{Tuple{Int}, Real, Nothing},
        shape = UnitRange{Int}[])
    return Num(f(Symbolics.unwrap(i.sym)))
end

"""
    DoubleIndexedVariable(name::Symbol, i::Index, j::Index; identical::Bool=true) -> Num

Create a symbolic double-indexed variable `name(i, j)`.
If `identical=false`, returns `Num(0)` when `i == j`.
"""
function DoubleIndexedVariable(name::Symbol, i::Index, j::Index;
    identical::Bool = true)
    if !identical && i == j
        return Num(0)
    end
    f = SymbolicUtils.Sym{SymbolicUtils.SymReal}(name;
        type = SymbolicUtils.FnType{Tuple{Int, Int}, Real, Nothing},
        shape = UnitRange{Int}[])
    return Num(f(Symbolics.unwrap(i.sym), Symbolics.unwrap(j.sym)))
end
