"""
    @qnumbers

Convenience macro for the construction of operators.

Examples
========
```julia
h = FockSpace(:fock)
@qnumbers a::Destroy(h)

h = FockSpace(:one) ⊗ FockSpace(:two)
@qnumbers a::Destroy(h, 1) b::Destroy(h, 2)
```
"""
macro qnumbers(qs...)
    ex = Expr(:block)
    qnames = []
    for q in qs
        @assert q isa Expr && q.head == :(::)
        name = q.args[1]
        @assert name isa Symbol
        push!(qnames, name)
        f = q.args[2]
        @assert f isa Expr && f.head == :call
        op_type = f.args[1]
        op_args = f.args[2:end]
        name_quoted = Expr(:quote, name)
        # Insert name as second argument: Op(hilbert, name, extra_args...)
        construction = Expr(:call, esc(op_type), esc(op_args[1]), name_quoted, map(esc, op_args[2:end])...)
        push!(ex.args, :($(esc(name)) = $(construction)))
    end
    push!(ex.args, Expr(:tuple, map(esc, qnames)...))
    return ex
end
