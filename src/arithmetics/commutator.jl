"""
    commutator(a, b) -> QAdd

Compute the commutator ``[a, b] = a*b - b*a``.

For indexed expressions, the diagonal substitutions performed by `*` (via
`_accumulate_with_diag!`) automatically produce the correct partial-collapse
contributions when summation indices share a Hilbert subspace with the other
factor — no separate index-collapse is needed in this function.

# Examples
```julia
h = FockSpace(:f)
@qnumbers a::Destroy(h)
commutator(a, a')    # [a, a†] = 1
commutator(a', a)    # [a†, a] = -1
```

See also [`normal_order`](@ref).
"""
function commutator end

const _ZERO_QADD = QAdd(QTermDict(), Index[])
_zero_qadd() = _ZERO_QADD

# Scalars commute with everything
commutator(::Number, ::Number) = _zero_qadd()
commutator(::Number, ::QField) = _zero_qadd()
commutator(::QField, ::Number) = _zero_qadd()

# QSym, QSym: short-circuit when on different sites or identical
function commutator(a::QSym, b::QSym)
    _same_site(a, b) || return _zero_qadd()
    isequal(a, b) && return _zero_qadd()
    return a * b - b * a
end

# QAdd ↔ QSym defer to `*`, which performs the partial-collapse diagonal
# substitutions via `_accumulate_with_diag!`.
commutator(a::QAdd, b::QSym) = a * b - b * a
commutator(a::QSym, b::QAdd) = a * b - b * a

function commutator(a::QAdd, b::QAdd)
    isequal(a, b) && return _zero_qadd()
    return a * b - b * a
end

"""
    anticommutator(a, b)

Compute the anticommutator ``\\{a, b\\} = a*b + b*a``.

Returns a [`QAdd`](@ref) when either argument is a [`QField`](@ref); for two
plain numbers returns `2*a*b`.

# Examples
```julia
h = FockSpace(:f)
@qnumbers a::Destroy(h)
anticommutator(a, a')   # {a, a†} = 1 + 2 a†a
```

See also [`commutator`](@ref).
"""
anticommutator(a, b) = a * b + b * a
