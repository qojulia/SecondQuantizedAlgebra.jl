"""
    QField

Supertype of every operator expression.

Operators ([`QSym`](@ref)), their polynomials ([`QAdd`](@ref)) and expressions with formal
functions ([`QExpr`](@ref)) are all `QField`s. They support `+`, `-`, `*`, `/` by numbers,
integer powers, `adjoint`, and comparison with `==` and `isequal`. A method written for
`QField` accepts any of them.

# Examples
```jldoctest
julia> a = Destroy(FockSpace(:cavity), :a);

julia> (a, a' * a, cos(a + a')) isa NTuple{3, SecondQuantizedAlgebra.QField}
true
```

See also [`QSym`](@ref), [`QAdd`](@ref), [`QExpr`](@ref).
"""
abstract type QField end

"""
    QSym <: QField

Supertype of single operators.

Every operator, whatever its kind, is an [`Op`](@ref), the only subtype; its kind is read
with [`optype`](@ref).

See also [`QField`](@ref), [`Op`](@ref).
"""
abstract type QSym <: QField end

Base.zero(::T) where {T <: QField} = zero(T)
Base.zero(::Type{<:QField}) = 0

"""
    SiteCmp

Result of `site_compare`, the site relation of two operators in the canonical sort.

`Less` and `Greater` order operators on distinct sites. `Equal` marks a provably shared
site, where the reduce and commute passes may act. `Undetermined` marks two free indices
that may or may not coincide, where no rule may act. Both of the latter leave the pair in
place; they are separate values so that a comparison cannot return one where the other
is meant.
"""
@enum SiteCmp::UInt8 Less Equal Undetermined Greater

"""
    ReduceKind

Outcome tag for `reduce_pair`:
- `NoReduction`: pair does not reduce; `op` and `factor` slots are ignored.
- `ScalarReduction`: pair contracts to the scalar `factor`; both ops disappear.
- `OpReduction`: pair reduces to `op * factor`; `op` replaces the first input.
"""
@enum ReduceKind::UInt8 NoReduction ScalarReduction OpReduction

function position_of(pred::F, v::Union{AbstractVector, Tuple}) where {F}
    for i in eachindex(v)
        pred(v[i]) && return i
    end
    return 0
end
