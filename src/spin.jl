"""
    SpinSpace <: HilbertSpace

Hilbert space for collective spin operators (Spin-S).
"""
struct SpinSpace <: HilbertSpace
    name::Symbol
    spin::Rational{Int}
    function SpinSpace(name::Symbol, spin::Rational{Int})
        spin > 0 || throw(ArgumentError("Spin must be positive, got $spin"))
        denominator(spin) in (1, 2) || throw(ArgumentError("Spin must be integer or half-integer, got $spin"))
        new(name, spin)
    end
    SpinSpace(name::Symbol, spin::Integer) = SpinSpace(name, spin // 1)
end
Base.:(==)(a::SpinSpace, b::SpinSpace) = a.name == b.name && a.spin == b.spin
Base.hash(a::SpinSpace, h::UInt) = hash(:SpinSpace, hash(a.name, hash(a.spin, h)))

"""
    Spin <: QSym

Angular momentum operator (Sx, Sy, Sz) on a [`SpinSpace`](@ref).
Axis: 1=x, 2=y, 3=z.
"""
struct Spin <: QSym
    name::Symbol
    axis::Int
    space_index::Int
    function Spin(name::Symbol, axis::Int, space_index::Int)
        1 <= axis <= 3 || throw(ArgumentError("Spin axis must be 1, 2, or 3, got $axis"))
        new(name, axis, space_index)
    end
end

# Construction from Hilbert spaces
Spin(h::SpinSpace, name::Symbol, axis::Int) = Spin(name, axis, 1)
function Spin(h::ProductSpace, name::Symbol, axis::Int, idx::Int)
    1 <= idx <= length(h.spaces) || throw(ArgumentError("Index $idx out of range"))
    h.spaces[idx] isa SpinSpace || throw(ArgumentError("Space at index $idx is not a SpinSpace"))
    return Spin(name, axis, idx)
end

# Adjoint — Hermitian
Base.adjoint(op::Spin) = op

# Equality
Base.isequal(a::Spin, b::Spin) = a.name == b.name && a.axis == b.axis && a.space_index == b.space_index
Base.:(==)(a::Spin, b::Spin) = isequal(a, b)

# Hashing
Base.hash(a::Spin, h::UInt) = hash(:Spin, hash(a.name, hash(a.axis, hash(a.space_index, h))))

# Ladder (not applicable)
ladder(::Spin) = 0
