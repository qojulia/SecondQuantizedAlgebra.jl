"""
    PhaseSpace <: HilbertSpace

[`HilbertSpace`](@ref) defining the Phase space for the position and momentum operators.
See also: [`Position`](@ref), [`Momentum`](@ref)
"""
struct PhaseSpace{S} <: ConcreteHilbertSpace
    name::S
end
Base.:(==)(h1::T, h2::T) where {T<:PhaseSpace} = (h1.name==h2.name)

"""
    Position <: QSym

Position operator on a [`PhaseSpace`](@ref).
"""
struct Position{H<:HilbertSpace,S,A,M} <: QSym
    hilbert::H
    name::S
    aon::A
    metadata::M
    function Position{H,S,A,M}(hilbert::H, name::S, aon::A, metadata::M) where {H,S,A,M}
        @assert has_hilbert(PhaseSpace, hilbert, aon)
        new(hilbert, name, aon, metadata)
    end
end

"""
    Momentum <: QSym

Momentum operator on a [`PhaseSpace`](@ref).
"""
struct Momentum{H<:HilbertSpace,S,A,M} <: QSym
    hilbert::H
    name::S
    aon::A
    metadata::M
    function Momentum{H,S,A,M}(hilbert::H, name::S, aon::A, metadata::M) where {H,S,A,M}
        @assert has_hilbert(PhaseSpace, hilbert, aon)
        new(hilbert, name, aon, metadata)
    end
end

for T in (:Position, :Momentum)
    @eval Base.isequal(a::$T, b::$T) =
        isequal(a.hilbert, b.hilbert) && isequal(a.name, b.name) && isequal(a.aon, b.aon)
end

for f in [:Position, :Momentum]
    @eval $(f)(hilbert::H, name::S, aon::A; metadata::M=NO_METADATA) where {H,S,A,M} = $(f){
        H,S,A,M
    }(
        hilbert, name, aon, metadata
    )
    @eval $(f)(hilbert::PhaseSpace, name; metadata=NO_METADATA) = $(f)(
        hilbert, name, 1; metadata
    )
    @eval function $(f)(
        hilbert::H, name::S, aon::A; metadata::M=NO_METADATA
    ) where {H<:ProductSpace,S,A<:Int,M}
        if hilbert.spaces[aon] isa ClusterSpace
            hilbert.spaces[get_i(aon)].op_name[] = name
            op = $(f)(hilbert.spaces[aon].original_space, name; metadata)
            return _cluster(hilbert, op, aon) #todo OK?
        else
            return $(f){H,S,A,M}(hilbert, name, aon, metadata)
        end
    end
    @eval function $(f)(hilbert::ProductSpace, name; metadata=NO_METADATA)
        i = findall(
            x->isa(x, PhaseSpace) || isa(x, ClusterSpace{<:PhaseSpace}), hilbert.spaces
        )
        if length(i)==1
            return $(f)(hilbert, name, i[1]; metadata)
        else
            isempty(i) &&
                error("Can only create $($(f)) on PhaseSpace! Not included in $(hilbert)")
            length(i)>1 && error(
                "More than one PhaseSpace in $(hilbert)! Specify on which Hilbert space $($(f)) should be created with $($(f))(hilbert,name,i)!",
            )
        end
    end
    @eval function Base.hash(op::T, h::UInt) where {T<:($(f))}
        hash($(f), hash(op.hilbert, hash(op.name, hash(op.aon, h))))
    end
end

Base.adjoint(op::Position) = op
Base.adjoint(op::Momentum) = op

# Commutation relation in simplification
function *(p::Momentum, x::Position)
    check_hilbert(p, x)
    aon_p = acts_on(p)
    aon_x = acts_on(x)
    if aon_p == aon_x
        return x*p - 1im
    elseif aon_p < aon_x
        return QMul(1, [p, x])
    else
        return QMul(1, [x, p])
    end
end
# ismergeable(::Momentum,::Position) = true

# TODO: test if faster; delete if and elseif in *-function above? See also fock.jl
function ismergeable(p::Momentum, x::Position)
    aon_p = acts_on(p)
    aon_x = acts_on(x)
    return aon_p == aon_x
end
