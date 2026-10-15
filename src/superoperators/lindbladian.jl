"""
    Lindbladian{O}

A Lindblad generator that keeps its Hamiltonian and decay channels.

The generator is `L = -i[H, ·] + Σ γ 𝒟[J, K]`. `O` is `QAdd` for a polynomial model and `QExpr` when the Hamiltonian or a jump contains a
formal operator function. Build a generator with [`hamiltonian_action`](@ref),
[`dissipator`](@ref) and [`liouvillian`](@ref), and combine generators with `+`, scalar `*`
and [`Σ`](@ref). `L(ρ)` applies the generator, and `transpose(L)(A)` is the Heisenberg drift
of the observable `A`.

A `Lindbladian` is a form, not a claim of physicality: `H` may be non-Hermitian and the
rates complex or negative. It always preserves the trace. It preserves Hermiticity when `H`
and the rate matrix are Hermitian, and it generates a completely positive evolution when
the rate matrix is also positive semidefinite.

Equality is structural. Equal spellings of one channel compare equal, but a cross term is
not expanded, so two generators that act alike are compared with
`SuperOperator(L₁) == SuperOperator(L₂)`.

# Examples
```jldoctest
julia> h = FockSpace(:c); a = Destroy(h, :a);

julia> L = hamiltonian_action(a' * a) + 2 * dissipator(a)
-i[a' * a, ·] + 2 * 𝒟[a]

julia> transpose(L)(a)
(-1 + -1im) * a
```

See also [`SuperOperator`](@ref), [`to_numeric_master`](@ref).
"""
struct Lindbladian{O <: Union{QAdd, QExpr}}
    hamiltonian::O
    channels::Dict{Channel{O}, CNum}
end

"""
    hamiltonian(L::Lindbladian)
    channels(L::Lindbladian)

The Hamiltonian of a Lindbladian, and a vector of its `Channel => rate` pairs in
[`channel_order_key`](@ref) order.

# Examples
```jldoctest
julia> using SecondQuantizedAlgebra: hamiltonian, channels, jump

julia> h = FockSpace(:c); a = Destroy(h, :a);

julia> L = hamiltonian_action(a' * a) + 2 * dissipator(a);

julia> hamiltonian(L)
a' * a

julia> ch, rate = first(channels(L)); jump(ch), rate
(a, 2)
```

See also [`Channel`](@ref), [`jump`](@ref).
"""
hamiltonian(L::Lindbladian) = L.hamiltonian
function channels(L::Lindbladian)
    entries = collect(L.channels)
    sort!(entries; by = entry -> channel_order_key(entry.first))
    return entries
end
@doc (@doc hamiltonian) channels

"""
    LindbladianBuilder{O}()

An accumulator that builds a `Lindbladian{O}` from many terms in linear time.

Add generators with [`accumulate!`](@ref) and finish with [`build`](@ref). `sum` over an
array of Lindbladians uses it already.

# Examples
```jldoctest
julia> using SecondQuantizedAlgebra: LindbladianBuilder, QAdd, accumulate!, build

julia> h = FockSpace(:c); a = Destroy(h, :a);

julia> b = LindbladianBuilder{QAdd}();

julia> accumulate!(b, dissipator(a)); accumulate!(b, dissipator(a));

julia> build(b)
2 * 𝒟[a]
```

See also [`SuperBuilder`](@ref).
"""
struct LindbladianBuilder{O <: Union{QAdd, QExpr}, HB}
    hamiltonian::HB
    channels::Dict{Channel{O}, CNum}
end
LindbladianBuilder{O}() where {O <: Union{QAdd, QExpr}} =
    LindbladianBuilder{O, typeof(hamiltonian_accumulator(O))}(
    hamiltonian_accumulator(O), Dict{Channel{O}, CNum}(),
)

function accumulate!(b::LindbladianBuilder{O}, L::Lindbladian{O}, scale::CNum = CNUM_ONE) where {O}
    iszero_cnum(scale) && return b
    H = L.hamiltonian
    add_hamiltonian!(b.hamiltonian, isequal(scale, CNUM_ONE) ? H : scale_jump(H, scale))
    for (ch, γ) in L.channels
        add_rate!(b.channels, ch, mul_cnum(γ, scale))
    end
    return b
end

build(b::LindbladianBuilder{O}) where {O} =
    Lindbladian{O}(built_hamiltonian(b.hamiltonian), copy(b.channels))

"""
    add_channel!(b::LindbladianBuilder, J, K, γ, scope, ne)

Add `γ Σ_scope 𝒟[J, K]` (under `ne`) in canonical form: jump constants become the
Hamiltonian correction `G = (iγ/2)(β̄ J - α K')` (normalization rule 2), native content of
the jumps moves into the rate (rule 3), and unmentioned binders are counted (rule I5).
"""
function add_channel!(
        b::LindbladianBuilder{O}, J::O, K::O, γ::CNum, scope::Vector{Index},
        ne::Vector{NonEqualPair},
    ) where {O <: Union{QAdd, QExpr}}
    iszero_cnum(γ) && return b
    check_jump(J)
    check_jump(K)
    J0, α = split_constant(J)
    K0, β = split_constant(K)
    if !(iszero(α) && iszero(β))
        G = (adjoint(β) * J0 - α * adjoint(K0)) * mul_cnum(CNUM_HALF, mul_cnum(CNUM_IM, γ))
        iszero(G) || add_hamiltonian!(b.hamiltonian, sum_under_scope(G, scope, ne))
    end
    (iszero(J0) || iszero(K0)) && return b
    J0, K0, rate = move_content_to_rate(J0, K0, γ)
    kept, cne, rate = count_unmentioned_binders(J0, K0, rate, scope, ne)
    iszero_cnum(rate) && return b
    add_rate!(b.channels, Channel{O}(J0, K0, kept, cne), rate)
    return b
end

function add_rate!(channels::Dict{Channel{O}, CNum}, ch::Channel{O}, rate::CNum) where {O}
    total = add_cnum(get(channels, ch, CNUM_ZERO), rate)
    iszero_cnum(total) ? delete!(channels, ch) : (channels[ch] = total)
    return channels
end

Base.zero(::Type{Lindbladian}) = zero(Lindbladian{QAdd})
Base.zero(::Type{Lindbladian{QAdd}}) = Lindbladian{QAdd}(zero_qadd(), Dict{Channel{QAdd}, CNum}())
Base.zero(::Type{Lindbladian{QExpr}}) = Lindbladian{QExpr}(qexpr_zero(), Dict{Channel{QExpr}, CNum}())
Base.zero(::L) where {L <: Lindbladian} = zero(L)
Base.iszero(L::Lindbladian) = iszero(L.hamiltonian) && isempty(L.channels)

@noinline function lindbladian_one_error()
    throw(
        ArgumentError(
            "the identity map is not a Lindbladian (every Lindbladian satisfies " *
                "`transpose(L)(1) == 0`); use `one(SuperOperator)`",
        ),
    )
end
Base.one(::Type{<:Lindbladian}) = lindbladian_one_error()
Base.one(::Lindbladian) = lindbladian_one_error()
Base.isone(::Lindbladian) = false

Base.isequal(a::Lindbladian{O}, b::Lindbladian{O}) where {O} =
    isequal(a.hamiltonian, b.hamiltonian) && isequal(a.channels, b.channels)
Base.:(==)(a::Lindbladian{O}, b::Lindbladian{O}) where {O} = isequal(a, b)
Base.hash(L::Lindbladian, h::UInt) = hash(:Lindbladian, hash(L.hamiltonian, hash(L.channels, h)))

as_generator_operator(x::QAdd) = x
as_generator_operator(x::QExpr) = x
as_generator_operator(x::QSym) = single_qadd(CNUM_ONE, Op[x])
as_generator_operator(x::Coefficient) = single_qadd(to_cnum(x), Op[])

function lindbladian_from(H::O) where {O <: Union{QAdd, QExpr}}
    return Lindbladian{O}(drop_identity(H), Dict{Channel{O}, CNum}())
end

"""
    hamiltonian_action(H)

The coherent part of a Lindblad generator with Hamiltonian `H`.

It acts as `ρ ↦ -i[H, ρ]`. The identity part of `H` is dropped, since it commutes with every state. A summed `H` keeps
its summation scope.

# Examples
```jldoctest
julia> h = FockSpace(:c); a = Destroy(h, :a);

julia> hamiltonian_action(a' * a + 1)
-i[a' * a, ·]
```

See also [`dissipator`](@ref), [`Lindbladian`](@ref).
"""
hamiltonian_action(H::Union{QField, Coefficient}) = lindbladian_from(as_generator_operator(H))

jump_type(::QAdd, ::QAdd) = QAdd
jump_type(::QField, ::QField) = QExpr

promote_jump(::Type{QAdd}, x::QAdd) = x
promote_jump(::Type{QExpr}, x::QAdd) = qexpr_from_qadd(x)
promote_jump(::Type{QExpr}, x::QExpr) = x

"""
    dissipator(J)
    dissipator(J, K)
    dissipator(J, Γ, i => j, ...)

A decay channel of a Lindblad generator with jump operator `J`.

`dissipator(J)` acts as `𝒟[J](ρ) = J ρ J' - (J'J ρ + ρ J'J)/2`, and `dissipator(J, K)` is the cross channel `𝒟[J, K](ρ) = J ρ K' - (K'J ρ + ρ K'J)/2`. Rates are
scalar prefactors, `κ * dissipator(a)`. A numeric factor of a jump moves into the rate, and
a constant in a jump moves into the Hamiltonian as the drive it generates.

`dissipator(J, Γ, i => j)` is the correlated family `Σᵢ Σⱼ Γ(i, j) 𝒟[Jᵢ, Jⱼ]`, where `i` is
a free index of `J` and `j` a second index on the same space and range. The diagonal
`i = j` is split off as `Γ(i, i) 𝒟[Jᵢ]`. Several pairs give a family over several indices.
A jump that carries its own `Σ` throws.

# Examples
```jldoctest
julia> h = FockSpace(:c); a = Destroy(h, :a);

julia> dissipator(a)(a' * a)
1 + 2 * a' * a

julia> dissipator(2a) == 4 * dissipator(a)
true
```

See also [`hamiltonian_action`](@ref), [`liouvillian`](@ref), [`Lindbladian`](@ref).
"""
dissipator(J::Union{QField, Coefficient}) = dissipator_pair(as_generator_operator(J), as_generator_operator(J))
dissipator(J::Union{QField, Coefficient}, K::Union{QField, Coefficient}) =
    dissipator_pair(as_generator_operator(J), as_generator_operator(K))

function dissipator_pair(J::QField, K::QField)
    O = jump_type(J, K)
    return channel_lindbladian(O, promote_jump(O, J), promote_jump(O, K), CNUM_ONE, Index[], EMPTY_NE)
end

function channel_lindbladian(
        ::Type{O}, J::O, K::O, γ::CNum, scope::Vector{Index}, ne::Vector{NonEqualPair},
    ) where {O <: Union{QAdd, QExpr}}
    b = LindbladianBuilder{O}()
    add_channel!(b, J, K, γ, scope, ne)
    return build(b)
end

function check_family_pair(J::QAdd, i::Index, j::Index)
    i.space_index == j.space_index || throw(
        ArgumentError("`dissipator(J, Γ, i => j)`: $(index_name(i)) and $(index_name(j)) act on different spaces"),
    )
    i.range_id == j.range_id || throw(
        ArgumentError("`dissipator(J, Γ, i => j)`: $(index_name(i)) and $(index_name(j)) have different ranges"),
    )
    mentions_index(J, i) || throw(
        ArgumentError("`dissipator(J, Γ, i => j)`: $(index_name(i)) is not a free index of the jump"),
    )
    mentions_index(J, j) && throw(
        ArgumentError("`dissipator(J, Γ, i => j)`: the partner index $(index_name(j)) already occurs in the jump"),
    )
    return nothing
end

function dissipator(J::Union{QSym, QAdd}, Γ::Coefficient, pair::Pair{Index, Index}, more::Pair{Index, Index}...)
    Jq = as_generator_operator(J)
    isempty(Jq.indices) || summed_side_error("`dissipator`")
    family = Pair{Index, Index}[pair, more...]
    for (i, j) in family
        check_family_pair(Jq, i, j)
    end
    specs = Tuple{Dict{Index, Index}, CNum, Vector{Index}, Vector{NonEqualPair}}[
        (Dict{Index, Index}(), to_cnum(Γ), Index[], EMPTY_NE),
    ]
    for (i, j) in family
        next = similar(specs, 0)
        for (partners, γ, sc, ne) in specs
            off = copy(partners)
            off[i] = j
            push!(next, (off, γ, vcat(sc, Index[i, j]), merge_ne_pair(ne, i, j)))
            push!(next, (partners, change_index(γ, j, i), vcat(sc, Index[i]), ne))
        end
        specs = next
    end
    b = LindbladianBuilder{QAdd}()
    for (partners, γ, sc, ne) in specs
        K = isempty(partners) ? Jq : change_index(Jq, partners)
        add_channel!(b, Jq, K, γ, sc, ne)
    end
    return build(b)
end

"""
    liouvillian(H, J; rates = nothing)

The Lindblad generator of the Hamiltonian `H` and the jump operators `J`, in the
QuantumOptics convention.

The generator is `-i[H, ·] + Σᵢⱼ Γ[i, j] 𝒟[J[i], J[j]]`. `rates` is `nothing` for unit
rates, a vector of rates, or a matrix `Γ` whose entry `Γ[i, j]` multiplies
`J[i] ρ J[j]'`. An empty `J` gives the coherent part only. The name is public but not
exported, since QuantumOpticsBase exports a numerical `liouvillian`.

# Examples
```jldoctest
julia> using SecondQuantizedAlgebra: liouvillian

julia> h = FockSpace(:c); a = Destroy(h, :a);

julia> liouvillian(a' * a, [a, a']; rates = [1, 2])
-i[a' * a, ·] + 𝒟[a] + 2 * 𝒟[a']
```

See also [`dissipator`](@ref), [`hamiltonian_action`](@ref).
"""
liouvillian(H::Union{QField, Coefficient}, J::AbstractVector; rates = nothing) =
    hamiltonian_action(H) + channel_family(J, rates)

channel_family(J::AbstractVector, ::Nothing) = sum(dissipator, J; init = zero(Lindbladian))

function channel_family(J::AbstractVector, rates::AbstractVector)
    length(rates) == length(J) || throw(
        ArgumentError("liouvillian: $(length(rates)) rates for $(length(J)) jump operators"),
    )
    return sum(k -> rates[k] * dissipator(J[k]), eachindex(J, rates); init = zero(Lindbladian))
end

function channel_family(J::AbstractVector, rates::AbstractMatrix)
    n = length(J)
    size(rates) == (n, n) || throw(
        ArgumentError("liouvillian: rate matrix of size $(size(rates)) for $n jump operators"),
    )
    L = zero(Lindbladian)
    for a in 1:n, b in 1:n
        iszero(rates[a, b]) && continue
        L = L + rates[a, b] * (a == b ? dissipator(J[a]) : dissipator(J[a], J[b]))
    end
    return L
end

channel_family(::AbstractVector, ::Any) =
    throw(ArgumentError("liouvillian: rates must be a vector or a matrix"))

liouvillian(H::Union{QField, Coefficient}) = hamiltonian_action(H)

function formal_lindbladian(L::Lindbladian{QAdd})
    chans = Dict{Channel{QExpr}, CNum}()
    for (ch, γ) in L.channels
        J = qexpr_from_qadd(ch.jump)
        K = is_diagonal(ch) ? J : qexpr_from_qadd(ch.partner)
        chans[Channel{QExpr}(J, K, copy(ch.scope), copy_ne(ch.ne))] = γ
    end
    return Lindbladian{QExpr}(qexpr_from_qadd(L.hamiltonian), chans)
end
formal_lindbladian(L::Lindbladian{QExpr}) = L

Base.:+(a::Lindbladian{QAdd}, b::Lindbladian{QAdd}) = add_lindbladians(a, b)
Base.:+(a::Lindbladian{QExpr}, b::Lindbladian{QExpr}) = add_lindbladians(a, b)
function add_lindbladians(a::Lindbladian{O}, b::Lindbladian{O}) where {O}
    builder = LindbladianBuilder{O}()
    accumulate!(builder, a)
    accumulate!(builder, b)
    return build(builder)
end
Base.:+(a::Lindbladian{QAdd}, b::Lindbladian{QExpr}) = formal_lindbladian(a) + b
Base.:+(a::Lindbladian{QExpr}, b::Lindbladian{QAdd}) = a + formal_lindbladian(b)
Base.:+(a::Lindbladian) = a

function scale(L::Lindbladian{O}, c::CNum) where {O}
    iszero_cnum(c) && return zero(Lindbladian{O})
    isequal(c, CNUM_ONE) && return L
    builder = LindbladianBuilder{O}()
    accumulate!(builder, L, c)
    return build(builder)
end
Base.:-(L::Lindbladian) = scale(L, CNUM_NEG1)
Base.:-(a::Lindbladian, b::Lindbladian) = a + (-b)
Base.:*(L::Lindbladian, c::Coefficient) = scale(L, to_cnum(c))
Base.:*(c::Coefficient, L::Lindbladian) = scale(L, to_cnum(c))
Base.:/(L::Lindbladian, b::Number) =
    b isa Integer && !(b isa Bool) && !iszero(b) ? L * (1 // b) : L * inv(b)
Base.:/(L::Lindbladian, b::Coefficient) = L * inv(b)
Base.://(L::Lindbladian, b::Integer) = L * (1 // b)
Base.://(L::Lindbladian, b::Coefficient) = L / b

function SuperOperator(L::Lindbladian{QAdd})
    out = SuperTermDict()
    indices = lower_hamiltonian!(out, L.hamiltonian)
    for (ch, γ) in L.channels
        lower_channel!(out, ch.jump, ch.partner, γ, ch.scope, ch.ne)
        indices = merge_unique(indices, ch.scope)
    end
    return superop_pruned(out, copy(indices))
end

function lower_hamiltonian!(out::SuperTermDict, H::QAdd)
    for (t, c) in H.arguments
        addto_key!(out, SuperTerm(copy(t.ops), Op[], copy_ne(t.ne)), mul_cnum(CNUM_NEG_IM, c))
        addto_key!(out, SuperTerm(Op[], copy(t.ops), copy_ne(t.ne)), mul_cnum(CNUM_IM, c))
    end
    return H.indices
end

function lower_channel!(
        out::SuperTermDict, J::QAdd, K::QAdd, γ::CNum, scope::Vector{Index},
        ne::Vector{NonEqualPair},
    )
    minus_half = mul_cnum(CNUM_NEG1, mul_cnum(CNUM_HALF, γ))
    for (tj, cj) in J.arguments, (tk, ck) in K.arguments
        emit_sandwich!(
            out, tj.ops, reversed_adjoint(tk.ops), mul_cnum(γ, mul_cnum(cj, conj_cnum(ck))),
            merge_ne(ne, merge_ne(tj.ne, tk.ne)), scope,
        )
    end
    for (t, c) in adjoint(K) * J
        term_ne = merge_ne(ne, t.ne)
        emit_sandwich!(out, t.ops, EMPTY_OPS, mul_cnum(minus_half, c), term_ne, scope)
        emit_sandwich!(out, EMPTY_OPS, t.ops, mul_cnum(minus_half, c), term_ne, scope)
    end
    return out
end

SuperOperator(S::SuperOperator) = S

Base.transpose(L::Lindbladian{QAdd}) = transpose(SuperOperator(L))
Base.adjoint(L::Lindbladian{QAdd}) = adjoint(SuperOperator(L))
dagger(L::Lindbladian) = adjoint(L)
qadjoint(L::Lindbladian) = adjoint(L)

(L::Lindbladian{QAdd})(ρ::Union{QSym, QAdd, Coefficient}) = SuperOperator(L)(ρ)

Base.:+(L::Lindbladian{QAdd}, c::Coefficient) = SuperOperator(L) + c
Base.:+(c::Coefficient, L::Lindbladian{QAdd}) = SuperOperator(L) + c
Base.:-(L::Lindbladian{QAdd}, c::Coefficient) = SuperOperator(L) - c
Base.:-(c::Coefficient, L::Lindbladian{QAdd}) = c - SuperOperator(L)

Base.inv(::Lindbladian{QAdd}) = negative_power_error()

function Base.:^(L::Lindbladian{QAdd}, n::Integer)
    n >= 0 || negative_power_error()
    return SuperOperator(L)^n
end

function Σ(L::Lindbladian{QAdd}, i::Index, non_equal::Vector{Index} = Index[])
    i.range_id == 0 && throw(ArgumentError("Σ: index $(index_name(i)) has no summation range"))
    i in L.hamiltonian.indices && shared_binder_error(i)
    b = LindbladianBuilder{QAdd}()
    iszero(L.hamiltonian) || add_hamiltonian!(b.hamiltonian, Σ(L.hamiltonian, i, non_equal))
    for (ch, γ) in L.channels
        i in ch.scope && shared_binder_error(i)
        ne_aug = merge_ne_pairs(ch.ne, i, non_equal)
        scope_i = vcat(ch.scope, Index[i])
        if !channel_mentions(ch.jump, ch.partner, γ, i)
            add_channel!(b, ch.jump, ch.partner, γ, scope_i, ne_aug)
            continue
        end
        ops = Op[o for q in (ch.jump, ch.partner) for t in keys(q.arguments) for o in t.ops]
        diag = coinciding_pairs(i, ops, ne_aug)
        add_channel!(b, ch.jump, ch.partner, γ, scope_i, merge_ne_pairs(ne_aug, diag))
        for (_, ext) in diag
            add_channel!(
                b, change_index(ch.jump, i, ext),
                change_index(ch.partner, i, ext), change_index(γ, i, ext), copy(ch.scope),
                diagonal_ne(ne_aug, diag, i, ext),
            )
        end
    end
    return build(b)
end

Σ(L::Lindbladian, i::Index, j::Index, rest::Index...) = Σ(Σ(L, i), j, rest...)
