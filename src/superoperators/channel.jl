"""
    Channel{O}

One decay channel of a [`Lindbladian`](@ref), without its rate.

A channel `Σ 𝒟[J, K]` has a jump `J`, a partner `K` and the summation indices of its sum.
An independent channel `𝒟[J]` has `K == J`. `channels(L)` pairs each channel with its rate,
and constants and numeric factors of a jump are moved out of the channel when it is built.

# Examples
```jldoctest
julia> using SecondQuantizedAlgebra: channels

julia> h = FockSpace(:c); a = Destroy(h, :a);

julia> first(channels(dissipator(2a)))
𝒟[a] => 4
```

See also [`jump`](@ref), [`hamiltonian`](@ref).
"""
struct Channel{O <: Union{QAdd, QExpr}}
    jump::O
    partner::O
    scope::Vector{Index}
    ne::Vector{NonEqualPair}
    hash::UInt
    function Channel{O}(
            jump::O, partner::O, scope::Vector{Index}, ne::Vector{NonEqualPair},
        ) where {O <: Union{QAdd, QExpr}}
        h = hash(:Channel, hash(jump, hash(partner, hash(scope, hash(ne, zero(UInt))))))
        return new{O}(jump, partner, scope, ne, h)
    end
end

Base.isequal(a::Channel, b::Channel) = a.hash == b.hash &&
    isequal(a.jump, b.jump) && isequal(a.partner, b.partner) &&
    isequal(a.scope, b.scope) && isequal(a.ne, b.ne)
Base.:(==)(a::Channel, b::Channel) = isequal(a, b)
Base.hash(c::Channel, h::UInt) = hash(c.hash, h)

"""
    jump(c::Channel)
    partner(c::Channel)
    scope(c::Channel)

The jump operator, the partner and the summation indices of a decay channel.

For a channel `Σ 𝒟[J, K]` these are `J`, `K` and the indices of the sum. The partner of an
independent channel `𝒟[J]` is `J` itself.

# Examples
```jldoctest
julia> using SecondQuantizedAlgebra: channels, jump, partner, scope

julia> h = FockSpace(:c); a = Destroy(h, :a);

julia> ch, rate = first(channels(dissipator(a, a')));

julia> jump(ch), partner(ch), isempty(scope(ch))
(a, a', true)
```

See also [`Channel`](@ref), [`hamiltonian`](@ref).
"""
jump(c::Channel) = c.jump
partner(c::Channel) = c.partner
scope(c::Channel) = c.scope
@doc (@doc jump) partner
@doc (@doc jump) scope

is_diagonal(c::Channel) = c.partner === c.jump || isequal(c.partner, c.jump)

mentions_index(q::QExpr, idx::Index) = idx in get_indices(q)

function split_constant(J::QAdd)
    constant = QTermDict()
    rest = QTermDict()
    for (t, c) in J.arguments
        isempty(t.ops) ? (constant[copy_key(t)] = c) : (rest[copy_key(t)] = c)
    end
    isempty(constant) && return (J, zero_qadd())
    return (QAdd(rest, copy(J.indices)), QAdd(constant, Index[]))
end

function split_constant(J::QExpr)
    qexpr_is_scalar(J) && return (qexpr_zero(), single_qadd(J.coeff, Op[]))
    if has_polynomial_body(J)
        body, constant = split_constant(only(J.args)::QAdd)
        return (qexpr_from_qadd(body), constant)
    end
    J.kind == QEXPR_ADD || return (J, zero_qadd())
    args = QExprArg[]
    constant = zero_qadd()
    for arg in J.args
        if arg isa QAdd
            body, c = split_constant(arg)
            constant = constant + c
            iszero(body) || push!(args, body)
        else
            push!(args, arg)
        end
    end
    iszero(constant) && return (J, constant)
    return (qexpr_sum(args), constant)
end

drop_identity(H::QAdd) = first(split_constant(H))
drop_identity(H::QExpr) = first(split_constant(H))

function native_content(c::CNum)
    is_native(c) && return c
    tail = c.tail
    if tail isa Poly && length(tail.terms) == 1
        s = scalar_coeff(monomial_scalar(only(tail.terms)))
        u = c * inv(s)
        return isequal(u * conj(u), CNUM_ONE) ? c : s
    end
    return CNUM_ONE
end

function lead_coefficient(J::QAdd)
    lead = argmin(term_order_key, keys(J.arguments))
    return J.arguments[lead]
end

function jump_content(J::QAdd)
    iszero(J) && return CNUM_ONE
    return native_content(lead_coefficient(J))
end

function jump_content(J::QExpr)
    iszero(J) && return CNUM_ONE
    has_polynomial_body(J) && return jump_content(only(J.args)::QAdd)
    J.kind == QEXPR_ADD && return CNUM_ONE
    return native_content(J.coeff)
end

scale_jump(J::QAdd, c::CNum) = J * c
scale_jump(J::QExpr, c::CNum) = qexpr_scale(J, c)

hamiltonian_accumulator(::Type{QAdd}) = QAddBuilder()
hamiltonian_accumulator(::Type{QExpr}) = QExprArg[]

add_hamiltonian!(b::QAddBuilder, H::QAdd) = accumulate!(b, H)
add_hamiltonian!(b::Vector{QExprArg}, H::QExprArg) = (iszero(H) || push!(b, H); b)

built_hamiltonian(b::QAddBuilder) = drop_identity(build(b))
built_hamiltonian(b::Vector{QExprArg}) = drop_identity(qexpr_sum(b))

"""
    sum_under_scope(G::QAdd, scope, ne) -> QAdd

`Σ_scope G` under the constraints `ne`, with every binder that `G` does not mention counted
into the coefficient.
"""
function sum_under_scope(G::QAdd, scope::Vector{Index}, ne::Vector{NonEqualPair})
    isempty(scope) && isempty(ne) && return G
    isempty(G.indices) || summed_side_error("a dissipator")
    out = QTermDict()
    for (t, c) in G.arguments
        temp = QTermDict()
        canonicalize!(temp, copy(t.ops), c, merge_ne(t.ne, ne))
        settle_scope!(out, temp, scope)
    end
    return QAdd(out, drop_unused_indices(out, copy(scope)))
end

@noinline function formal_scope_error()
    throw(
        ArgumentError(
            "a formal (QExpr) jump or Hamiltonian cannot carry a bound summation scope; " *
                "lower it with `taylor(expr, 0:n)` first",
        ),
    )
end

function sum_under_scope(G::QExpr, scope::Vector{Index}, ::Vector{NonEqualPair})
    isempty(scope) || formal_scope_error()
    return G
end

function channel_mentions(J, K, γ::CNum, idx::Index)
    return mentions_index(J, idx) || mentions_index(K, idx) ||
        depends_on_index_term(γ, EMPTY_OPS, idx)
end

function check_jump(J::QAdd)
    isempty(J.indices) || summed_side_error("a dissipator")
    return J
end
function check_jump(J::QExpr)
    for arg in J.args
        arg isa QAdd && (isempty(arg.indices) || summed_side_error("a dissipator"))
    end
    return J
end

"""
    move_content_to_rate(J, K, γ) -> (J, K, rate)

Normalization rule 3: divide the native content out of each jump and multiply it into the
rate, `γ 𝒟[cJ, dK] = γ c d̄ 𝒟[J, K]`. A float rescaling that would overflow or underflow a
coefficient of the lowered channel is skipped, so the channel keeps its jumps and rate.
"""
function move_content_to_rate(J::Union{QAdd, QExpr}, K::Union{QAdd, QExpr}, γ::CNum)
    cj = jump_content(J)
    ck = jump_content(K)
    J0, K0, rate = J, K, γ
    if !isequal(cj, CNUM_ONE)
        J0 = scale_jump(J, inv(cj))
        rate = rate * cj
    end
    if !isequal(ck, CNUM_ONE)
        K0 = scale_jump(K, inv(ck))
        rate = rate * conj(ck)
    end
    lowered_well_scaled(J0, K0, rate) || return (J, isequal(J, K) ? J : K, γ)
    return (J0, isequal(J0, K0) ? J0 : K0, rate)
end

well_scaled(z::Union{AbstractFloat, Complex{<:AbstractFloat}}) =
    isfinite(z) && (iszero(z) || abs(z) >= floatmin(Float64))
well_scaled(c::CNum) = !is_native(c) || well_scaled(native_scalar(c))
well_scaled(_) = true

jump_coefficients(J::QAdd) = collect(values(J.arguments))
jump_coefficients(J::QExpr) =
    has_polynomial_body(J) ? jump_coefficients(only(J.args)::QAdd) : CNum[J.coeff]

function lowered_well_scaled(J::Union{QAdd, QExpr}, K::Union{QAdd, QExpr}, rate::CNum)
    (well_scaled(rate) && !iszero(rate)) || return false
    for a in jump_coefficients(J), b in jump_coefficients(K)
        well_scaled(rate * (a * conj(b))) || return false
    end
    return true
end

"""
    count_unmentioned_binders(J, K, rate, scope, ne) -> (kept, ne, rate)

Normalization rule I5: keep the constraints that touch the channel, and count every binder
the channel does not mention into the rate when its constrained count is known.
"""
function count_unmentioned_binders(
        J::Union{QAdd, QExpr}, K::Union{QAdd, QExpr}, rate::CNum, scope::Vector{Index},
        ne::Vector{NonEqualPair},
    )
    channel_ne = NonEqualPair[]
    for p in canonical_ne(ne)
        x, y = p
        referenced = x in scope || y in scope ||
            channel_mentions(J, K, rate, x) || channel_mentions(J, K, rate, y)
        referenced && push!(channel_ne, p)
    end
    cne = isempty(channel_ne) ? EMPTY_NE : channel_ne
    mentions(γ, idx) = channel_mentions(J, K, γ, idx)
    rate, cne = eliminate_binders(mentions, rate, cne, scope)
    kept = Index[]
    for idx in scope
        idx in kept && continue
        (mentions(rate, idx) || ne_mentions(cne, idx)) && push!(kept, idx)
    end
    return (sort!(kept), cne, rate)
end

ne_mentions(ne::Vector{NonEqualPair}, idx::Index) = any(p -> p[1] == idx || p[2] == idx, ne)
