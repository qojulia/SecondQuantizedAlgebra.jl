linear_transform(basis::Vector{Op}, linear::Matrix{CNum}) =
    canonical_transform(AffineAction(basis, linear, fill(CNUM_ZERO, length(basis))))

generated_frame(G::QAdd, θ::Real, t::Num) = timed_transform(exact_generator_transform(G, θ), gauge(G, θ, t), t)

function shift_transform(basis::Vector{Op}, shift::Vector{CNum})
    slots = Dict{SiteKey, Int}()
    bases = Vector{Op}[]
    shifts = Vector{CNum}[]
    for (operator, c) in zip(basis, shift)
        k = get!(slots, site_key(operator), length(bases) + 1)
        k > length(bases) && (push!(bases, Op[]); push!(shifts, CNum[]))
        push!(bases[k], operator)
        push!(shifts[k], c)
    end
    blocks = AffineBlock[
        AffineBlock(infer_affine_structure(b), b, identity_linear(length(b)), c)
            for (b, c) in zip(bases, shifts)
    ]
    return canonical_transform(validated_action(blocks, ParamRelation[]))
end

function canonical_partners(basis::Vector{Op})
    partners = zeros(Int, length(basis))
    inverse_weights = fill(CNUM_ZERO, length(basis))
    for (i, operator) in enumerate(basis), (j, other) in enumerate(basis)
        for (_, weight) in commutator(other, operator)
            partners[i] = j
            inverse_weights[i] = inv(weight)
        end
    end
    return partners, inverse_weights
end

function displacement_gauge(
        basis::Vector{Op}, shift::Vector{CNum}, velocity::Vector{CNum}, scalar_velocity::Vector{CNum},
    )
    partners, inverse_weights = canonical_partners(basis)
    terms = Tuple{CNum, Vector{Op}}[]
    scalar = CNUM_ZERO
    for (i, partner) in enumerate(partners)
        push!(terms, ((CNUM_IM * inverse_weights[i]) * velocity[i], Op[basis[partner]]))
        scalar = scalar + (CNUM_NEG_IM * CNUM_HALF * inverse_weights[i]) * (shift[i] * scalar_velocity[partner])
    end
    push!(terms, (scalar, Op[]))
    return rule_qadd(terms)
end

displacement_gauge(basis::Vector{Op}, shift::Vector{CNum}, velocity::Vector{CNum}) =
    displacement_gauge(basis, shift, velocity, velocity)

"""
    Displace(a, α)
    Displace(x, p, dx, dp)
    Displace(a, α, t)
    Displace(x, p, dx, dp, t)

Displacement of a bosonic mode or of a pair of quadratures.

The mode `a` is shifted by the complex amplitude `α`, ``a \\mapsto a + \\alpha``. A pair of
quadratures is shifted by real amounts, ``x \\mapsto x + dx`` and ``p \\mapsto p + dp``.
With a final time variable `t`, the amplitudes may depend on time and the transformation
carries the [`gauge_term`](@ref) of the moving frame.

# Examples
```jldoctest
julia> a = Destroy(FockSpace(:cavity), :a);

julia> @variables α::Number;

julia> conjugate(a' * a, Displace(a, α))
conj(α)*α + conj(α) * a + α * a' + a' * a

julia> h = PhaseSpace(:mode); x = Position(h, :x); p = Momentum(h, :p);

julia> @variables dx dp;

julia> conjugate(x * p, Displace(x, p, dx, dp))
dp*dx + dp * x + dx * p + x * p
```

See also [`DisplacementFrame`](@ref), [`Rotation`](@ref), [`Squeeze`](@ref).
"""
mode_basis(a::Op, what::AbstractString) = (d = fock_or_throw(a, what); Op[d, adjoint(d)])
mode_shift(c::CNum) = CNum[c, conj(c)]

Displace(a::Op, α::Coefficient) = shift_transform(mode_basis(a, "`Displace`"), mode_shift(to_cnum(α)))

function Displace(a::Op, α::Coefficient, t::Num)
    basis = mode_basis(a, "`Displace`")
    tt = time_or_throw(t)
    c = to_cnum(α)
    gauge = displacement_gauge(basis, mode_shift(c), mode_shift(dt(c, tt)))
    return timed_transform(shift_transform(basis, mode_shift(c)), gauge, tt)
end

function number_generator(a::Op)
    d = fock_or_throw(a, "`Rotation`")
    return adjoint(d) * d
end

"""
    Rotation(a, θ)
    Rotation(a, b, θ)
    Rotation(x, p, θ)
    Rotation(S, axis::Int, θ)
    Rotation(σ, W::AbstractMatrix)
    Rotation(a, θ, t)
    Rotation(a, b, θ, t)
    Rotation(x, p, θ, t)
    Rotation(S, axis::Int, θ, t)
    Rotation(σ, W::AbstractMatrix, t)

Rotation of a mode, of a pair of modes or quadratures, of a spin, or of the levels of a finite-level system.

- `Rotation(a, θ)` turns the phase of a mode, ``a \\mapsto e^{-i\\theta} a``.
- `Rotation(a, b, θ)` mixes two modes as a beam splitter,
  ``a \\mapsto \\cos\\theta\\, a + \\sin\\theta\\, b`` and
  ``b \\mapsto -\\sin\\theta\\, a + \\cos\\theta\\, b``.
- `Rotation(x, p, θ)` rotates a pair of quadratures in phase space in the same way.
- `Rotation(S, axis, θ)` rotates a spin or Pauli triple by `θ` about `axis` (1, 2 or 3);
  `S` is any component of the triple.
- `Rotation(σ, W)` changes the basis of a finite-level system,
  ``\\sigma^{ij} \\mapsto W^\\dagger \\sigma^{ij} W``, acting on all transitions of
  the system; `σ` is any one of them. `W` is a square unitary matrix with one row per
  level; its size is checked, its unitarity is not.

With a final time variable `t`, the angle or matrix may depend on time and the
transformation carries the [`gauge_term`](@ref) of the moving frame.

# Examples
```jldoctest
julia> a = Destroy(FockSpace(:cavity), :a);

julia> @variables θ;

julia> conjugate(a, Rotation(a, θ))
exp(-im*θ) * a

julia> Sx = Spin(SpinSpace(:S), :S, 1);

julia> conjugate(Sx, Rotation(Sx, 3, θ))
cos(θ) * Sx - sin(θ) * Sy
```

See also [`RotatingFrame`](@ref), [`Squeeze`](@ref), [`Displace`](@ref).
"""
Rotation(a::Op, θ::Real) = exact_generator_transform(number_generator(a), θ)
Rotation(a::Op, θ::Num, t::Num) = generated_frame(number_generator(a), θ, time_or_throw(t))

"""
    Squeeze(a, r, ϕ = 0)
    Squeeze(a, b, r)
    Squeeze(x, p, r)
    Squeeze(a, r, ϕ, t)
    Squeeze(a, b, r, t)
    Squeeze(x, p, r, t)

Squeezing of a mode, of a pair of modes, or of a pair of quadratures.

- `Squeeze(a, r, ϕ)` squeezes one mode,
  ``a \\mapsto \\cosh r\\, a + e^{i\\phi} \\sinh r\\, a^\\dagger``.
- `Squeeze(a, b, r)` squeezes two modes jointly,
  ``a \\mapsto \\cosh r\\, a + \\sinh r\\, b^\\dagger`` and the same with `a` and `b`
  exchanged.
- `Squeeze(x, p, r)` rescales a pair of quadratures,
  ``x \\mapsto e^{r} x`` and ``p \\mapsto e^{-r} p``.

With a final time variable `t`, the parameters may depend on time and the transformation
carries the [`gauge_term`](@ref) of the moving frame. A time-dependent single-mode squeeze
always gives the phase, `Squeeze(a, r, 0, t)`, since the third argument is `ϕ`.

# Examples
```jldoctest
julia> a = Destroy(FockSpace(:cavity), :a);

julia> @variables r;

julia> conjugate(a, Squeeze(a, r))
cosh(r) * a + sinh(r) * a'

julia> conjugate(a, inv(Squeeze(a, r)))
cosh(r) * a - sinh(r) * a'
```

See also [`Bogoliubov`](@ref), [`Rotation`](@ref).
"""
function Squeeze(a::Op, r::Real, ϕ::Real = 0)
    phase_is_time(r, ϕ) && unitary_error(
        "the third argument of `Squeeze(a, r, ϕ)` is the squeezing phase, but `$ϕ` also " *
            "appears in `r`; a time-dependent squeeze is written `Squeeze(a, r, ϕ, t)`",
    )
    return single_mode_squeeze(a, r, ϕ)
end

phase_is_time(::Real, ::Real) = false
function phase_is_time(r::Num, ϕ::Num)
    raw = SymbolicUtils.unwrap(ϕ)
    SymbolicUtils.issym(raw) || return false
    return any(v -> isequal(v, raw), get_variables(r))
end

function single_mode_squeeze(a::Op, r::Real, ϕ::Real)
    d = fock_or_throw(a, "`Squeeze`")
    rate = CNUM_HALF * CNUM_IM * phase(ϕ)
    G = rule_qadd((rate, Op[adjoint(d), adjoint(d)]), (conj(rate), Op[d, d]))
    return exact_generator_transform(G, r)
end

function squeeze_gauge(d::Op, r::Real, ϕ::Real, t::Num)
    sh = to_cnum(sinh(r))
    ch = to_cnum(cosh(r))
    imaginary_dr = CNUM_IM * dt(r, t)
    dϕ = dt(ϕ, t)
    diagonal = dϕ * (sh * sh)
    mixing = dϕ * (sh * ch)
    rotation = phase(ϕ)
    return rule_qadd(
        (diagonal, Op[adjoint(d), d]),
        (CNUM_HALF * diagonal, Op[]),
        (
            CNUM_HALF * (rotation * (mixing - imaginary_dr)),
            Op[adjoint(d), adjoint(d)],
        ),
        (
            CNUM_HALF * (conj(rotation) * (mixing + imaginary_dr)),
            Op[d, d],
        ),
    )
end

function Squeeze(a::Op, r::Real, ϕ::Real, t::Num)
    d = fock_or_throw(a, "`Squeeze`")
    tt = time_or_throw(t)
    return timed_transform(single_mode_squeeze(a, r, ϕ), squeeze_gauge(d, r, ϕ, tt), tt)
end

function two_modes(a::Op, b::Op, what::AbstractString)
    x = fock_or_throw(a, what)
    y = fock_or_throw(b, what)
    site_key(x) == site_key(y) && unitary_error("$what needs two distinct modes")
    return x, y
end

function phase_pair(x::Op, p::Op, what::AbstractString)
    (is_position(x) && is_momentum(p)) || unitary_error(
        "$what expects a `(Position, Momentum)` pair in that order",
    )
    site_key(x) == site_key(p) || unitary_error(
        "$what expects both quadratures on the same site",
    )
    return (x, p)
end

function rotation_generator(a::Op, b::Op)
    if is_phase_space(a)
        phase_pair(a, b, "`Rotation`")
        return (a * a + b * b) * (1 // 2)
    end
    x, y = two_modes(a, b, "`Rotation`")
    return im * (adjoint(x) * y - adjoint(y) * x)
end

function squeeze_generator(a::Op, b::Op)
    if is_phase_space(a)
        phase_pair(a, b, "`Squeeze`")
        return (a * b + b * a) * (1 // 2)
    end
    x, y = two_modes(a, b, "`Squeeze`")
    return im * (adjoint(x) * adjoint(y) - y * x)
end

Rotation(a::Op, b::Op, θ::Real) = exact_generator_transform(rotation_generator(a, b), θ)
Rotation(a::Op, b::Op, θ::Real, t::Num) = generated_frame(rotation_generator(a, b), θ, time_or_throw(t))
Squeeze(a::Op, b::Op, r::Real) = exact_generator_transform(squeeze_generator(a, b), r)
Squeeze(a::Op, b::Op, r::Real, t::Num) = generated_frame(squeeze_generator(a, b), r, time_or_throw(t))

function Displace(x::Op, p::Op, dx::Real, dp::Real)
    phase_pair(x, p, "`Displace`")
    return shift_transform(Op[x, p], CNum[to_cnum(dx), to_cnum(dp)])
end

function Displace(x::Op, p::Op, dx::Real, dp::Real, t::Num)
    phase_pair(x, p, "`Displace`")
    tt = time_or_throw(t)
    shift = CNum[to_cnum(dx), to_cnum(dp)]
    gauge = displacement_gauge(Op[x, p], shift, CNum[dt(c, tt) for c in shift])
    return timed_transform(shift_transform(Op[x, p], shift), gauge, tt)
end

function spin_generator(S::Op, axis::Integer)
    (is_pauli(S) || is_spin(S)) || unitary_error(
        "`Rotation` expects a Pauli or Spin operator, got $(S.kind)",
    )
    1 <= axis <= 3 || unitary_error("axis must be 1, 2, or 3, got $axis")
    return single_qadd(is_pauli(S) ? CNUM_HALF : CNUM_ONE, Op[with_levels(S, Int(axis), 0)])
end

Rotation(S::Op, axis::Integer, θ::Real) = exact_generator_transform(spin_generator(S, axis), θ)
Rotation(S::Op, axis::Integer, θ::Real, t::Num) = generated_frame(spin_generator(S, axis), θ, time_or_throw(t))

function nlevel_or_throw(σ::Op, W::AbstractMatrix)
    is_transition(σ) || unitary_error(
        "`Rotation` expects an ordinary `Transition` operator, got $(σ.kind)",
    )
    n = Int(σ.nlev)
    size(W) == (n, n) || unitary_error(
        "`W` must be $n×$n for a $n-level space, got $(size(W))",
    )
    return n
end

function nlevel_rotation(σ::Op, W::AbstractMatrix)
    n = nlevel_or_throw(σ, W)
    coefficients = CNum[to_cnum(W[i, j]) for i in 1:n, j in 1:n]
    linear = fill(CNUM_ZERO, n * n, n * n)
    for i in 1:n, j in 1:n, k in 1:n, l in 1:n
        coefficient = conj(coefficients[i, k]) * coefficients[j, l]
        iszero(coefficient) || (linear[(i - 1) * n + j, (k - 1) * n + l] = coefficient)
    end
    return linear_transform(site_generators(σ), linear), coefficients
end

Rotation(σ::Op, W::AbstractMatrix) = first(nlevel_rotation(σ, W))

function nlevel_gauge(σ::Op, W::Matrix{Coeff}, t::Num)
    n = size(W, 1)
    conjugated_derivatives = Matrix{Coeff}(undef, n, n)
    for j in 1:n, k in 1:n
        conjugated_derivatives[k, j] = conj(dt(W[k, j], t))
    end
    gauge = QTermDict()
    for j in 1:n, l in 1:n
        coefficient = CNUM_ZERO
        for k in 1:n
            coefficient = coefficient + conjugated_derivatives[k, j] * W[k, l]
        end
        coefficient = CNUM_IM * coefficient
        iszero(coefficient) ||
            addto!(gauge, Op[with_levels(σ, j, l)], coefficient)
    end
    return QAdd(gauge, EMPTY_INDICES)
end

function Rotation(σ::Op, W::AbstractMatrix, t::Num)
    tt = time_or_throw(t)
    U, coefficients = nlevel_rotation(σ, W)
    return timed_transform(U, nlevel_gauge(σ, coefficients, tt), tt)
end
