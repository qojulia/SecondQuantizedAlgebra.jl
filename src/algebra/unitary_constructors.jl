linear_transform(basis::Vector{Op}, linear::Matrix{CNum}, relations::Vector{ParamRelation} = ParamRelation[]) =
    canonical_transform(AffineAction(basis, linear, fill(CNUM_ZERO, length(basis)), relations))

fock_displacement(d::Op, c::CNum) =
    canonical_transform(AffineAction(Op[d, adjoint(d)], identity_linear(2), CNum[c, conj(c)]))

function fock_displacement_gauge(
        d::Op, c::CNum, derivative::CNum, derivative_adjoint::CNum,
        raising_coefficient::CNum, lowering_coefficient::CNum,
    )
    scalar = CNUM_NEG_IM * (CNUM_HALF * (conj(c) * derivative - c * derivative_adjoint))
    return rule_qadd(
        (raising_coefficient, Op[adjoint(d)]),
        (lowering_coefficient, Op[d]),
        (scalar, Op[]),
    )
end

fock_displacement_gauge(d::Op, c::CNum, derivative::CNum) = fock_displacement_gauge(
    d, c, derivative, conj(derivative), CNUM_NEG_IM * derivative, CNUM_IM * conj(derivative),
)

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
function Displace(a::Op, α::Coefficient)
    d = fock_or_throw(a, "`Displace`")
    return fock_displacement(d, to_cnum(α))
end

function Displace(a::Op, α::Coefficient, t::Num)
    d = fock_or_throw(a, "`Displace`")
    tt = time_or_throw(t)
    c = to_cnum(α)
    gauge = fock_displacement_gauge(d, c, dt(c, tt))
    return timed_transform(fock_displacement(d, c), gauge, tt)
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
function Rotation(a::Op, θ::Real)
    d = fock_or_throw(a, "`Rotation`")
    return linear_transform(Op[d, adjoint(d)], CNum[conj_phase(θ) CNUM_ZERO; CNUM_ZERO phase(θ)])
end

function Rotation(a::Op, θ::Num, t::Num)
    d = fock_or_throw(a, "`Rotation`")
    tt = time_or_throw(t)
    return timed_transform(Rotation(a, θ), gauge(adjoint(d) * d, θ, tt), tt)
end

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
    ch = to_cnum(cosh(r))
    sh = phase(ϕ) * to_cnum(sinh(r))
    return linear_transform(Op[d, adjoint(d)], CNum[ch sh; conj(sh) conj(ch)], ParamRelation[hyp_rel(r)])
end

function squeeze_gauge(d::Op, r::Real, ϕ::Real, t::Num)
    sh = to_cnum(sinh(r))
    ch = to_cnum(cosh(r))
    imaginary_dr = CNUM_IM * dt(r, t)
    dϕ = dt(ϕ, t)
    diagonal = dϕ * (sh * sh)
    mixing = dϕ * (sh * ch)
    return rule_qadd(
        (diagonal, Op[adjoint(d), d]),
        (CNUM_HALF * diagonal, Op[]),
        (
            CNUM_HALF * (phase(ϕ) * (mixing - imaginary_dr)),
            Op[adjoint(d), adjoint(d)],
        ),
        (
            CNUM_HALF * (conj_phase(ϕ) * (mixing + imaginary_dr)),
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

function fock_mode_rotation(a::Op, b::Op, θ::Real)
    x, y = two_modes(a, b, "`Rotation`")
    c = to_cnum(cos(θ))
    s = to_cnum(sin(θ))
    return linear_transform(
        Op[x, y, adjoint(x), adjoint(y)],
        CNum[
            c s CNUM_ZERO CNUM_ZERO
            -s c CNUM_ZERO CNUM_ZERO
            CNUM_ZERO CNUM_ZERO conj(c) conj(s)
            CNUM_ZERO CNUM_ZERO -conj(s) conj(c)
        ],
        ParamRelation[trig_rel(θ)],
    )
end

function quadrature_rotation(x::Op, p::Op, θ::Real)
    phase_pair(x, p, "`Rotation`")
    c = to_cnum(cos(θ))
    s = to_cnum(sin(θ))
    return linear_transform(Op[x, p], CNum[c s; -s c], ParamRelation[trig_rel(θ)])
end

Rotation(a::Op, b::Op, θ::Real) =
    is_phase_space(a) ? quadrature_rotation(a, b, θ) : fock_mode_rotation(a, b, θ)

function Rotation(a::Op, b::Op, θ::Real, t::Num)
    tt = time_or_throw(t)
    if is_phase_space(a)
        U = quadrature_rotation(a, b, θ)
        generator = (a * a + b * b) * (1 // 2)
        return timed_transform(U, gauge(generator, θ, tt), tt)
    end
    x, y = two_modes(a, b, "`Rotation`")
    generator = im * (adjoint(x) * y - adjoint(y) * x)
    return timed_transform(fock_mode_rotation(x, y, θ), gauge(generator, θ, tt), tt)
end

function fock_two_mode_squeeze(a::Op, b::Op, r::Real)
    x, y = two_modes(a, b, "`Squeeze`")
    u = to_cnum(cosh(r))
    v = to_cnum(sinh(r))
    return linear_transform(
        Op[x, y, adjoint(x), adjoint(y)],
        CNum[
            u CNUM_ZERO CNUM_ZERO v
            CNUM_ZERO u v CNUM_ZERO
            CNUM_ZERO conj(v) conj(u) CNUM_ZERO
            conj(v) CNUM_ZERO CNUM_ZERO conj(u)
        ],
        ParamRelation[hyp_rel(r)],
    )
end

function quadrature_squeeze(x::Op, p::Op, r::Real)
    phase_pair(x, p, "`Squeeze`")
    up = to_cnum(exp(r))
    down = to_cnum(inv(exp(r)))
    return linear_transform(Op[x, p], CNum[up CNUM_ZERO; CNUM_ZERO down])
end

Squeeze(a::Op, b::Op, r::Real) =
    is_phase_space(a) ? quadrature_squeeze(a, b, r) : fock_two_mode_squeeze(a, b, r)

function Squeeze(a::Op, b::Op, r::Real, t::Num)
    tt = time_or_throw(t)
    if is_phase_space(a)
        U = quadrature_squeeze(a, b, r)
        generator = (a * b + b * a) * (1 // 2)
        return timed_transform(U, gauge(generator, r, tt), tt)
    end
    x, y = two_modes(a, b, "`Squeeze`")
    generator = im * (adjoint(x) * adjoint(y) - y * x)
    return timed_transform(fock_two_mode_squeeze(x, y, r), gauge(generator, r, tt), tt)
end

quadrature_displacement(x::Op, p::Op, cx::CNum, cp::CNum) =
    canonical_transform(AffineAction(Op[x, p], identity_linear(2), CNum[cx, cp]))

function quadrature_displacement_gauge(
        x::Op, p::Op, cx::CNum, cp::CNum,
        derivative_x::CNum, derivative_p::CNum,
        x_coefficient::CNum, p_coefficient::CNum,
    )
    scalar = -CNUM_HALF * (cp * derivative_x - cx * derivative_p)
    return rule_qadd(
        (x_coefficient, Op[x]), (p_coefficient, Op[p]), (scalar, Op[]),
    )
end

function Displace(x::Op, p::Op, dx::Real, dp::Real)
    phase_pair(x, p, "`Displace`")
    return quadrature_displacement(x, p, to_cnum(dx), to_cnum(dp))
end

function Displace(x::Op, p::Op, dx::Real, dp::Real, t::Num)
    phase_pair(x, p, "`Displace`")
    tt = time_or_throw(t)
    cx = to_cnum(dx)
    cp = to_cnum(dp)
    derivative_x = dt(cx, tt)
    derivative_p = dt(cp, tt)
    gauge = quadrature_displacement_gauge(
        x, p, cx, cp, derivative_x, derivative_p,
        derivative_p, -derivative_x,
    )
    return timed_transform(quadrature_displacement(x, p, cx, cp), gauge, tt)
end

axis_or_throw(axis::Integer) = 1 <= axis <= 3 ? Int32(axis) :
    unitary_error("axis must be 1, 2, or 3, got $axis")

function triple_or_throw(S::Op)
    (is_pauli(S) || is_spin(S)) || unitary_error(
        "`Rotation` expects a Pauli or Spin operator, got $(S.kind)",
    )
    return S
end

function Rotation(S::Op, axis::Integer, θ::Real)
    triple_or_throw(S)
    main_axis = Int(axis_or_throw(axis))
    u = with_levels(S, mod1(main_axis + 1, 3), 0)
    v = with_levels(S, mod1(main_axis + 2, 3), 0)
    fixed = with_levels(S, main_axis, 0)
    c = to_cnum(cos(θ))
    s = to_cnum(sin(θ))
    return linear_transform(
        Op[u, v, fixed],
        CNum[c -s CNUM_ZERO; s c CNUM_ZERO; CNUM_ZERO CNUM_ZERO CNUM_ONE],
        ParamRelation[trig_rel(θ)],
    )
end

function Rotation(S::Op, axis::Integer, θ::Real, t::Num)
    triple_or_throw(S)
    tt = time_or_throw(t)
    fixed = with_levels(S, Int(axis_or_throw(axis)), 0)
    generator = scaled(S.kind === OP_PAULI ? CNUM_HALF : CNUM_ONE, fixed)
    return timed_transform(Rotation(S, axis, θ), gauge(generator, θ, tt), tt)
end

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
    nlevel_or_throw(σ, W)
    n = size(W, 1)
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
