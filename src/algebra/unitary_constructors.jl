include("unitary_affine.jl")
include("unitary_affine_provenance.jl")

# === Fock transformations ===

function fock_displacement(d::Op, c::CNum)
    action = AffineAction(
        Op[d, adjoint(d)],
        CNum[
            CNUM_ONE CNUM_ZERO
            CNUM_ZERO CNUM_ONE
        ],
        CNum[c, conj_cnum(c)],
    )
    return canonical_transform(action)
end

function fock_displacement_gauge(
        d::Op, c::CNum, derivative::CNum, derivative_adjoint::CNum,
        raising_coefficient::CNum, lowering_coefficient::CNum,
    )
    scalar = mul_cnum(
        CNUM_NEG_IM,
        mul_cnum(
            CNUM_HALF,
            add_cnum(
                mul_cnum(conj_cnum(c), derivative),
                neg_cnum(mul_cnum(c, derivative_adjoint)),
            ),
        ),
    )
    return rule_qadd(
        (raising_coefficient, Op[adjoint(d)]),
        (lowering_coefficient, Op[d]),
        (scalar, Op[]),
    )
end

"""
    Displace(a, α[, t])

Displace a Fock mode by the scalar amplitude `α`, so `a ↦ a + α`. The timed form
stores the complete c-number gauge of the moving displacement.
"""
function Displace(a::Op, α::Coefficient)
    d = fock_or_throw(a, "`Displace`")
    return fock_displacement(d, to_cnum(α))
end

function Displace(a::Op, α::Coefficient, t::Num)
    d = fock_or_throw(a, "`Displace`")
    tt = time_or_throw(t)
    c = to_cnum(α)
    derivative = dt(c, tt)
    derivative_adjoint = conj_cnum(derivative)
    gauge = fock_displacement_gauge(
        d, c, derivative, derivative_adjoint,
        mul_cnum(CNUM_NEG_IM, derivative),
        mul_cnum(CNUM_IM, derivative_adjoint),
    )
    return timed_transform(fock_displacement(d, c), gauge, tt)
end

"""
    Rotation(a, θ)
    Rotation(a, θ, t)

Rotate a Fock mode by `exp(-im*θ*a'a)`, so `a ↦ exp(-im*θ)*a`. The positional timed
form requires a symbolic moving angle and a symbolic time variable.
"""
function Rotation(a::Op, θ::Real)
    d = fock_or_throw(a, "`Rotation`")
    action = AffineAction(
        Op[d, adjoint(d)],
        CNum[
            conj_phase(θ) CNUM_ZERO
            CNUM_ZERO phase(θ)
        ],
        CNum[CNUM_ZERO, CNUM_ZERO],
    )
    return canonical_transform(action)
end

function Rotation(a::Op, θ::Num, t::Num)
    d = fock_or_throw(a, "`Rotation`")
    tt = time_or_throw(t)
    return timed_transform(Rotation(a, θ), gauge(adjoint(d) * d, θ, tt), tt)
end

"""
    Squeeze(a, r, ϕ = 0)
    Squeeze(a, r, ϕ, t)

Single-mode squeezing with `a ↦ cosh(r)*a + exp(im*ϕ)*sinh(r)*a'`. The timed form
supports both a moving magnitude and a moving phase.
"""
function Squeeze(a::Op, r::Real, ϕ::Real = 0)
    d = fock_or_throw(a, "`Squeeze`")
    ch = to_cnum(cosh(r))
    sh = mul_cnum(phase(ϕ), to_cnum(sinh(r)))
    action = AffineAction(
        Op[d, adjoint(d)],
        CNum[
            ch sh
            conj_cnum(sh) conj_cnum(ch)
        ],
        CNum[CNUM_ZERO, CNUM_ZERO];
        relations = ParamRelation[hyp_rel(r)],
    )
    return canonical_transform(action)
end

function squeeze_gauge(d::Op, r::Real, ϕ::Real, t::Num)
    sh = to_cnum(sinh(r))
    ch = to_cnum(cosh(r))
    imaginary_dr = mul_cnum(CNUM_IM, dt(r, t))
    dϕ = dt(ϕ, t)
    diagonal = mul_cnum(dϕ, mul_cnum(sh, sh))
    mixing = mul_cnum(dϕ, mul_cnum(sh, ch))
    return rule_qadd(
        (diagonal, Op[adjoint(d), d]),
        (mul_cnum(CNUM_HALF, diagonal), Op[]),
        (
            mul_cnum(
                CNUM_HALF,
                mul_cnum(phase(ϕ), add_cnum(mixing, neg_cnum(imaginary_dr))),
            ),
            Op[adjoint(d), adjoint(d)],
        ),
        (
            mul_cnum(
                CNUM_HALF,
                mul_cnum(conj_phase(ϕ), add_cnum(mixing, imaginary_dr)),
            ),
            Op[d, d],
        ),
    )
end

function Squeeze(a::Op, r::Real, ϕ::Real, t::Num)
    d = fock_or_throw(a, "`Squeeze`")
    tt = time_or_throw(t)
    return timed_transform(Squeeze(a, r, ϕ), squeeze_gauge(d, r, ϕ, tt), tt)
end

# === Two-mode and phase-space transformations ===

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
    return nothing
end

function beamsplitter(a::Op, b::Op, θ::Real)
    x, y = two_modes(a, b, "`BeamSplitter`")
    c = to_cnum(cos(θ))
    s = to_cnum(sin(θ))
    negative_s = neg_cnum(s)
    action = AffineAction(
        Op[x, y, adjoint(x), adjoint(y)],
        CNum[
            c s CNUM_ZERO CNUM_ZERO
            negative_s c CNUM_ZERO CNUM_ZERO
            CNUM_ZERO CNUM_ZERO conj_cnum(c) conj_cnum(s)
            CNUM_ZERO CNUM_ZERO neg_cnum(conj_cnum(s)) conj_cnum(c)
        ],
        CNum[CNUM_ZERO, CNUM_ZERO, CNUM_ZERO, CNUM_ZERO];
        relations = ParamRelation[trig_rel(θ)],
    )
    return canonical_transform(action)
end

"""Mix two Fock modes by a passive beam-splitter rotation."""
BeamSplitter(a::Op, b::Op, θ::Real) = beamsplitter(a, b, θ)

function BeamSplitter(a::Op, b::Op, θ::Real, t::Num)
    tt = time_or_throw(t)
    x, y = two_modes(a, b, "`BeamSplitter`")
    generator = im * (adjoint(x) * y - adjoint(y) * x)
    return timed_transform(beamsplitter(x, y, θ), gauge(generator, θ, tt), tt)
end

function quadrature_rotation(x::Op, p::Op, θ::Real)
    phase_pair(x, p, "`Rotation`")
    c = to_cnum(cos(θ))
    s = to_cnum(sin(θ))
    action = AffineAction(
        Op[x, p],
        CNum[
            c s
            neg_cnum(s) c
        ],
        CNum[CNUM_ZERO, CNUM_ZERO];
        relations = ParamRelation[trig_rel(θ)],
    )
    return canonical_transform(action)
end

"""Mix two Fock modes (beamsplitter) or rotate a canonical quadrature pair."""
Rotation(a::Op, b::Op, θ::Real) =
    is_phase_space(a) ? quadrature_rotation(a, b, θ) : BeamSplitter(a, b, θ)

function Rotation(a::Op, b::Op, θ::Real, t::Num)
    tt = time_or_throw(t)
    if is_phase_space(a)
        U = quadrature_rotation(a, b, θ)
        generator = (a * a + b * b) * (1 // 2)
        return timed_transform(U, gauge(generator, θ, tt), tt)
    end
    return BeamSplitter(a, b, θ, tt)
end

function two_mode_squeeze(a::Op, b::Op, r::Real)
    x, y = two_modes(a, b, "`TwoModeSqueeze`")
    u = to_cnum(cosh(r))
    v = to_cnum(sinh(r))
    action = AffineAction(
        Op[x, y, adjoint(x), adjoint(y)],
        CNum[
            u CNUM_ZERO CNUM_ZERO v
            CNUM_ZERO u v CNUM_ZERO
            CNUM_ZERO conj_cnum(v) conj_cnum(u) CNUM_ZERO
            conj_cnum(v) CNUM_ZERO CNUM_ZERO conj_cnum(u)
        ],
        CNum[CNUM_ZERO, CNUM_ZERO, CNUM_ZERO, CNUM_ZERO];
        relations = ParamRelation[hyp_rel(r)],
    )
    return canonical_transform(action)
end

"""Apply a two-mode bosonic squeezing transformation."""
TwoModeSqueeze(a::Op, b::Op, r::Real) = two_mode_squeeze(a, b, r)

function TwoModeSqueeze(a::Op, b::Op, r::Real, t::Num)
    tt = time_or_throw(t)
    x, y = two_modes(a, b, "`TwoModeSqueeze`")
    generator = im * (adjoint(x) * adjoint(y) - y * x)
    return timed_transform(two_mode_squeeze(x, y, r), gauge(generator, r, tt), tt)
end

function quadrature_squeeze(x::Op, p::Op, r::Real)
    phase_pair(x, p, "`Squeeze`")
    up = to_cnum(exp(r))
    down = to_cnum(inv(exp(r)))
    action = AffineAction(
        Op[x, p],
        CNum[
            up CNUM_ZERO
            CNUM_ZERO down
        ],
        CNum[CNUM_ZERO, CNUM_ZERO],
    )
    return canonical_transform(action)
end

"""Squeeze two Fock modes or a canonical quadrature pair."""
Squeeze(a::Op, b::Op, r::Real) =
    is_phase_space(a) ? quadrature_squeeze(a, b, r) : TwoModeSqueeze(a, b, r)

function Squeeze(a::Op, b::Op, r::Real, t::Num)
    tt = time_or_throw(t)
    if is_phase_space(a)
        U = quadrature_squeeze(a, b, r)
        generator = (a * b + b * a) * (1 // 2)
        return timed_transform(U, gauge(generator, r, tt), tt)
    end
    return TwoModeSqueeze(a, b, r, tt)
end

function quadrature_displacement(x::Op, p::Op, cx::CNum, cp::CNum)
    action = AffineAction(
        Op[x, p],
        CNum[
            CNUM_ONE CNUM_ZERO
            CNUM_ZERO CNUM_ONE
        ],
        CNum[cx, cp],
    )
    return canonical_transform(action)
end

function quadrature_displacement_gauge(
        x::Op, p::Op, cx::CNum, cp::CNum,
        derivative_x::CNum, derivative_p::CNum,
        x_coefficient::CNum, p_coefficient::CNum,
    )
    scalar = mul_cnum(
        neg_cnum(CNUM_HALF),
        add_cnum(
            mul_cnum(cp, derivative_x),
            neg_cnum(mul_cnum(cx, derivative_p)),
        ),
    )
    return rule_qadd(
        (x_coefficient, Op[x]), (p_coefficient, Op[p]), (scalar, Op[]),
    )
end

"""Displace a canonical quadrature pair by real scalar shifts."""
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
        derivative_p, neg_cnum(derivative_x),
    )
    return timed_transform(quadrature_displacement(x, p, cx, cp), gauge, tt)
end

# === Pauli and spin transformations ===

axis_or_throw(axis::Integer) = 1 <= axis <= 3 ? Int32(axis) :
    unitary_error("axis must be 1, 2, or 3, got $axis")
axis_op(o::Op, axis::Integer) =
    Op(o.kind, o.name_id, o.space_index, o.index, Int32(axis), 0, 0, 0)

function triple_or_throw(S::Op)
    (is_pauli(S) || is_spin(S)) || unitary_error(
        "`Rotation` expects a Pauli or Spin operator, got $(S.kind)",
    )
    return nothing
end

"""Rotate a Pauli or spin triple by `θ` around axis 1, 2, or 3."""
function Rotation(S::Op, axis::Integer, θ::Real)
    triple_or_throw(S)
    main_axis = Int(axis_or_throw(axis))
    u = axis_op(S, mod1(main_axis + 1, 3))
    v = axis_op(S, mod1(main_axis + 2, 3))
    fixed = axis_op(S, main_axis)
    c = to_cnum(cos(θ))
    s = to_cnum(sin(θ))
    action = AffineAction(
        Op[u, v, fixed],
        CNum[
            c neg_cnum(s) CNUM_ZERO
            s c CNUM_ZERO
            CNUM_ZERO CNUM_ZERO CNUM_ONE
        ],
        CNum[CNUM_ZERO, CNUM_ZERO, CNUM_ZERO];
        relations = ParamRelation[trig_rel(θ)],
    )
    return canonical_transform(action)
end

function Rotation(S::Op, axis::Integer, θ::Real, t::Num)
    triple_or_throw(S)
    tt = time_or_throw(t)
    fixed = axis_op(S, Int(axis_or_throw(axis)))
    generator = S.kind === OP_PAULI ? scaled(CNUM_HALF, fixed) :
        scaled(CNUM_ONE, fixed)
    return timed_transform(Rotation(S, axis, θ), gauge(generator, θ, tt), tt)
end

# === Ordinary N-level basis transformations ===

transition_op(o::Op, i::Integer, j::Integer) =
    Op(OP_TRANSITION, o.name_id, o.space_index, o.index, Int32(i), Int32(j), o.g, o.nlev)

function nlevel_or_throw(σ::Op, W::AbstractMatrix)
    is_transition(σ) || unitary_error(
        "`BasisRotation` expects an ordinary `Transition` operator, got $(σ.kind)",
    )
    n = Int(σ.nlev)
    size(W) == (n, n) || unitary_error(
        "`W` must be $n×$n for a $n-level space, got $(size(W))",
    )
    return n
end

function coefficient_matrix(W::AbstractMatrix)
    n, m = size(W)
    converted = Matrix{Coeff}(undef, n, m)
    for j in 1:m, i in 1:n
        converted[i, j] = to_cnum(W[i, j])
    end
    return converted
end

function exact_unitary_or_throw(W::Matrix{Coeff}, Wdagger::Matrix{Coeff})
    n = size(W, 1)
    scratch = ParamRelation[]
    for j in 1:n, k in 1:n
        residual = j == k ? CNUM_NEG1 : CNUM_ZERO
        for i in 1:n
            residual = add_cnum(
                residual, mul_cnum(Wdagger[j, i], W[i, k]),
            )
        end
        reduced = reduce_all(residual, ParamRelation[], false, scratch)
        iszero_cnum(reduced) || unitary_error(
            "`W` must be provably unitary exactly; entry ($j, $k) of `W'W - I` is " *
                "`$(to_num(reduced))`",
        )
    end
    return W
end

function matrix_unit_action(σ::Op, W::Matrix{Coeff}, Wdagger::Matrix{Coeff})
    n = size(W, 1)
    dimension = n * n
    basis = Op[transition_op(σ, i, j) for i in 1:n for j in 1:n]
    linear = fill(CNUM_ZERO, dimension, dimension)
    for i in 1:n, j in 1:n
        row = (i - 1) * n + j
        for k in 1:n, l in 1:n
            coefficient = mul_cnum(Wdagger[k, i], W[j, l])
            iszero_cnum(coefficient) && continue
            column = (k - 1) * n + l
            linear[row, column] = coefficient
        end
    end
    return AffineAction(
        UnitaryLinearAction(), basis, linear, fill(CNUM_ZERO, dimension),
    )
end

function dagger_coefficients(W::Matrix{Coeff})
    n = size(W, 1)
    out = Matrix{Coeff}(undef, n, n)
    for j in 1:n, i in 1:n
        out[i, j] = conj_cnum(W[j, i])
    end
    return out
end

function nlevel_rotation(σ::Op, W::AbstractMatrix)
    nlevel_or_throw(σ, W)
    coefficients = coefficient_matrix(W)
    dagger = dagger_coefficients(coefficients)
    exact_unitary_or_throw(coefficients, dagger)
    U = canonical_transform(matrix_unit_action(σ, coefficients, dagger))
    return U, coefficients
end

"""
    BasisRotation(σ, W)
    BasisRotation(σ, W, t)

Rotate an ordinary N-level basis by a matrix `W` satisfying `W'W = I`. The timed form
derives the Hamiltonian gauge `im*Ẇ'W` entrywise with respect to `t`.
"""
function BasisRotation(σ::Op, W::AbstractMatrix)
    U, _ = nlevel_rotation(σ, W)
    return U
end

function nlevel_gauge(σ::Op, W::Matrix{Coeff}, t::Num)
    n = size(W, 1)
    conjugated_derivatives = Matrix{Coeff}(undef, n, n)
    for j in 1:n, k in 1:n
        conjugated_derivatives[k, j] = conj_cnum(dt(W[k, j], t))
    end
    gauge = QTermDict()
    for j in 1:n, l in 1:n
        coefficient = CNUM_ZERO
        for k in 1:n
            coefficient = add_cnum(
                coefficient, mul_cnum(conjugated_derivatives[k, j], W[k, l]),
            )
        end
        coefficient = mul_cnum(CNUM_IM, coefficient)
        iszero_cnum(coefficient) ||
            addto!(gauge, Op[transition_op(σ, j, l)], coefficient)
    end
    return QAdd(gauge, EMPTY_INDICES)
end

function BasisRotation(σ::Op, W::AbstractMatrix, t::Num)
    tt = time_or_throw(t)
    U, coefficients = nlevel_rotation(σ, W)
    return timed_transform(U, nlevel_gauge(σ, coefficients, tt), tt)
end

Rotation(σ::Op, W::AbstractMatrix) = BasisRotation(σ, W)
Rotation(σ::Op, W::AbstractMatrix, t::Num) = BasisRotation(σ, W, t)

include("unitary_displacement_frames.jl")
