# === Fock transformations ===

"""
    Displace(a, α[, t])

Displace a Fock mode by the scalar amplitude `α`, so `a ↦ a + α`. The timed form
stores the complete c-number gauge of the moving displacement.
"""
function Displace(a::Op, α::Coefficient)
    d = fock_or_throw(a, "`Displace`")
    c = to_cnum(α)
    return static_transform(
        with_adjoint(d, rule_qadd((CNUM_ONE, Op[d]), (c, Op[]))),
        with_adjoint(d, rule_qadd((CNUM_ONE, Op[d]), (neg_cnum(c), Op[]))),
    )
end

function Displace(a::Op, α::Coefficient, t::Num)
    d = fock_or_throw(a, "`Displace`")
    tt = time_or_throw(t)
    c = to_cnum(α)
    derivative = dt(c, tt)
    scalar = mul_cnum(
        CNUM_NEG_IM,
        mul_cnum(
            CNUM_HALF,
            add_cnum(
                mul_cnum(conj_cnum(c), derivative),
                neg_cnum(mul_cnum(c, conj_cnum(derivative))),
            ),
        ),
    )
    gauge = rule_qadd(
        (mul_cnum(CNUM_NEG_IM, derivative), Op[adjoint(d)]),
        (mul_cnum(CNUM_IM, conj_cnum(derivative)), Op[d]),
        (scalar, Op[]),
    )
    return timed_transform(Displace(a, α), gauge, tt)
end

"""
    Rotation(a, θ)
    Rotation(a, θ, t)

Rotate a Fock mode by `exp(-im*θ*a'a)`, so `a ↦ exp(-im*θ)*a`. The positional timed
form requires a symbolic moving angle and a symbolic time variable.
"""
function Rotation(a::Op, θ::Real)
    d = fock_or_throw(a, "`Rotation`")
    return static_transform(
        with_adjoint(d, scaled(conj_phase(θ), d)),
        with_adjoint(d, scaled(phase(θ), d)),
    )
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
    return static_transform(
        with_adjoint(d, rule_qadd((ch, Op[d]), (sh, Op[adjoint(d)]))),
        with_adjoint(d, rule_qadd((ch, Op[d]), (neg_cnum(sh), Op[adjoint(d)]))),
        ParamRelation[hyp_rel(r)],
    )
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
    x, y = two_modes(a, b, "`Rotation`")
    c = to_cnum(cos(θ))
    s = to_cnum(sin(θ))
    negative_s = neg_cnum(s)
    return static_transform(
        merge(
            with_adjoint(x, rule_qadd((c, Op[x]), (s, Op[y]))),
            with_adjoint(y, rule_qadd((c, Op[y]), (negative_s, Op[x]))),
        ),
        merge(
            with_adjoint(x, rule_qadd((c, Op[x]), (negative_s, Op[y]))),
            with_adjoint(y, rule_qadd((c, Op[y]), (s, Op[x]))),
        ),
        ParamRelation[trig_rel(θ)],
    )
end

function quadrature_rotation(x::Op, p::Op, θ::Real)
    phase_pair(x, p, "`Rotation`")
    c = to_cnum(cos(θ))
    s = to_cnum(sin(θ))
    negative_s = neg_cnum(s)
    return static_transform(
        pair_rules(
            x, p, rule_qadd((c, Op[x]), (s, Op[p])),
            rule_qadd((c, Op[p]), (negative_s, Op[x])),
        ),
        pair_rules(
            x, p, rule_qadd((c, Op[x]), (negative_s, Op[p])),
            rule_qadd((c, Op[p]), (s, Op[x])),
        ),
        ParamRelation[trig_rel(θ)],
    )
end

"""Mix two Fock modes (beamsplitter) or rotate a canonical quadrature pair."""
Rotation(a::Op, b::Op, θ::Real) =
    is_phase_space(a) ? quadrature_rotation(a, b, θ) : beamsplitter(a, b, θ)

function Rotation(a::Op, b::Op, θ::Real, t::Num)
    tt = time_or_throw(t)
    if is_phase_space(a)
        U = quadrature_rotation(a, b, θ)
        generator = (a * a + b * b) * (1 // 2)
    else
        x, y = two_modes(a, b, "`Rotation`")
        U = beamsplitter(x, y, θ)
        generator = im * (adjoint(x) * y - adjoint(y) * x)
    end
    return timed_transform(U, gauge(generator, θ, tt), tt)
end

function two_mode_squeeze(a::Op, b::Op, r::Real)
    x, y = two_modes(a, b, "`Squeeze`")
    u = to_cnum(cosh(r))
    v = to_cnum(sinh(r))
    negative_v = neg_cnum(v)
    return static_transform(
        merge(
            with_adjoint(x, rule_qadd((u, Op[x]), (v, Op[adjoint(y)]))),
            with_adjoint(y, rule_qadd((u, Op[y]), (v, Op[adjoint(x)]))),
        ),
        merge(
            with_adjoint(x, rule_qadd((u, Op[x]), (negative_v, Op[adjoint(y)]))),
            with_adjoint(y, rule_qadd((u, Op[y]), (negative_v, Op[adjoint(x)]))),
        ),
        ParamRelation[hyp_rel(r)],
    )
end

function quadrature_squeeze(x::Op, p::Op, r::Real)
    phase_pair(x, p, "`Squeeze`")
    up = to_cnum(exp(r))
    down = to_cnum(inv(exp(r)))
    return static_transform(
        pair_rules(x, p, scaled(up, x), scaled(down, p)),
        pair_rules(x, p, scaled(down, x), scaled(up, p)),
    )
end

"""Squeeze two Fock modes or a canonical quadrature pair."""
Squeeze(a::Op, b::Op, r::Real) =
    is_phase_space(a) ? quadrature_squeeze(a, b, r) : two_mode_squeeze(a, b, r)

function Squeeze(a::Op, b::Op, r::Real, t::Num)
    tt = time_or_throw(t)
    if is_phase_space(a)
        U = quadrature_squeeze(a, b, r)
        generator = (a * b + b * a) * (1 // 2)
    else
        x, y = two_modes(a, b, "`Squeeze`")
        U = two_mode_squeeze(x, y, r)
        generator = im * (adjoint(x) * adjoint(y) - y * x)
    end
    return timed_transform(U, gauge(generator, r, tt), tt)
end

"""Displace a canonical quadrature pair by real scalar shifts."""
function Displace(x::Op, p::Op, dx::Real, dp::Real)
    phase_pair(x, p, "`Displace`")
    cx = to_cnum(dx)
    cp = to_cnum(dp)
    return static_transform(
        pair_rules(
            x, p, rule_qadd((CNUM_ONE, Op[x]), (cx, Op[])),
            rule_qadd((CNUM_ONE, Op[p]), (cp, Op[])),
        ),
        pair_rules(
            x, p, rule_qadd((CNUM_ONE, Op[x]), (neg_cnum(cx), Op[])),
            rule_qadd((CNUM_ONE, Op[p]), (neg_cnum(cp), Op[])),
        ),
    )
end

function Displace(x::Op, p::Op, dx::Real, dp::Real, t::Num)
    tt = time_or_throw(t)
    cx = to_cnum(dx)
    cp = to_cnum(dp)
    derivative_x = dt(cx, tt)
    derivative_p = dt(cp, tt)
    scalar = mul_cnum(
        neg_cnum(CNUM_HALF),
        add_cnum(
            mul_cnum(cp, derivative_x),
            neg_cnum(mul_cnum(cx, derivative_p)),
        ),
    )
    gauge = rule_qadd(
        (derivative_p, Op[x]), (neg_cnum(derivative_x), Op[p]), (scalar, Op[]),
    )
    return timed_transform(Displace(x, p, dx, dp), gauge, tt)
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
    negative_s = neg_cnum(s)
    return static_transform(
        Dict{Op, QAdd}(
            u => rule_qadd((c, Op[u]), (negative_s, Op[v])),
            v => rule_qadd((c, Op[v]), (s, Op[u])),
            fixed => scaled(CNUM_ONE, fixed),
        ),
        Dict{Op, QAdd}(
            u => rule_qadd((c, Op[u]), (s, Op[v])),
            v => rule_qadd((c, Op[v]), (negative_s, Op[u])),
            fixed => scaled(CNUM_ONE, fixed),
        ),
        ParamRelation[trig_rel(θ)],
    )
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
        "`Rotation` expects an ordinary `Transition` operator, got $(σ.kind)",
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

function matrix_unit_rules(
        σ::Op, W::Matrix{Coeff}, Wdagger::Matrix{Coeff} = dagger_coefficients(W),
    )
    n = size(W, 1)
    operator_terms = QTerm[
        QTerm(Op[transition_op(σ, k, l)], EMPTY_NE) for k in 1:n, l in 1:n
    ]
    nonzero_rows = [
        [(column, W[row, column]) for column in 1:n if !iszero_cnum(W[row, column])]
            for row in 1:n
    ]
    rules = Dict{Op, QAdd}()
    for i in 1:n, j in 1:n
        terms = QTermDict()
        for (k, _) in nonzero_rows[i], (l, wjl) in nonzero_rows[j]
            coefficient = mul_cnum(Wdagger[k, i], wjl)
            iszero_cnum(coefficient) ||
                addto_key!(terms, operator_terms[k, l], coefficient)
        end
        rules[transition_op(σ, i, j)] = QAdd(terms, EMPTY_INDICES)
    end
    return rules
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
    U = static_transform(
        matrix_unit_rules(σ, coefficients, dagger),
        matrix_unit_rules(σ, dagger, coefficients),
    )
    return U, coefficients
end

"""
    Rotation(σ, W)
    Rotation(σ, W, t)

Rotate an ordinary N-level basis by a matrix `W` satisfying `W'W = I`. The timed form
derives the Hamiltonian gauge `im*Ẇ'W` entrywise with respect to `t`.
"""
function Rotation(σ::Op, W::AbstractMatrix)
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

function Rotation(σ::Op, W::AbstractMatrix, t::Num)
    tt = time_or_throw(t)
    U, coefficients = nlevel_rotation(σ, W)
    return timed_transform(U, nlevel_gauge(σ, coefficients, tt), tt)
end
