# === Fock transformations ===

"""
    Displace(a, α[, t])

Displace a Fock mode by the scalar amplitude `α`, so `a ↦ a + α`. The timed form
stores the complete c-number gauge of the moving displacement.
"""
function Displace(a::Op, α::Coefficient)
    d = _fock_or_throw(a, "`Displace`")
    c = _to_cnum(α)
    return _static_transform(
        _with_adjoint(d, _rule_qadd((_CNUM_ONE, Op[d]), (c, Op[]))),
        _with_adjoint(d, _rule_qadd((_CNUM_ONE, Op[d]), (_neg_cnum(c), Op[]))),
    )
end

function Displace(a::Op, α::Coefficient, t::Num)
    d = _fock_or_throw(a, "`Displace`")
    tt = _time_or_throw(t)
    c = _to_cnum(α)
    derivative = _dt(c, tt)
    scalar = _mul_cnum(
        _CNUM_NEG_IM,
        _mul_cnum(
            _CNUM_HALF,
            _add_cnum(
                _mul_cnum(_conj_cnum(c), derivative),
                _neg_cnum(_mul_cnum(c, _conj_cnum(derivative))),
            ),
        ),
    )
    gauge = _rule_qadd(
        (_mul_cnum(_CNUM_NEG_IM, derivative), Op[adjoint(d)]),
        (_mul_cnum(_CNUM_IM, _conj_cnum(derivative)), Op[d]),
        (scalar, Op[]),
    )
    return _timed_transform(Displace(a, α), gauge, tt)
end

"""
    Rotation(a, θ)
    Rotation(a, θ, t)

Rotate a Fock mode by `exp(-im*θ*a'a)`, so `a ↦ exp(-im*θ)*a`. The positional timed
form requires a symbolic moving angle and a symbolic time variable.
"""
function Rotation(a::Op, θ::Real)
    d = _fock_or_throw(a, "`Rotation`")
    return _static_transform(
        _with_adjoint(d, _scaled(_conj_phase(θ), d)),
        _with_adjoint(d, _scaled(_phase(θ), d)),
    )
end

function Rotation(a::Op, θ::Num, t::Num)
    d = _fock_or_throw(a, "`Rotation`")
    tt = _time_or_throw(t)
    return _timed_transform(Rotation(a, θ), _gauge(adjoint(d) * d, θ, tt), tt)
end

"""
    Squeeze(a, r, ϕ = 0)
    Squeeze(a, r, ϕ, t)

Single-mode squeezing with `a ↦ cosh(r)*a + exp(im*ϕ)*sinh(r)*a'`. The timed form
supports both a moving magnitude and a moving phase.
"""
function Squeeze(a::Op, r::Real, ϕ::Real = 0)
    d = _fock_or_throw(a, "`Squeeze`")
    ch = _to_cnum(cosh(r))
    sh = _mul_cnum(_phase(ϕ), _to_cnum(sinh(r)))
    return _static_transform(
        _with_adjoint(d, _rule_qadd((ch, Op[d]), (sh, Op[adjoint(d)]))),
        _with_adjoint(d, _rule_qadd((ch, Op[d]), (_neg_cnum(sh), Op[adjoint(d)]))),
        ParamRelation[_hyp_rel(r)],
    )
end

function _squeeze_gauge(d::Op, r::Real, ϕ::Real, t::Num)
    sh = _to_cnum(sinh(r))
    ch = _to_cnum(cosh(r))
    imaginary_dr = _mul_cnum(_CNUM_IM, _dt(r, t))
    dϕ = _dt(ϕ, t)
    diagonal = _mul_cnum(dϕ, _mul_cnum(sh, sh))
    mixing = _mul_cnum(dϕ, _mul_cnum(sh, ch))
    return _rule_qadd(
        (diagonal, Op[adjoint(d), d]),
        (_mul_cnum(_CNUM_HALF, diagonal), Op[]),
        (
            _mul_cnum(
                _CNUM_HALF,
                _mul_cnum(_phase(ϕ), _add_cnum(mixing, _neg_cnum(imaginary_dr))),
            ),
            Op[adjoint(d), adjoint(d)],
        ),
        (
            _mul_cnum(
                _CNUM_HALF,
                _mul_cnum(_conj_phase(ϕ), _add_cnum(mixing, imaginary_dr)),
            ),
            Op[d, d],
        ),
    )
end

function Squeeze(a::Op, r::Real, ϕ::Real, t::Num)
    d = _fock_or_throw(a, "`Squeeze`")
    tt = _time_or_throw(t)
    return _timed_transform(Squeeze(a, r, ϕ), _squeeze_gauge(d, r, ϕ, tt), tt)
end

# === Two-mode and phase-space transformations ===

function _two_modes(a::Op, b::Op, what::AbstractString)
    x = _fock_or_throw(a, what)
    y = _fock_or_throw(b, what)
    site_key(x) == site_key(y) && _unitary_error("$what needs two distinct modes")
    return x, y
end

function _phase_pair(x::Op, p::Op, what::AbstractString)
    (is_position(x) && is_momentum(p)) || _unitary_error(
        "$what expects a `(Position, Momentum)` pair in that order",
    )
    site_key(x) == site_key(p) || _unitary_error(
        "$what expects both quadratures on the same site",
    )
    return nothing
end

function _beamsplitter(a::Op, b::Op, θ::Real)
    x, y = _two_modes(a, b, "`Rotation`")
    c = _to_cnum(cos(θ))
    s = _to_cnum(sin(θ))
    negative_s = _neg_cnum(s)
    return _static_transform(
        merge(
            _with_adjoint(x, _rule_qadd((c, Op[x]), (s, Op[y]))),
            _with_adjoint(y, _rule_qadd((c, Op[y]), (negative_s, Op[x]))),
        ),
        merge(
            _with_adjoint(x, _rule_qadd((c, Op[x]), (negative_s, Op[y]))),
            _with_adjoint(y, _rule_qadd((c, Op[y]), (s, Op[x]))),
        ),
        ParamRelation[_trig_rel(θ)],
    )
end

function _quadrature_rotation(x::Op, p::Op, θ::Real)
    _phase_pair(x, p, "`Rotation`")
    c = _to_cnum(cos(θ))
    s = _to_cnum(sin(θ))
    negative_s = _neg_cnum(s)
    return _static_transform(
        _pair_rules(
            x, p, _rule_qadd((c, Op[x]), (s, Op[p])),
            _rule_qadd((c, Op[p]), (negative_s, Op[x])),
        ),
        _pair_rules(
            x, p, _rule_qadd((c, Op[x]), (negative_s, Op[p])),
            _rule_qadd((c, Op[p]), (s, Op[x])),
        ),
        ParamRelation[_trig_rel(θ)],
    )
end

"""Mix two Fock modes (beamsplitter) or rotate a canonical quadrature pair."""
Rotation(a::Op, b::Op, θ::Real) =
    _is_phase_space(a) ? _quadrature_rotation(a, b, θ) : _beamsplitter(a, b, θ)

function Rotation(a::Op, b::Op, θ::Real, t::Num)
    tt = _time_or_throw(t)
    if _is_phase_space(a)
        U = _quadrature_rotation(a, b, θ)
        generator = (a * a + b * b) * (1 // 2)
    else
        x, y = _two_modes(a, b, "`Rotation`")
        U = _beamsplitter(x, y, θ)
        generator = im * (adjoint(x) * y - adjoint(y) * x)
    end
    return _timed_transform(U, _gauge(generator, θ, tt), tt)
end

function _two_mode_squeeze(a::Op, b::Op, r::Real)
    x, y = _two_modes(a, b, "`Squeeze`")
    u = _to_cnum(cosh(r))
    v = _to_cnum(sinh(r))
    negative_v = _neg_cnum(v)
    return _static_transform(
        merge(
            _with_adjoint(x, _rule_qadd((u, Op[x]), (v, Op[adjoint(y)]))),
            _with_adjoint(y, _rule_qadd((u, Op[y]), (v, Op[adjoint(x)]))),
        ),
        merge(
            _with_adjoint(x, _rule_qadd((u, Op[x]), (negative_v, Op[adjoint(y)]))),
            _with_adjoint(y, _rule_qadd((u, Op[y]), (negative_v, Op[adjoint(x)]))),
        ),
        ParamRelation[_hyp_rel(r)],
    )
end

function _quadrature_squeeze(x::Op, p::Op, r::Real)
    _phase_pair(x, p, "`Squeeze`")
    up = _to_cnum(exp(r))
    down = _to_cnum(inv(exp(r)))
    return _static_transform(
        _pair_rules(x, p, _scaled(up, x), _scaled(down, p)),
        _pair_rules(x, p, _scaled(down, x), _scaled(up, p)),
    )
end

"""Squeeze two Fock modes or a canonical quadrature pair."""
Squeeze(a::Op, b::Op, r::Real) =
    _is_phase_space(a) ? _quadrature_squeeze(a, b, r) : _two_mode_squeeze(a, b, r)

function Squeeze(a::Op, b::Op, r::Real, t::Num)
    tt = _time_or_throw(t)
    if _is_phase_space(a)
        U = _quadrature_squeeze(a, b, r)
        generator = (a * b + b * a) * (1 // 2)
    else
        x, y = _two_modes(a, b, "`Squeeze`")
        U = _two_mode_squeeze(x, y, r)
        generator = im * (adjoint(x) * adjoint(y) - y * x)
    end
    return _timed_transform(U, _gauge(generator, r, tt), tt)
end

"""Displace a canonical quadrature pair by real scalar shifts."""
function Displace(x::Op, p::Op, dx::Real, dp::Real)
    _phase_pair(x, p, "`Displace`")
    cx = _to_cnum(dx)
    cp = _to_cnum(dp)
    return _static_transform(
        _pair_rules(
            x, p, _rule_qadd((_CNUM_ONE, Op[x]), (cx, Op[])),
            _rule_qadd((_CNUM_ONE, Op[p]), (cp, Op[])),
        ),
        _pair_rules(
            x, p, _rule_qadd((_CNUM_ONE, Op[x]), (_neg_cnum(cx), Op[])),
            _rule_qadd((_CNUM_ONE, Op[p]), (_neg_cnum(cp), Op[])),
        ),
    )
end

function Displace(x::Op, p::Op, dx::Real, dp::Real, t::Num)
    tt = _time_or_throw(t)
    cx = _to_cnum(dx)
    cp = _to_cnum(dp)
    derivative_x = _dt(cx, tt)
    derivative_p = _dt(cp, tt)
    scalar = _mul_cnum(
        _neg_cnum(_CNUM_HALF),
        _add_cnum(
            _mul_cnum(cp, derivative_x),
            _neg_cnum(_mul_cnum(cx, derivative_p)),
        ),
    )
    gauge = _rule_qadd(
        (derivative_p, Op[x]), (_neg_cnum(derivative_x), Op[p]), (scalar, Op[]),
    )
    return _timed_transform(Displace(x, p, dx, dp), gauge, tt)
end

# === Pauli and spin transformations ===

_axis_or_throw(axis::Integer) = 1 <= axis <= 3 ? Int32(axis) :
    _unitary_error("axis must be 1, 2, or 3, got $axis")
_axis_op(o::Op, axis::Integer) =
    Op(o.kind, o.name_id, o.space_index, o.index, Int32(axis), 0, 0, 0)

function _triple_or_throw(S::Op)
    (is_pauli(S) || is_spin(S)) || _unitary_error(
        "`Rotation` expects a Pauli or Spin operator, got $(S.kind)",
    )
    return nothing
end

"""Rotate a Pauli or spin triple by `θ` around axis 1, 2, or 3."""
function Rotation(S::Op, axis::Integer, θ::Real)
    _triple_or_throw(S)
    main_axis = Int(_axis_or_throw(axis))
    u = _axis_op(S, mod1(main_axis + 1, 3))
    v = _axis_op(S, mod1(main_axis + 2, 3))
    fixed = _axis_op(S, main_axis)
    c = _to_cnum(cos(θ))
    s = _to_cnum(sin(θ))
    negative_s = _neg_cnum(s)
    return _static_transform(
        Dict{Op, QAdd}(
            u => _rule_qadd((c, Op[u]), (negative_s, Op[v])),
            v => _rule_qadd((c, Op[v]), (s, Op[u])),
            fixed => _scaled(_CNUM_ONE, fixed),
        ),
        Dict{Op, QAdd}(
            u => _rule_qadd((c, Op[u]), (s, Op[v])),
            v => _rule_qadd((c, Op[v]), (negative_s, Op[u])),
            fixed => _scaled(_CNUM_ONE, fixed),
        ),
        ParamRelation[_trig_rel(θ)],
    )
end

function Rotation(S::Op, axis::Integer, θ::Real, t::Num)
    _triple_or_throw(S)
    tt = _time_or_throw(t)
    fixed = _axis_op(S, Int(_axis_or_throw(axis)))
    generator = S.kind === OP_PAULI ? _scaled(_CNUM_HALF, fixed) :
        _scaled(_CNUM_ONE, fixed)
    return _timed_transform(Rotation(S, axis, θ), _gauge(generator, θ, tt), tt)
end

# === Ordinary N-level basis transformations ===

_transition_op(o::Op, i::Integer, j::Integer) =
    Op(OP_TRANSITION, o.name_id, o.space_index, o.index, Int32(i), Int32(j), o.g, o.nlev)

function _nlevel_or_throw(σ::Op, W::AbstractMatrix)
    is_transition(σ) || _unitary_error(
        "`Rotation` expects an ordinary `Transition` operator, got $(σ.kind)",
    )
    n = Int(σ.nlev)
    size(W) == (n, n) || _unitary_error(
        "`W` must be $n×$n for a $n-level space, got $(size(W))",
    )
    return n
end

function _coefficient_matrix(W::AbstractMatrix)
    n, m = size(W)
    converted = Matrix{Coeff}(undef, n, m)
    for j in 1:m, i in 1:n
        converted[i, j] = _to_cnum(W[i, j])
    end
    return converted
end

function _exact_unitary_or_throw(W::Matrix{Coeff}, Wdagger::Matrix{Coeff})
    n = size(W, 1)
    scratch = ParamRelation[]
    for j in 1:n, k in 1:n
        residual = j == k ? _CNUM_NEG1 : _CNUM_ZERO
        for i in 1:n
            residual = _add_cnum(
                residual, _mul_cnum(Wdagger[j, i], W[i, k]),
            )
        end
        reduced = _reduce_all(residual, ParamRelation[], false, scratch)
        _iszero_cnum(reduced) || _unitary_error(
            "`W` must be provably unitary exactly; entry ($j, $k) of `W'W - I` is " *
                "`$(to_num(reduced))`",
        )
    end
    return W
end

function _matrix_unit_rules(
        σ::Op, W::Matrix{Coeff}, Wdagger::Matrix{Coeff} = _dagger_coefficients(W),
    )
    n = size(W, 1)
    operator_terms = QTerm[
        QTerm(Op[_transition_op(σ, k, l)], _EMPTY_NE) for k in 1:n, l in 1:n
    ]
    nonzero_rows = [
        [(column, W[row, column]) for column in 1:n if !_iszero_cnum(W[row, column])]
        for row in 1:n
    ]
    rules = Dict{Op, QAdd}()
    for i in 1:n, j in 1:n
        terms = QTermDict()
        for (k, _) in nonzero_rows[i], (l, wjl) in nonzero_rows[j]
            coefficient = _mul_cnum(Wdagger[k, i], wjl)
            _iszero_cnum(coefficient) ||
                _addto_key!(terms, operator_terms[k, l], coefficient)
        end
        rules[_transition_op(σ, i, j)] = QAdd(terms, _EMPTY_INDICES)
    end
    return rules
end

function _dagger_coefficients(W::Matrix{Coeff})
    n = size(W, 1)
    out = Matrix{Coeff}(undef, n, n)
    for j in 1:n, i in 1:n
        out[i, j] = _conj_cnum(W[j, i])
    end
    return out
end

function _nlevel_rotation(σ::Op, W::AbstractMatrix)
    _nlevel_or_throw(σ, W)
    coefficients = _coefficient_matrix(W)
    dagger = _dagger_coefficients(coefficients)
    _exact_unitary_or_throw(coefficients, dagger)
    U = _static_transform(
        _matrix_unit_rules(σ, coefficients, dagger),
        _matrix_unit_rules(σ, dagger, coefficients),
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
    U, _ = _nlevel_rotation(σ, W)
    return U
end

function _nlevel_gauge(σ::Op, W::Matrix{Coeff}, t::Num)
    n = size(W, 1)
    conjugated_derivatives = Matrix{Coeff}(undef, n, n)
    for j in 1:n, k in 1:n
        conjugated_derivatives[k, j] = _conj_cnum(_dt(W[k, j], t))
    end
    gauge = QTermDict()
    for j in 1:n, l in 1:n
        coefficient = _CNUM_ZERO
        for k in 1:n
            coefficient = _add_cnum(
                coefficient, _mul_cnum(conjugated_derivatives[k, j], W[k, l]),
            )
        end
        coefficient = _mul_cnum(_CNUM_IM, coefficient)
        _iszero_cnum(coefficient) ||
            _addto!(gauge, Op[_transition_op(σ, j, l)], coefficient)
    end
    return QAdd(gauge, _EMPTY_INDICES)
end

function Rotation(σ::Op, W::AbstractMatrix, t::Num)
    tt = _time_or_throw(t)
    U, coefficients = _nlevel_rotation(σ, W)
    return _timed_transform(U, _nlevel_gauge(σ, coefficients, tt), tt)
end
