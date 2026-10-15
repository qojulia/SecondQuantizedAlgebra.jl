function positive_frequency_or_throw(ω::Real, what::AbstractString)
    c = to_cnum(ω)
    if is_native(c)
        z = native_scalar(c)
        (iszero(imag(z)) && real(z) > 0) || throw(
            ArgumentError("$what needs a positive frequency, got $ω"),
        )
    end
    return ω
end

frequency_root(x::Real) = sqrt(Num(x))

"""
    quadratures(a, ω = 1) -> (x, p)

Position and momentum of the oscillator whose ladder operator at frequency `ω` is the Fock
mode `a`, with unit mass and ``\\hbar = 1``:

```math
x = \\frac{a + a^\\dagger}{\\sqrt{2\\omega}}, \\qquad
p = \\frac{i\\omega(a^\\dagger - a)}{\\sqrt{2\\omega}}, \\qquad [x, p] = i.
```

Writing a Hamiltonian in `x` and `p` and choosing the frequency here selects its ladder
basis: `p^2/2 + ω^2*x^2/2` is `ω*(a'*a + 1/2)` at the oscillator's own frequency and
acquires `a^2 + a'^2` terms at any other reference frequency. Integer and rational `ω`
give exact radicals; a symbolic `ω` is assumed positive. [`ladder_operator`](@ref) is
the inverse map.

# Examples

```jldoctest
julia> @qnumbers a::Destroy(FockSpace(:a));

julia> x, p = SecondQuantizedAlgebra.quadratures(a, 2);

julia> x
1//2 * a + 1//2 * a'

julia> simplify(p * p / 2 + 2^2 * x * x / 2)
1 + 2 * a' * a
```

See also [`ladder_operator`](@ref), [`FrequencyChange`](@ref), [`Position`](@ref).
"""
function quadratures(a::Op, ω::Real = 1)
    d = fock_or_throw(a, "`quadratures`")
    positive_frequency_or_throw(ω, "`quadratures`")
    root = frequency_root(2 * ω)
    x = (d + adjoint(d)) / root
    p = (im * ω) * (adjoint(d) - d) / root
    return x, p
end

"""
    ladder_operator(x, p, ω = 1) -> QField

Ladder operator at frequency `ω` of a canonical pair ``[x, p] = i``, with unit mass and
``\\hbar = 1``:

```math
a_\\omega = \\frac{\\omega x + i p}{\\sqrt{2\\omega}}.
```

`x` and `p` may be [`Position`](@ref) and [`Momentum`](@ref) operators or the quadratures
of a Fock mode from [`quadratures`](@ref). The second case expresses one ladder basis of
an oscillator in another: `ladder_operator(quadratures(a, ω0)..., ωr)` is
``u a - v a^\\dagger`` with ``u = (\\omega_0 + \\omega_r)/(2\\sqrt{\\omega_0\\omega_r})``,
``v = (\\omega_0 - \\omega_r)/(2\\sqrt{\\omega_0\\omega_r})``. The canonical pair is the
caller's contract; it is not checked. The result is simplified.

# Examples

```jldoctest
julia> @qnumbers a::Destroy(FockSpace(:a));

julia> SecondQuantizedAlgebra.ladder_operator(SecondQuantizedAlgebra.quadratures(a, 2)..., 1//2)
5//4 * a - 3//4 * a'
```

See also [`quadratures`](@ref), [`FrequencyChange`](@ref).
"""
function ladder_operator(x::QField, p::QField, ω::Real = 1)
    positive_frequency_or_throw(ω, "`ladder_operator`")
    return simplify((ω * x + im * p) / frequency_root(2 * ω))
end

"""
    FrequencyChange(a, ω0, ωr)

Change the ladder basis of one oscillator from frequency `ω0` to the reference frequency
`ωr`. The Fock mode `a` is the ladder operator at `ω0` before the transformation and at
`ωr` after it, with the shared quadratures of [`quadratures`](@ref):

```math
a \\mapsto u a + v a^\\dagger, \\qquad
u = \\frac{\\omega_0 + \\omega_r}{2\\sqrt{\\omega_0\\omega_r}}, \\quad
v = \\frac{\\omega_0 - \\omega_r}{2\\sqrt{\\omega_0\\omega_r}}.
```

This is [`Squeeze`](@ref)`(a, log(ω0/ωr)/2)` with exact algebraic coefficients in place
of `cosh` and `sinh` of a logarithm. Integer and rational frequencies give exact
radicals; symbolic frequencies are assumed positive.

# Examples

```jldoctest
julia> @qnumbers a::Destroy(FockSpace(:a));

julia> conjugate(a, SecondQuantizedAlgebra.FrequencyChange(a, 2, 1//2))
5//4 * a + 3//4 * a'
```

See also [`quadratures`](@ref), [`ladder_operator`](@ref), [`Bogoliubov`](@ref).
"""
function FrequencyChange(a::Op, ω0::Real, ωr::Real)
    d = fock_or_throw(a, "`FrequencyChange`")
    positive_frequency_or_throw(ω0, "`FrequencyChange`")
    positive_frequency_or_throw(ωr, "`FrequencyChange`")
    root = frequency_root(2 * ω0) * frequency_root(2 * ωr)
    u = to_cnum((ω0 + ωr) / root)
    v = to_cnum((ω0 - ωr) / root)
    return exact_bogoliubov(Op[d], CNum[u v; conj_cnum(v) conj_cnum(u)])
end
