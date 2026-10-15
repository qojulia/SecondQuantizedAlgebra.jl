function harmonics(S::SuperOperator, ω::Real, t::Num)
    frequency = harmonic_frequency_coeff(ω, t)
    parts = Dict{Int, SuperTermDict}()
    for (term, c) in S.arguments, (n, cn) in coefficient_harmonics(c, Num(ω), frequency, t)
        addto_key!(get!(SuperTermDict, parts, n), copy_key(term), cn)
    end
    return Dict{Int, SuperOperator}(n => superop_pruned(d, copy(S.indices)) for (n, d) in parts)
end

function harmonic(S::SuperOperator, ω::Real, t::Num, n::Integer)
    frequency = harmonic_frequency_coeff(ω, t)
    order = Int(n)
    return map_coefficients(c -> get(coefficient_harmonics(c, Num(ω), frequency, t), order, CNUM_ZERO), S)
end

function harmonics(L::Lindbladian{QAdd}, ω::Real, t::Num)
    builders = Dict{Int, typeof(LindbladianBuilder{QAdd}())}()
    builder(n) = get!(LindbladianBuilder{QAdd}, builders, n)
    for (n, Hn) in harmonics(L.hamiltonian, ω, t)
        add_hamiltonian!(builder(n).hamiltonian, Hn)
    end
    frequency = harmonic_frequency_coeff(ω, t)
    for (ch, γ) in L.channels
        jumps = harmonics(ch.jump, ω, t)
        partners = is_diagonal(ch) ? jumps : harmonics(ch.partner, ω, t)
        for (r, γr) in coefficient_harmonics(γ, Num(ω), frequency, t),
                (m, Jm) in jumps, (l, Kl) in partners
            add_channel!(builder(r + m - l), Jm, Kl, γr, copy(ch.scope), ch.ne)
        end
    end
    out = Dict{Int, Lindbladian{QAdd}}()
    for (n, b) in builders
        Ln = build(b)
        iszero(Ln) || (out[n] = Ln)
    end
    return out
end

harmonic(L::Lindbladian{QAdd}, ω::Real, t::Num, n::Integer) =
    get(() -> zero(Lindbladian{QAdd}), harmonics(L, ω, t), Int(n))

"""
    secular(x, ω, t)

The time-independent (rotating-wave) part of `x`, `harmonic(x, ω, t, 0)`. Every
oscillating component must be an exact integer harmonic of `ω`; otherwise an
`ArgumentError` is thrown. For a `Lindbladian` the result is again a `Lindbladian`. Its
rate matrices are time averages of those of `L`, so they stay positive semidefinite when
those of `L` are at all times. See [`harmonics`](@ref).
"""
secular(x::Union{Number, Coeff, QSym, QAdd, SuperOperator, Lindbladian{QAdd}}, ω::Real, t::Num) =
    harmonic(x, ω, t, 0)
