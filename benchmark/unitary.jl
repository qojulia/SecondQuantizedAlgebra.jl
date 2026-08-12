function benchmark_unitary!(SUITE)
    @variables θ::Real ω::Real t::Real α::Real

    h = FockSpace(:a) ⊗ FockSpace(:b)
    a = Destroy(h, :a, 1)
    b = Destroy(h, :b, 2)
    static = Rotation(a, θ)
    timed = Rotation(a, ω * t, t)
    quadratic = a' * a + a' * b + b' * a
    hamiltonian = ω * a' * a + α * (a + a')

    group = SUITE["Unitary transforms"]
    group["construct"]["rotation"] = @benchmarkable Rotation($a, $θ)
    group["construct"]["displacement"] = @benchmarkable Displace($a, $α)
    group["construct"]["squeeze"] = @benchmarkable Squeeze($a, $θ)
    group["conjugate"]["leaf"] = @benchmarkable conjugate($a, $static)
    group["conjugate"]["quadratic"] = @benchmarkable conjugate($quadratic, $static)
    group["conjugate"]["Hamiltonian"] = @benchmarkable conjugate($hamiltonian, $static)
    group["transform"]["static"] = @benchmarkable transform($hamiltonian, $static)
    group["transform"]["timed"] = @benchmarkable transform($hamiltonian, $timed)
    group["inverse"] = @benchmarkable inv($static)
    group["compose"]["static-static"] =
        @benchmarkable $static * Rotation($a, $θ / 2)
    group["compose"]["static-timed"] = @benchmarkable $static * $timed
    group["compose"]["timed-static"] = @benchmarkable $timed * $static
    group["compose"]["timed-timed"] = @benchmarkable $timed * $timed

    hn = NLevelSpace(:atom, 2)
    σ = Transition(hn, :σ, 1, 2)
    W = [0 1; 1 0]
    Wt = [cos(ω * t) -sin(ω * t); sin(ω * t) cos(ω * t)]
    group["N-level"]["static"] = @benchmarkable Rotation($σ, $W)
    group["N-level"]["timed"] = @benchmarkable Rotation($σ, $Wt, $t)

    phase = SecondQuantizedAlgebra.expim(ω * t)
    phase_expr = phase * a + conj(phase) * a
    group["phase"]["construct"] =
        @benchmarkable SecondQuantizedAlgebra.expim($ω * $t)
    group["phase"]["cancel"] = @benchmarkable $phase * conj($phase)
    group["phase"]["substitute"] = @benchmarkable substitute($phase * $a, Dict($t => 1.0))
    group["phase"]["numeric"] = @benchmarkable substitute($phase * $a, Dict($ω => 2.0, $t => 1.0))
    group["phase"]["to exponential"] =
        @benchmarkable SecondQuantizedAlgebra.exponential_form(cos($ω * $t) * $a)
    group["phase"]["to trigonometric"] =
        @benchmarkable SecondQuantizedAlgebra.trigonometric_form($phase_expr)
    group["phase display"]["terminal"] = @benchmarkable sprint(show, $phase_expr)
    group["phase display"]["LaTeX"] = @benchmarkable latexify($phase_expr)

    plain = (α^2 + ω * t + 1) * a
    trig = (cos(ω * t)^2 + sin(ω * t)^2) * a
    hyperbolic = (cosh(θ)^2 - sinh(θ)^2) * a
    high_power = (cos(θ)^20 + sin(θ)^20) * a
    group["coefficient reduction"]["no trig"] = @benchmarkable simplify($plain)
    group["coefficient reduction"]["quadratic"] = @benchmarkable simplify($trig)
    group["coefficient reduction"]["hyperbolic"] = @benchmarkable simplify($hyperbolic)
    group["coefficient reduction"]["high power"] = @benchmarkable simplify($high_power)
    return SUITE
end
