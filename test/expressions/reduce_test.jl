using Test
using SecondQuantizedAlgebra
using Symbolics: Symbolics, @variables, Num
import SymbolicUtils
import SecondQuantizedAlgebra: ParamRelation,
    to_cnum, to_num, is_poly, iszero_cnum, reduce_trig, reduce_cnum, trig_relations,
    reduce_via_transient, sym_trig_relations!, from_raw, simplify

# `cos^2+sin^2` and `cosh^2-sinh^2` are identities, so they are folded on the parameter
# polynomial tier instead of by the CAS: two orders of magnitude cheaper and, unlike
# `Symbolics.simplify`, still complete above degree 2. These tests pin what the reduction
# folds, what the cancellation gate deliberately leaves alone, and the tier contract that
# keeps the CAS reachable for everything else.

# Every coefficient here is real, and `isequal(::Complex{Num}, ::Num)` is ambiguous, so
# compare real parts.
realpart(c) = real(to_num(c))
reduced(ex) = realpart(reduce_trig(to_cnum(ex)))
relations(ex) = trig_relations(to_cnum(ex).tail.terms)

@testset "Parameter relations" begin
    @variables θ φ r ω t g κ

    @testset "trigonometric folding" begin
        @test isequal(reduced(cos(θ)^2 + sin(θ)^2), Num(1))
        @test isequal(reduced(cos(θ)^4 + 2cos(θ)^2 * sin(θ)^2 + sin(θ)^4), Num(1))
        @test isequal(reduced(cosh(r)^2 - sinh(r)^2), Num(1))
        @test isequal(reduced(cosh(r)^4 - 2cosh(r)^2 * sinh(r)^2 + sinh(r)^4), Num(1))
        # A coefficient carrying both a trigonometric and a hyperbolic pair.
        @test isequal(reduced(cos(θ)^2 * cosh(r)^2 + sin(θ)^2 * cosh(r)^2 - sinh(r)^2), Num(1))
    end

    @testset "two independent pairs" begin
        # What `Rotation(a, θ) * Rotation(b, φ)` leaves on `a'*a*b'*b`: one relation per
        # pair, so a single relation is not enough.
        ex = (cos(θ)^2 + sin(θ)^2) * (cos(φ)^2 + sin(φ)^2)
        @test length(relations(ex)) == 2
        @test isequal(reduced(ex), Num(1))
    end

    @testset "cancellation gate" begin
        # The gate keeps a rewrite only when it shortens the term list, so an expansion
        # that does not collapse is left in its original form.
        @test isequal(reduced(cos(θ)^3 * sin(θ)), cos(θ)^3 * sin(θ))
        @test isequal(reduced(cos(θ)^4 + sin(θ)^4), cos(θ)^4 + sin(θ)^4)
        # Ungated, the same input does expand: the two modes genuinely differ.
        c = to_cnum(cos(θ)^4 + sin(θ)^4)
        rels = trig_relations(c.tail.terms)
        @test !isequal(realpart(reduce_cnum(c, rels, false)), realpart(c))
        @test isequal(realpart(reduce_cnum(c, rels, true)), realpart(c))
    end

    @testset "no relation, no rewrite" begin
        # A radical of `cos` carries no Pythagorean identity.
        @test isequal(reduced(sqrt(cos(θ))), sqrt(cos(θ)))
        # A `cos` whose partner is absent, and a partner over a different argument.
        @test isempty(relations(cos(θ)^2 + sin(φ)^2))
        @test isequal(reduced(cos(θ)^2 + sin(φ)^2), cos(θ)^2 + sin(φ)^2)
        # Coefficients with no trigonometric factor at all are returned untouched.
        @test isequal(reduced(g * κ^2), g * κ^2)
    end

    @testset "conj-wrapped pair" begin
        # `conj_atom` wraps `Number`-symtype atoms, so a conjugated coefficient can hold
        # only `(conj(cos u), conj(sin u))`; `conj(cos u)^2 + conj(sin u)^2 = conj(1) = 1`.
        @variables u::Number
        uu = SymbolicUtils.unwrap(u)
        @test length(relations(conj(cos(uu))^2 + conj(sin(uu))^2)) == 1
        @test isequal(reduced(conj(cos(uu))^2 + conj(sin(uu))^2), Num(1))
        # A mismatched wrapping is not a pair.
        @test isempty(relations(conj(cos(uu))^2 + sin(uu)^2))
    end

    @testset "composite argument via the transient swap" begin
        # `cos(ω*t)` is not an atom, so its coefficient never reaches the polynomial tier
        # and the CAS folds neither degree 2 nor degree 4 there.
        rel = ParamRelation(cos(ω * t), sin(ω * t), -1)
        deg2 = to_cnum(cos(ω * t)^2 + sin(ω * t)^2)
        deg4 = to_cnum(cos(ω * t)^4 + 2cos(ω * t)^2 * sin(ω * t)^2 + sin(ω * t)^4)
        @test !is_poly(deg2)
        @test !is_poly(deg4)
        @test isequal(realpart(reduce_cnum(deg2, [rel], true)), Num(1))
        @test isequal(realpart(reduce_cnum(deg4, [rel], true)), Num(1))
        # The CAS leaves the degree-4 form alone, which is why this route exists.
        @test !isequal(Symbolics.simplify(realpart(deg4)), Num(1))
        # A coefficient with no reducible power comes back unchanged.
        prod = to_cnum(cos(ω * t) * sin(ω * t))
        @test isequal(realpart(reduce_cnum(prod, [rel], true)), realpart(prod))

        # The automatic reducer follows the same route, including above degree two.
        @test isequal(realpart(reduce_trig(deg4)), Num(1))

        # A user parameter with the old implementation's fixed scratch name must remain a
        # user parameter. Fresh stand-ins are checked against symbols already in the input.
        collision = Symbolics.variable(:__sqa_rel_1)
        with_collision = to_cnum(
            collision * (cos(ω * t)^2 + sin(ω * t)^2),
        )
        @test isequal(realpart(reduce_trig(with_collision)), collision)
    end

    @testset "guarded high-power refusal" begin
        original = to_cnum(cos(θ)^140 + sin(θ)^140)
        # The exact binomial coefficients no longer fit in `Int`; refusing the rewrite is
        # safer than overflowing a coefficient or partially rewriting the expression.
        @test isequal(reduce_trig(original), original)
    end

    @testset "zero is dropped, not stored" begin
        # A canonical `Poly` never reports itself zero, so a reduction reaching zero has to
        # rebuild as a native zero or `iszero`/`isequal` on the result would lie.
        c = to_cnum(cos(θ)^4 + 2cos(θ)^2 * sin(θ)^2 + sin(θ)^4 - 1)
        @test !iszero_cnum(c)
        @test iszero_cnum(reduce_trig(c))
        # ... and the enclosing term disappears from the sum.
        h = FockSpace(:f)
        a = Destroy(h, :a)
        q = (cos(θ)^4 + 2cos(θ)^2 * sin(θ)^2 + sin(θ)^4 - 1) * a' * a
        @test length(q) == 1
        @test iszero(simplify(q))
        @test length(simplify(q)) == 0
    end

    @testset "simplify folds identities" begin
        h = FockSpace(:f)
        a = Destroy(h, :a)
        @test isequal(simplify((cos(θ)^2 + sin(θ)^2) * a' * a), a' * a)
        @test isequal(simplify((cosh(r)^2 - sinh(r)^2) * a), 1 * a)
    end

    @testset "tier contract is unchanged" begin
        # Widening the atom test to cover trigonometric calls of composite arguments would
        # have moved these onto the polynomial tier, where the CAS is never reached, and
        # lost every rewrite beyond the Pythagorean identity.
        for x in (exp(g + κ), sin(g * κ), sqrt(g * κ))
            @test !is_poly(to_cnum(x))
        end
        @test isequal(Symbolics.simplify(cos(θ) * sin(θ)), sin(2θ) / 2)
    end

    @testset "ParamRelation construction" begin
        r1 = ParamRelation(cos(θ), sin(θ), -1)
        @test r1.sign == -1
        @test isequal(r1.hi, cos(θ))
        # Wrapping is idempotent: a `Num` member is not double-wrapped.
        @test SymbolicUtils.unwrap(r1.hi) === SymbolicUtils.unwrap(cos(θ))
        @test ParamRelation(cosh(θ), sinh(θ), 1).sign == 1
    end
end

@testset "Raw coefficient reduction seams" begin
    @variables θ ω t
    relation = ParamRelation(cos(θ), sin(θ), -1)
    polynomial = to_cnum(cos(θ)^2 + sin(θ)^2)
    @test isequal(
        realpart(reduce_via_transient(polynomial, [relation], true)),
        Num(1),
    )

    raw = from_raw(SymbolicUtils.unwrap(cos(ω * t) + sin(ω * t)))
    raw_relations = ParamRelation[]
    @test sym_trig_relations!(raw_relations, raw) === raw_relations
    @test length(raw_relations) == 1

    polynomial_relations = ParamRelation[]
    @test sym_trig_relations!(
        polynomial_relations,
        to_cnum(cos(θ) + sin(θ)),
    ) === polynomial_relations
    @test length(polynomial_relations) == 1
end
