using Test
using SecondQuantizedAlgebra
using SecondQuantizedAlgebra: _site_compare, _can_commute, _reduce_pair,
    SiteCmp, Less, Equal, Undetermined, Greater, _EMPTY_NE

function _is_canonical(t)
    ops = t.ops
    ne = t.ne
    for i in 1:(length(ops) - 1)
        cmp = _site_compare(ops[i], ops[i + 1], ne)
        if cmp === Greater
            return false
        end
        if cmp === Equal
            _can_commute(ops[i], ops[i + 1]) || return false
            _reduce_pair(ops[i], ops[i + 1]) === nothing || return false
        end
    end
    return true
end

@testset "Canonical-form invariant: `*`" begin
    h = FockSpace(:f) ⊗ NLevelSpace(:a, 3)
    a = Destroy(h, :a, 1)
    σ12 = Transition(h, :σ, 1, 2, 2)
    σ21 = Transition(h, :σ, 2, 1, 2)
    q = a * adjoint(a) * σ12 * σ21
    for (t, _) in q
        @test _is_canonical(t)
    end
end

@testset "Canonical-form invariant: normal_order" begin
    h = FockSpace(:f)
    a = Destroy(h, :a)
    q = a * adjoint(a) * a
    for (t, _) in normal_order(q)
        @test _is_canonical(t)
    end
end

@testset "Canonical-form invariant: commutator" begin
    h = FockSpace(:f)
    a = Destroy(h, :a)
    q = commutator(a, adjoint(a))
    for (t, _) in q
        @test _is_canonical(t)
    end
end

@testset "Canonical-form invariant: substitute" begin
    h = FockSpace(:f)
    a = Destroy(h, :a)
    q = adjoint(a) * a
    sub = SecondQuantizedAlgebra.SymbolicUtils.substitute(q, Dict{Any, Any}(a => 0.5))
    for (t, _) in sub
        @test _is_canonical(t)
    end
end

@testset "Canonical-form invariant: expand_completeness" begin
    h = NLevelSpace(:a, 3)
    σ11 = Transition(h, :σ, 1, 1)
    exp = expand_completeness(σ11)
    for (t, _) in exp
        @test _is_canonical(t)
    end
end
