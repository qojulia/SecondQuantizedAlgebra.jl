using Test
using SecondQuantizedAlgebra
using SecondQuantizedAlgebra: _site_compare, _can_commute, _commute_pair,
    _reduce_pair, _ground_state_expand, SiteCmp, Less, Greater, Equal, Undetermined,
    _CNUM_ONE, _CNUM_ZERO, _EMPTY_NE

@testset "Operator hooks smoke" begin
    h = FockSpace(:f)
    a = Destroy(h, :a)
    ad = adjoint(a)
    # Same-site cross-type pairs return Equal; canonical direction is in _can_commute.
    @test _site_compare(a, ad, _EMPTY_NE) === Equal
    @test _site_compare(ad, a, _EMPTY_NE) === Equal
    @test _can_commute(a, ad) === false
    @test _can_commute(ad, a) === true
    sw = _commute_pair(a, ad)
    @test sw[1] === ad && sw[2] === a

    ha = NLevelSpace(:atom, 2)
    σ12 = Transition(ha, :σ, 1, 2)
    σ21 = Transition(ha, :σ, 2, 1)
    @test _site_compare(σ12, σ21, _EMPTY_NE) === Equal
    @test _can_commute(σ12, σ21) === false
    red = _reduce_pair(σ21, σ12)   # σ²¹·σ¹² = σ²²
    @test red isa Tuple
    @test red[1].i == 2 && red[1].j == 2
    @test _ground_state_expand(Transition(ha, :σ, 1, 1)) !== nothing
end
