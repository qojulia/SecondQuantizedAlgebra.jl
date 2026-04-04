using SecondQuantizedAlgebra
using Test

@testset "cluster" begin
    # Setup cluster system
    order = 4
    N_c = 3  # number of clusters
    N = [Parameter(Symbol(:N_, i)) for i in 1:N_c]

    hf = FockSpace(:cavity)
    ha = [NLevelSpace(Symbol(:atoms, j), 3) for j in 1:N_c]
    ha_c = [ClusterSpace(ha[j], N[j], order) for j in 1:N_c]
    h = ⊗(hf, ha_c...)

    # Define the fundamental operators
    a = Destroy(h, :a, 1)
    S(i, j, c) = Transition(h, Symbol(:σ, c), i, j, 1+c)  # c=cluster

    @testset "Cluster Detection" begin
        @test SecondQuantizedAlgebra.has_cluster(h)
        @test !(SecondQuantizedAlgebra.has_cluster(hf))
    end

    @testset "Cluster Operator Properties" begin
        @test S(2, 2, 1) ≠ S(2, 2, 2) ≠ S(2, 2, 3)
        @test length(S(2, 2, 1)) == length(S(2, 2, 2)) == length(S(2, 2, 3)) == order
    end

    @testset "Cluster Operator Commutation" begin
        @test isequal(S(2, 2, 2)[1]*S(2, 2, 1)[2], S(2, 2, 1)[2]*S(2, 2, 2)[1])
    end

    @testset "Acts On Properties" begin
        @test acts_on(S(2, 2, 1)[1]) == SecondQuantizedAlgebra.ClusterAon(2, 1)
        @test acts_on(S(2, 2, 2)[2]) == SecondQuantizedAlgebra.ClusterAon(3, 2)
        @test acts_on(S(2, 2, 1)[1]) < acts_on(S(2, 2, 1)[2]) < acts_on(S(2, 2, 2)[1])
    end
end
