@testset "onset" begin
    dummydesign = gen_debug_design(;n_subj=300,n_item=1000)
    @testset "UniformOnset" begin
        unifOnset = UniformOnset(;offset=100,width=50)

        # test random numbers are UniformOnset
        rand_vec = UnfoldSim.rand_onsets(StableRNG(1),unifOnset,dummydesign)
        @test size(rand_vec) == (1000,300)
        @test minimum(rand_vec) ≈ 100
        @test maximum(rand_vec) ≈ 150
    end
    @testset "LogNormalOnset" begin
        # test basics
        logNormalOnset = LogNormalOnset(;μ=4,σ=1)
        rand_vec = UnfoldSim.rand_onsets(StableRNG(1),logNormalOnset, dummydesign)
        @test size(rand_vec) == (1000,300)
        @test isapprox(mean(log.(rand_vec)), 4;atol=0.01)
        @test minimum(rand_vec) > 0
        @test isapprox(std(log.(rand_vec)),1;atol=0.01)

        # test offset
        logNormalOnset = LogNormalOnset(;μ=4,σ=1,offset=100)
        rand_vec = UnfoldSim.rand_onsets(StableRNG(1),logNormalOnset, dummydesign)
        @test minimum(rand_vec) > 100
        # test Truncated
        logNormalOnset =
            LogNormalOnset(; μ = 4, σ = 1, truncate_lower = 10, truncate_upper = 100)
        rand_vec = UnfoldSim.rand_onsets(StableRNG(1), logNormalOnset, dummydesign)
        @test maximum(rand_vec) <= 100
        @test minimum(rand_vec) >= 10
    end
    @testset "gen_onsets" begin
        # test accumulate always increasing
        unifOnset = UniformOnset(;offset=0,width=50)

        accumOnset = UnfoldSim.generate(StableRNG(1),unifOnset,gen_debug_simulation(onset=unifOnset))
        @test  all(diff(accumOnset,dims=1) .>= 0)
    end

end
