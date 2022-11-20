@testset "onset" begin

    @testset "UniformOnset" begin
        unifOnset = UniformOnset(;offset=100,width=50)

        # test random numbers are UniformOnset
        rand_vec = UnfoldSim.rand_onsets(StableRNG(1),unifOnset,1000,300,nothing)
        @test size(rand_vec) == (1000,300)
        @test minimum(rand_vec) ≈ 100
        @test maximum(rand_vec) ≈ 150
    end
    @testset "LogNormalOnset" begin
        # test basics
        logNormalOnset = LogNormalOnset(;μ=4,σ=1)
        rand_vec = UnfoldSim.rand_onsets(StableRNG(1),logNormalOnset,1000,300,nothing)
        @test size(rand_vec) == (1000,300)
        @test isapprox(mean(log.(rand_vec)), 4;atol=0.01)
        @test minimum(rand_vec) > 0
        @test isapprox(std(log.(rand_vec)),1;atol=0.01)

        # test offset
        logNormalOnset = LogNormalOnset(;μ=4,σ=1,offset=100)
        rand_vec = UnfoldSim.rand_onsets(StableRNG(1),logNormalOnset,1000,300,nothing)
        @test minimum(rand_vec) > 100
        
        # test Truncated
        logNormalOnset = LogNormalOnset(;μ=4,σ=1,truncate_upper=100)
        rand_vec = UnfoldSim.rand_onsets(StableRNG(1),logNormalOnset,1000,300,nothing)
        @test maximum(rand_vec) <= 100
        @test minimum(rand_vec) >= 0 
    end
    @testset "gen_onsets" begin
        # test accumulate always increasing
        unifOnset = UniformOnset(;offset=100,width=50)

        accumOnset = gen_onsets(StableRNG(1),gen_debug_simulation(onset=unifOnset))
        @test  all(diff(accumOnset) .<= 0)
    end

end
