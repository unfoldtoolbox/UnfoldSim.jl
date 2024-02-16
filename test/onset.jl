@testset "onset" begin
    dummydesign = gen_debug_design(; n_subjects = 300, n_item = 1000)
    @testset "UniformOnset" begin
        uniform_onset = UniformOnset(; offset = 100, width = 50)

        # test random numbers are UniformOnset
        rand_vec = UnfoldSim.simulate_interonset_distances(
            StableRNG(1),
            uniform_onset,
            dummydesign,
        )
        @test size(rand_vec) == (1000, 300)
        @test minimum(rand_vec) ≈ 100
        @test maximum(rand_vec) ≈ 150
    end
    @testset "LogNormalOnset" begin
        # test basics
        lognormal_onset = LogNormalOnset(; μ = 4, σ = 1)
        rand_vec = UnfoldSim.simulate_interonset_distances(
            StableRNG(1),
            lognormal_onset,
            dummydesign,
        )
        @test size(rand_vec) == (1000, 300)
        @test isapprox(mean(log.(rand_vec)), 4; atol = 0.01)
        @test minimum(rand_vec) > 0
        @test isapprox(std(log.(rand_vec)), 1; atol = 0.01)

        # test offset
        lognormal_onset = LogNormalOnset(; μ = 4, σ = 1, offset = 100)
        rand_vec = UnfoldSim.simulate_interonset_distances(
            StableRNG(1),
            lognormal_onset,
            dummydesign,
        )
        @test minimum(rand_vec) > 100

        # test Truncated
        lognormal_onset =
            LogNormalOnset(; μ = 4, σ = 1, truncate_lower = 10, truncate_upper = 100)
        rand_vec = UnfoldSim.simulate_interonset_distances(
            StableRNG(1),
            lognormal_onset,
            dummydesign,
        )
        @test maximum(rand_vec) <= 100
        @test minimum(rand_vec) >= 10
    end
    @testset "sim_onsets" begin
        # test accumulate always increasing
        uniform_onset = UniformOnset(; offset = 0, width = 50)

        accumulated_onset = UnfoldSim.simulate_onsets(
            StableRNG(1),
            uniform_onset,
            gen_debug_simulation(onset = uniform_onset),
        )
        @test all(diff(accumulated_onset, dims = 1) .>= 0)
    end

end
