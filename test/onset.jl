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
        lognormal_onset = LogNormalOnset(; μ = 4, σ = 1, truncate_upper = 100)
        rand_vec = UnfoldSim.simulate_interonset_distances(
            StableRNG(1),
            lognormal_onset,
            dummydesign,
        )
        @test maximum(rand_vec) <= 100
        @test minimum(rand_vec) >= 0
    end
    @testset "sim_onsets" begin
        uniform_onset = UniformOnset(; offset = 0, width = 50)

        accumulated_onset = UnfoldSim.simulate_onsets(
            StableRNG(1),
            uniform_onset,
            gen_debug_simulation(onset = uniform_onset),
        )
        # test accumulate always increasing
        @test all(diff(accumulated_onset, dims = 1) .>= 0)

        # test that the first onset is at >=1 (not 0)
        @test accumulated_onset[1] >= 1
    end

    @testset "OnsetFormula" begin

        design =
            SingleSubjectDesign(conditions = Dict(:cond => ["A", "B"])) |>
            x -> RepeatDesign(x, 10000)


        o = UniformOnsetFormula(width_formula = @formula(0 ~ 1 + cond), width_β = [50, 20])
        events = generate_events(design)
        onsets = UnfoldSim.simulate_interonset_distances(StableRNG(1), o, design)
        @test minimum(onsets[1:2:end]) == 0
        @test maximum(onsets[1:2:end]) == 50
        @test minimum(onsets[2:2:end]) == 0
        @test maximum(onsets[2:2:end]) == 70

        o = UniformOnsetFormula(
            offset_formula = @formula(0 ~ 1 + cond),
            offset_β = [50, 20],
            width_β = [50],
        )
        events = generate_events(design)
        onsets = UnfoldSim.simulate_interonset_distances(StableRNG(1), o, design)
        @test minimum(onsets[1:2:end]) == 50
        @test maximum(onsets[1:2:end]) == 100
        @test minimum(onsets[2:2:end]) == 70
        @test maximum(onsets[2:2:end]) == 120


        o = LogNormalOnsetFormula(
            μ_formula = @formula(0 ~ 1 + cond),
            μ_β = [1, 1],
            σ_β = [1],
        )
        events = generate_events(design)
        onsets = UnfoldSim.simulate_interonset_distances(StableRNG(1), o, design)
        @test minimum(onsets[1:2:end]) == 0
        @test maximum(onsets[1:2:end]) < 150
        @test minimum(onsets[2:2:end]) == 0
        @test maximum(onsets[2:2:end]) > 300



    end

    @testset "ShiftOnset" begin
        design =
            SingleSubjectDesign(conditions = Dict(:cond => ["A", "B"])) |>
            x -> RepeatDesign(x, 100)

        o = UniformOnset(width = 50, offset = 10)

        without = UnfoldSim.simulate_interonset_distances(StableRNG(1), o, design)
        with = UnfoldSim.simulate_interonset_distances(
            StableRNG(1),
            ShiftOnsetByOne(o),
            design,
        )
        # ShiftOnsetByOne adds a 0 to the front, thereby the first "non-0" "real" simulated inter onset distance is used for the second event
        @test with[1] == 0

        @test without[1:(end-1)] == with[2:end]


    end
end
