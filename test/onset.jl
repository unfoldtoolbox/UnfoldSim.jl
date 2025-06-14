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

    @testset "Sequence/Drift Onset" begin
        rng = StableRNG(1)
        fs = 500
        p3 = LinearModelComponent(;
            basis = UnfoldSim.hanning(Int(0.5 * fs)),
            formula = @formula(0 ~ 1 + condition),
            β = [1.0, 0],
        )

        resp = LinearModelComponent(;
            basis = UnfoldSim.hanning(Int(0.5 * fs)),
            formula = @formula(0 ~ 1 + condition),
            β = [0.5, 0],
            offset = -10,
        )
        sequence_onset = SequenceOnset(
            Dict(
                'S' => UniformOnset(width = 0, offset = 80),
                'C' => DriftOnset(),
                'R' => UniformOnset(width = 0, offset = 120),
            ),
        )
        model_parameter = Dict(:drift_rate => "drift_rate")
        drift = UnfoldSim.DriftComponent(500, 500, KellyModel, model_parameter)
        components = Dict('S' => [p3], 'C' => [drift], 'R' => [resp])
        design_single = UnfoldSim.SingleSubjectDesign(
            conditions = Dict(:drift_rate => [0.5, 0.8], :condition => [1]),
        )
        design_seq = UnfoldSim.SequenceDesign(design_single, "SCR_")
        design_rep = UnfoldSim.RepeatDesign(design_seq, 10)
        simulation = UnfoldSim.Simulation(
            design_rep,
            components,
            sequence_onset,
            UnfoldSim.NoNoise(),
        )

        result_onsets = simulate_onsets(rng, sequence_onset, simulation)

        size(result_onsets) == (60,)
        result_onsets[1] == 121
        result_onsets[2] == 201
        result_onsets[3] == 906

        # Test DriftOnset combined with UniformOnset
        sequence_onset = SequenceOnset(
            Dict(
                'S' => UniformOnset(width = 0, offset = 80),
                'C' => (DriftOnset(), UniformOnset(width = 0, offset = 140)),
                'R' => UniformOnset(width = 0, offset = 120),
            ),
        )
        simulation = UnfoldSim.Simulation(
            design_rep,
            components,
            sequence_onset,
            UnfoldSim.NoNoise(),
        )

        result_onsets = simulate_onsets(rng, sequence_onset, simulation)
        size(result_onsets) == (60,)
        result_onsets[1] == 121
        result_onsets[2] == 201
        result_onsets[3] == 1046
    end

end
