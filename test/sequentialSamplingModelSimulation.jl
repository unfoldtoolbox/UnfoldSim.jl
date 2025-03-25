@testset "sequentialSamplingModelSimulation" begin
    fs = 500
    Δt = 1 / fs # time step
    tEnd = 1.0 # trial Duration
    time_vec = 0:Δt:tEnd # time base
    max_length = tEnd / Δt
    rng = StableRNG(1)
    @testset "KellyModel" begin
        assert_event_onset = 0.663
        assert_drift_rate = "drift_rate"
        km = KellyModel(event_onset = assert_event_onset, drift_rate = assert_drift_rate)

        @test km.event_onset == assert_event_onset
        @test km.drift_rate == assert_drift_rate
    end

    @testset "KellyModel_simulate_cpp" begin
        boundary = 1.0
        result_rt, result_trace = UnfoldSim.KellyModel_simulate_cpp(
            rng,
            KellyModel(boundary = boundary),
            time_vec,
            Δt,
        )
        @test size(result_rt) == ()
        @test size(result_trace) == (501,)
        @test isapprox(result_rt, 399.6903067274333, atol = 1e-8)
        @test any(result_trace .== 0)
        @test any(result_trace .>= boundary)

        result_sim_rt, result_sim_trace = UnfoldSim.SSM_Simulate(rng, KellyModel(), fs, max_length)
        @test result_rt == result_sim_rt
        @test result_trace == result_sim_trace
    end

    @testset "trace_sequential_sampling_model" begin
        boundary = 1.0
        model_parameter = Dict(:boundary => boundary);
        c = UnfoldSim.DriftComponent(
            max_length,
            fs,
            KellyModel,
            model_parameter,
        )
        design_single = UnfoldSim.SingleSubjectDesign(
            conditions = Dict(:drift_rate => [0.5, 0.8], :condition => [1]),
        )
        design_seq = UnfoldSim.SequenceDesign(design_single, "SCR_")

        result_rts, result_traces = UnfoldSim.trace_sequential_sampling_model(rng, c, design_seq)
        @test size(result_rts) == (6,)
        @test size(result_traces) == (501, 6)
        @test any(result_traces .>= 1.0)
    end

    @testset "SSM_Simulate" begin
        result_rt, result_trace = UnfoldSim.SSM_Simulate(deepcopy(rng), KellyModel(), fs, max_length)
        @test size(result_rt) == ()
        @test size(result_trace) == (501,)
        @test isapprox(result_rt, 399.6903067274333, atol = 1e-8)
        @test any(result_trace .== 0)
        @test any(result_trace .>= boundary)

        result_rt, result_trace = UnfoldSim.SSM_Simulate(deepcopy(rng), DDM(), fs, max_length)
        @test size(result_rt) == ()
        @test size(result_trace) == (501,)
        @test isapprox(result_rt, 223.00000000000003, atol = 1e-8)
        @test any(result_trace .== 0)

        result_rt, result_trace = UnfoldSim.SSM_Simulate(deepcopy(rng), LBA(), fs, max_length)
        @test size(result_rt) == ()
        @test size(result_trace) == (501,)
        @test isapprox(result_rt, 397.0, atol = 1e-8)
        @test any(result_trace .== 0)
    end
end
