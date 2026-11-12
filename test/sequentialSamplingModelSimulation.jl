@testset "sequentialSamplingModelSimulation" begin
    fs = 500
    Δt = 1 / fs # time step
    tEnd = 1.0 # trial Duration
    time_vec = 0:Δt:tEnd # time base
    max_length = Int(tEnd / Δt)
    rng = StableRNG(1)
    @testset "KellyModel" begin
        assert_sensor_delay = 0.663
        assert_drift_rate = "drift_rate"
        km = KellyModel(sensor_encoding_delay = assert_sensor_delay, drift_rate = assert_drift_rate)

        @test km.sensor_encoding_delay == assert_sensor_delay
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
        @test isapprox(result_rt, 291.53516402920155, atol = 1e-8)
        @test any(result_trace .== 0)
        @test any(result_trace .>= boundary)

        # same seed as the call above, so the results must match
        result_sim_rt, result_sim_trace = UnfoldSim.SSM_Simulate(StableRNG(1), KellyModel(), fs, max_length)
        @test result_rt == result_sim_rt
        @test result_trace[1:max_length] == result_sim_trace
    end

    @testset "simulate_drift_component" begin
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

        result_rts, result_traces = UnfoldSim.simulate_drift_component(rng, c, design_seq)
        @test size(result_rts) == (6,)
        @test size(result_traces) == (500, 6)
        @test any(result_traces .>= 1.0)
    end

    @testset "SSM_Simulate" begin
        boundary = 1.0
        result_rt, result_trace = UnfoldSim.SSM_Simulate(deepcopy(StableRNG(1)), KellyModel(), fs, max_length)
        @test size(result_rt) == ()
        @test size(result_trace) == (500,)
        @test isapprox(result_rt, 291.53516402920155, atol = 1e-8)
        @test any(result_trace .== 0)
        @test any(result_trace .>= boundary)

        result_rt, result_trace = UnfoldSim.SSM_Simulate(deepcopy(StableRNG(1)), DDM(), fs, max_length)
        @test size(result_rt) == ()
        @test size(result_trace) == (500,)
        @test isapprox(result_rt, 223.00000000000003, atol = 1e-8)
        @test any(result_trace .== 0)

        result_rt, result_trace = UnfoldSim.SSM_Simulate(deepcopy(StableRNG(1)), LBA(), fs, max_length)
        @test size(result_rt) == ()
        @test size(result_trace) == (500,)
        @test isapprox(result_rt, 397.0, atol = 1e-8)
        @test any(result_trace .== 0)
    end
end
