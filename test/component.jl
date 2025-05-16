@testset "component" begin
    @testset "componentfunction" begin

        design = UnfoldSim.SingleSubjectDesign(; conditions = Dict(:duration => 10:-1:5))

        mybasisfun = design -> (collect.(range.(1, generate_events(design).duration)))
        signal = LinearModelComponent(;
            basis = (mybasisfun, 15),
            formula = @formula(0 ~ 1),
            β = [1],
        )

        erp = UnfoldSim.simulate_component(StableRNG(1), signal, design)

        @test size(erp) == (15, 6)
        @test all(erp[11:15, :] .== 0)
        @test erp[1:9, 2] == collect(1.0:9)

        # test shorter cut
        signal = LinearModelComponent(;
            basis = (mybasisfun, 5),
            formula = @formula(0 ~ 1),
            β = [1],
        )

        erp = UnfoldSim.simulate_component(StableRNG(1), signal, design)
        @test size(erp) == (5, 6)
        @test !any(erp .== 0)



    end
    @testset "LMM" begin
        @test UnfoldSim.weight_σs(Dict(:subj => [1, 2]), 0.5, 1.0).subj ==
              LowerTriangular([0.5 0; 0 1.0])
        @test UnfoldSim.weight_σs(Dict(:subj => [1, 2]), 0.5, 2.0).subj ==
              LowerTriangular([0.25 0; 0 0.5])
        @test UnfoldSim.weight_σs(Dict(:item => [1], :subj => [1, 2]), 0.5, 1.0).subj ==
              LowerTriangular([0.5 0; 0 1.0])
        @test UnfoldSim.weight_σs(Dict(:subj => [1, 2, [1 0; 0 1]]), 0.5, 2.0).subj ==
              LowerTriangular([0.25 0; 0 0.5])
        @test UnfoldSim.weight_σs(Dict(:subj => [1, 2, [1 0.5; 0.5 1]]), 1.0, 1.0).subj ==
              create_re(1, 2; corrmat = [1 0.5; 0.5 1])
        @test UnfoldSim.weight_σs(Dict(:subj => [1, 2, [1 0.5; 0.5 1]]), 1.0, 2.0).subj ==
              create_re(1, 2; corrmat = [1 0.5; 0.5 1]) ./ 2
    end

    @testset "DriftComponent" begin
      # Test UnfoldSim.simulate_component(rng, c::DriftComponent, design::AbstractDesign)
      boundary = 1.0
      model_parameter = Dict(:boundary => boundary);
      c = DriftComponent(500, 500, KellyModel, model_parameter);
      design_single = UnfoldSim.SingleSubjectDesign(conditions = Dict(:condition => [1]));
      design_seq = UnfoldSim.SequenceDesign(design_single,"SCR_");
      result_traces = UnfoldSim.simulate_component(StableRNG(1),c,design_seq)
      
      @test size(result_traces) == (501, 3)
      @test any(result_traces .== 0)
      @test any(result_traces .>= boundary)

      # Test UnfoldSim.simulate_component(rng, c::DriftComponent, design::AbstractDesign)
      boundary = 1.0
      model_parameter = Dict(:boundary => boundary);
      c = DriftComponent(500, 500, KellyModel, model_parameter);
      design_single = UnfoldSim.SingleSubjectDesign(conditions = Dict(:drift_rate => [0.5, 0.8], :condition => [1]));
      design_seq = UnfoldSim.SequenceDesign(design_single,"SCR_");
      result_traces = UnfoldSim.simulate_component(StableRNG(1),c,design_seq)

      @test size(result_traces) == (501, 6)
      @test any(result_traces .== 0)
      @test any(result_traces .>= boundary)

      # Test calculate_response_times_for_ssm(rng, component::DriftComponent, design::AbstractDesign)
      model_parameter = Dict();
      c = DriftComponent(500, 500, KellyModel, model_parameter);
      design_single = UnfoldSim.SingleSubjectDesign(conditions = Dict(:drift_rate => [0.5, 0.8], :condition => [1]));
      design_seq = UnfoldSim.SequenceDesign(design_single,"SCR_");
      sub_design = UnfoldSim.SubselectDesign(design_seq, 'C')
      result_rts = UnfoldSim.calculate_response_times_for_ssm(StableRNG(1),c,sub_design)
      @test size(result_rts) == (2,)
      @test isapprox(result_rts, [399.6903067274333, 388.89617910657597], atol=1e-8)

      # Test get_model_parameter(rng, evt, d::Dict)
      rng = StableRNG(1)
      model_parameter = Dict(:drift_rate => "drift_rate");
      c = DriftComponent(500, 500, KellyModel, model_parameter);
      drift_rates = [0.5, 0.8]
      design_single = UnfoldSim.SingleSubjectDesign(conditions = Dict(:drift_rate => drift_rates, :condition => [1]));
      events = UnfoldSim.generate_events(rng, design_single)
      for (i, evt) in enumerate(eachrow(events))
            parameters = UnfoldSim.get_model_parameter(rng, evt, c.model_parameters)
            @test parameters[:drift_rate] == drift_rates[i]
      end
        
    @testset "get_basis" begin


      rng = StableRNG(1)
        design = UnfoldSim.SingleSubjectDesign(; conditions = Dict(:duration => 10:-1:5))
        mybasisfun =
            (rng, design) -> (collect.(range.(1, generate_events(rng, design).duration)))
        signal = LinearModelComponent(;
            basis = (mybasisfun, 15),
            formula = @formula(0 ~ 1),
            β = [1],
        )
        @test UnfoldSim.get_basis(deepcopy(rng), signal, design) ==
              UnfoldSim.get_basis(signal, design)

        shuffle_design = UnfoldSim.SingleSubjectDesign(;
            conditions = Dict(:duration => 10:-1:5),
            event_order_function = shuffle,
        )
        # with same seed => equal result
        @test UnfoldSim.get_basis(StableRNG(1), signal, shuffle_design) ==
              UnfoldSim.get_basis(StableRNG(1), signal, shuffle_design)
        # with different seed => unequal result
        @test UnfoldSim.get_basis(StableRNG(1), signal, shuffle_design) !=
              UnfoldSim.get_basis(StableRNG(2), signal, shuffle_design)

    end

    @testset "max/min offset" begin
        # test max/min offset
        smin10 = LinearModelComponent(;
            basis = [1, 2, 3],
            formula = @formula(0 ~ 1),
            β = [1],
            offset = -10,
        )
        splus5 = LinearModelComponent(;
            basis = [1, 2, 3],
            formula = @formula(0 ~ 1),
            β = [1],
            offset = 5,
        )
        @test UnfoldSim.get_offset(smin10) == -10
        @test UnfoldSim.maxoffset([smin10, splus5]) == 5
        @test UnfoldSim.minoffset([smin10, splus5]) == -10
        @test UnfoldSim.minoffset(Dict('A' => [smin10, splus5])) == -10
        @test UnfoldSim.maxoffset(Dict('A' => [smin10, smin10], 'B' => [splus5, splus5])) ==
              5
        # test that you can have a super large negative offset and dont run into errors (e.g. an event cannot even run in the issue to start before simulation time = 0)

        smin10000 = LinearModelComponent(;
            basis = [1, 2, 3],
            formula = @formula(0 ~ 1),
            β = [1],
            offset = -10_000,
        )
        design = UnfoldSim.SingleSubjectDesign(; conditions = Dict(:duration => 10:-1:5))
        d, e = simulate(design, smin10000, UniformOnset(50, 0))
        @test length(d) > 10_000
        @test e.latency[1] > 10_000

>>>>>>> fix#124
    end
end
