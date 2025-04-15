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
        @test d[e.latency[1]-10_000] == 1

        smax10000 = LinearModelComponent(;
            basis = [1, 2, 3],
            formula = @formula(0 ~ 1),
            β = [1],
            offset = +10_000,
        )
        d, e = simulate(design, smax10000, UniformOnset(50, 0))
        @test length(d) > 10_000
        @test e.latency[1] < 100
        @test d[e.latency[1]+10_000] == 1


        # if we go back -10_000 and front +10_000, we should get a signal measuring 20_000
        d, e = simulate(design, [smax10000, smin10000], UniformOnset(50, 0))
        @test length(d) > 20_000
        @test length(d) < 25_000 # earlier tests had the signal at 30_000, a bti too long
        @test d[e.latency[1]+10_000] == 1
        @test d[e.latency[1]-10_000] == 1



        smax10 = LinearModelComponent(;
            basis = [1, 2, 3],
            formula = @formula(0 ~ 1),
            β = [1],
            offset = +1000,
        )
        smax20 = LinearModelComponent(;
            basis = [1, 2, 3],
            formula = @formula(0 ~ 1),
            β = [1],
            offset = +2000,
        )

        d, e = simulate(design, [smax10, smax20], UniformOnset(50, 0))
        @test d[e.latency[1]+1000] == 1
        @test d[e.latency[1]+2000] == 1

        smin10 = LinearModelComponent(;
            basis = [1, 2, 3],
            formula = @formula(0 ~ 1),
            β = [1],
            offset = -1000,
        )
        smin20 = LinearModelComponent(;
            basis = [1, 2, 3],
            formula = @formula(0 ~ 1),
            β = [1],
            offset = -2000,
        )

        d, e = simulate(design, [smin10, smin20], UniformOnset(50, 0))
        @test d[e.latency[1]-1000] == 1
        @test d[e.latency[1]-2000] == 1
    end
end
