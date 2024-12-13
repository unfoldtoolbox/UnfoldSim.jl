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

end
