@testset "component" begin

    @testset "LMM" begin
        @test UnfoldSim.weight_σs(Dict(:subj=>[1,2]),0.5,1.).subj == LowerTriangular([0.5  0; 0 1.0])
        @test UnfoldSim.weight_σs(Dict(:subj=>[1,2]),0.5,2.).subj == LowerTriangular([0.25 0; 0 0.5])
        @test UnfoldSim.weight_σs(Dict(:item=>[1],:subj=>[1,2]),0.5,1.).subj == LowerTriangular([0.5  0; 0 1.0])
        @test UnfoldSim.weight_σs(Dict(:subj=>[1,2,[1 0; 0 1]]),0.5,2.).subj == LowerTriangular([0.25 0; 0 0.5])
        @test UnfoldSim.weight_σs(Dict(:subj=>[1,2,[1 0.5; 0.5 1]]),1.,1.).subj == create_re(1,2;corrmat=[1 0.5; 0.5 1])
        @test UnfoldSim.weight_σs(Dict(:subj=>[1,2,[1 0.5; 0.5 1]]),1.,2.).subj == create_re(1,2;corrmat=[1 0.5; 0.5 1])./2
    end

end
