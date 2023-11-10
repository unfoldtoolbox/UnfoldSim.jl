hart = UnfoldSim.headmodel()
@testset "hartmut" begin
    @test length(hart.electrodes["label"]) == 231-4
    @test length(hart.cortical["label"]) == 2004
    @test length(hart.artefactual["label"]) == 4260
end
@testset "leadfield/magnitude" begin
    lf = leadfield(hart)
    @test isa(lf,AbstractArray)
    mg = magnitude(hart,type="norm")
    mg_lf = magnitude(lf)
    @test mg == mg_lf

    or = orientation(hart)
    
    @test isa(lf,AbstractArray)
    mg_man = magnitude(lf,or)
    mg = magnitude(hart)
    @test mg_man == mg
    @test size(mg) == (231-4,2004)
    
    
end