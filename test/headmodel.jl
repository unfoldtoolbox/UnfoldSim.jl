hart = UnfoldSim.Hartmut()
@testset "hartmut" begin
    @test length(hart.electrodes["label"]) == 231 - 4
    @test length(hart.cortical["label"]) == 2004
    @test length(hart.artefactual["label"]) == 4260
end
@testset "leadfield/magnitude" begin
    lf = leadfield(hart)
    or = orientation(hart)
    
    mg = magnitude(hart)
    mg_man = magnitude(lf, or)
    
    @test isa(lf, AbstractArray)
    @test mg_man == mg
    @test size(mg) == (231 - 4, 2004)

end
