@testset "Check Sequences" begin
    @test isa(UnfoldSim.check_sequence("bla_"), String)
    @test isa(UnfoldSim.check_sequence("bla"), String)
    @test_throws AssertionError UnfoldSim.check_sequence("b_la_")
    @test_throws AssertionError UnfoldSim.check_sequence("b_la")
    @test_throws AssertionError UnfoldSim.check_sequence("_bla")

end
@test length(UnfoldSim.sequencestring(StableRNG(1), "A{10,10}")) == 10
@test length(UnfoldSim.sequencestring(StableRNG(1), "A{10,10}B")) == 11
@test length(UnfoldSim.sequencestring(StableRNG(1), "A{10,20}")) >= 10