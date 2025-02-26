using UnfoldSim
@testset "hanning" begin
    @test UnfoldSim.hanning(0.021, 0.04, 1000)[40] == 1.0
    @test UnfoldSim.hanning(0.011, 0.04, 1000)[40] == 1.0
    @test UnfoldSim.hanning(0.011, 0.02, 1000)[20] == 1.0
    @test_throws Exception UnfoldSim.hanning(0.011, 0.0, 1000)
end