using UnfoldSim
@testset "hanning" begin
    @test UnfoldSim.hanning(0.021, 0.04, 1000)[40] == 1.0
    @test UnfoldSim.hanning(0.011, 0.04, 1000)[40] == 1.0
    @test UnfoldSim.hanning(0.011, 0.02, 1000)[20] == 1.0
    @test isapprox(argmax(UnfoldSim.hanning(0.021, 0.04, 256)) / 256, 0.0390625)
    @test_throws Exception UnfoldSim.hanning(0.011, 0.0, 1000)
end

@testset "p100,N170,p300,n400" begin
    sfreq = 1000
    @test argmax(p100(; sfreq)) == 0.1 * sfreq
    @test argmin(n170(; sfreq)) == 0.169 * sfreq # Why not 0.17? Because the peak of the function is in between samples (169 and 170)
    @test argmax(p300(; sfreq)) == 0.3 * sfreq
    @test argmin(n400(; sfreq)) == 0.4 * sfreq
end
