using UnfoldSim
@testset "hanning" begin
    @test UnfoldSim.hanning(0.021, 0.04, 1000)[41] == 1.0 # why 41 not 40? beacuse round(0.5) = 0 and round(1.5) = 2 -- and we are living on the edge!
    @test UnfoldSim.hanning(0.011, 0.04, 1000)[40] == 1.0
    @test isapprox(UnfoldSim.hanning(0.021, 0.04, 256) ,0.0429688)
    @test UnfoldSim.hanning(0.011, 0.02, 1000)[20] == 1.0
    @test_throws Exception UnfoldSim.hanning(0.011, 0.0, 1000)
end

@testset "p100,N170,p300,n400" begin
    sfreq = 1000
    @test argmax(p100(; sfreq)) == 0.1 * sfreq
    @test argmin(n170(; sfreq)) == 0.17 * sfreq
    @test argmax(p300(; sfreq)) == 0.3 * sfreq
    @test argmin(n400(; sfreq)) == 0.4 * sfreq
end