@testset "noise" begin

    for n in [PinkNoise RedNoise WhiteNoise ExponentialNoise]

        noisevec_1 = simulate_noise(StableRNG(1), n(; noiselevel = 1), 123)
        @test size(noisevec_1) == (123,)
        noisevec_0 = simulate_noise(StableRNG(1), n(; noiselevel = 0), 123)
        @test all(noisevec_0 .== 0)
        noisevec_2 = simulate_noise(StableRNG(1), n(; noiselevel = 2), 123)
        @test all(2 .* noisevec_1 .≈ noisevec_2)
    end
    noisevec = simulate_noise(StableRNG(1), NoNoise(), 123)
    @test size(noisevec) == (123,)
end
