using StableRNGs
@testset "noise" begin

    for n in [PinkNoise RedNoise WhiteNoise]

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

@testset "simulate_noise - types & scaling" begin
    rng = StableRNG(1)
    n = 64

    # NoNoise -> zeros
    z = simulate_noise(deepcopy(rng), NoNoise(), n)
    @test size(z) == (n,)
    @test all(z .== 0.0)


    # WhiteNoise scaling: use identical RNG state (via deepcopy) to compare scaling
    a = simulate_noise(deepcopy(rng), WhiteNoise(noiselevel = 1), n)
    b = simulate_noise(deepcopy(rng), WhiteNoise(noiselevel = 2), n)
    @test all(isapprox.(2 .* a, b; atol = 0, rtol = 1e-12))

    # PinkNoise & RedNoise scaling (relies on deterministic t.func with same RNG copy)
    p1 = simulate_noise(deepcopy(rng), PinkNoise(noiselevel = 1), n)
    p2 = simulate_noise(deepcopy(rng), PinkNoise(noiselevel = 3), n)
    @test all(isapprox.(3 .* p1, p2; atol = 0, rtol = 1e-12))

    r1 = simulate_noise(deepcopy(rng), RedNoise(noiselevel = 1), n)
    r2 = simulate_noise(deepcopy(rng), RedNoise(noiselevel = 0.5), n)
    @test all(isapprox.(0.5 .* r1, r2; atol = 0, rtol = 1e-12))

    # ExponentialNoise: returns correct length and differs for different seeds
    e1 = simulate_noise(deepcopy(rng), ExponentialNoise(noiselevel = 1, ν = 1.2), n)
    e2 =
        simulate_noise(deepcopy(StableRNG(2)), ExponentialNoise(noiselevel = 1, ν = 1.2), n)
    @test size(e1) == (n,)
    @test size(e2) == (n,)
    @test !(all(isapprox.(e1, e2; atol = 0, rtol = 1e-12)))
end

@testset "add_noise! behaviour" begin
    rng = StableRNG(42)
    n = 32
    # assume default Simulation() constructor exists in package
    sim = Simulation(
        SingleSubjectDesign(),
        [LinearModelComponent(basis = [1, 2], formula = @formula(0 ~ 1), β = [1])],
        UniformOnset(),
        NoNoise(),
    )

    # add_noise! with NoNoise: leaves signal unchanged
    s = zeros(Float64, n)
    UnfoldSim.add_noise!(rng, NoNoise(), s, sim)
    @test all(s .== 0.0)

    # add_noise! with WhiteNoise: resulting signal equals original + simulate_noise(deepcopy(rng), ...)
    s2 = fill(1.0, n)
    # compute expected noise using same base rng state (deepcopy inside add_noise! will produce same sequence)
    expected_noise =
        simulate_noise(deepcopy(rng), WhiteNoise(noiselevel = 1), size(s2), sim)
    expected = copy(s2) .+ reshape(expected_noise, size(s2))
    UnfoldSim.add_noise!(rng, WhiteNoise(noiselevel = 1), s2, sim)
    @test all(isapprox.(s2, expected; atol = 0, rtol = 1e-12))

    # test the new tuple option for simulate_noise to be the same as no-tuple
    a = simulate_noise(deepcopy(rng), WhiteNoise(noiselevel = 1), (10, 2), sim)
    b = simulate_noise(deepcopy(rng), WhiteNoise(noiselevel = 1), 20, sim)
    @test a == b
    @test size(a) == (20,)
end