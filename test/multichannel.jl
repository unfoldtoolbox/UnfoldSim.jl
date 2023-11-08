using UnfoldSim
using StableRNGs
using Test
#---
hart = UnfoldSim.headmodel()
mg = magnitude(hart)


design = SingleSubjectDesign(conditions=Dict(:condA=>["levelA","levelB"])) |> x->RepeatDesign(x,5);
signal = LinearModelComponent(;basis=[0,1,2,3,0],formula = @formula(0~1+condA),β = [2,5]);
onset = UniformOnset(;width=20,offset=4);


@testset "multichannel" begin
    mc = UnfoldSim.MultichannelComponent(signal,[-2,-1,0,1,2,3,4])
    s = simulate(StableRNG(1),mc,design)
    @test size(s) == (7,5,10) # 7 chanels, 5 timepoints, 10 trials
    @ŧest all(s[3,:,:] .== 0) # 3 is 0 in the projection
    @test s[4,:,3] == [0,1,2,3,0].*2 # basis * β[1]

    # test different projection size, should error
    mcA = UnfoldSim.MultichannelComponent(signal,[-2,-1,0,1,2,3,4])
    mcB = UnfoldSim.MultichannelComponent(signal,[-2,-1,0,1,2,3,4,5])
    @test_throws AssertionError UnfoldSim.n_channels([mcA,mcB])
    @test UnfoldSim.n_channels([mcA,mcA]) == 7

    simulation = Simulation(design, mc,  onset, NoNoise());
    data,events = simulate(StableRNG(1),simulation)

    @test all(data[:,events.latency[1]+1] .≈ [-4. -2. 0 2. 4. 6. 8.]')

    @test data[1,:] == data[2,:]*2
    @test all(data[3,:] .== 0)
end

@testset "multichannel with helper" begin

    mc = UnfoldSim.MultichannelComponent(signal, hart=>"Left Middle Temporal Gyrus, posterior division")
    data,events = simulate(StableRNG(1),design, mc,  onset, NoNoise())
    @test size(data,1) == 231
end