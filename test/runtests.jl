using UnfoldSim
include("setup.jl")

@testset "UnfoldSim.jl" begin
    include("bases.jl")
    include("component.jl")
    include("design.jl")
    include("headmodel.jl")
    include("helper.jl")
    include("multichannel.jl")
    include("noise.jl")
    include("onset.jl")
    include("sequence.jl")
    include("simulation.jl")
end
