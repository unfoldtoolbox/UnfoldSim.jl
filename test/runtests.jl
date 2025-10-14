using UnfoldSim
include("setup.jl")

@testset "UnfoldSim.jl" begin
    include("bases.jl")
    include("component.jl")
    include("design.jl")
    include("noise.jl")
    include("onset.jl")
    include("simulation.jl")
    include("helper.jl")
    include("sequence.jl")
end
