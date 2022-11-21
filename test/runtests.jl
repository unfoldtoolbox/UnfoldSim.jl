using UnfoldSim
include("setup.jl")

@testset "UnfoldSim.jl" begin
    include("component.jl")
    include("design.jl")
    include("noise.jl")
    include("onset.jl")
    include("simulation.jl")
end
