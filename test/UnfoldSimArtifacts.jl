@testset "gazevec_from_angle_3d" begin
    @test UnfoldSim.gazevec_from_angle_3d(0, 0) ≈ [0; 1; 0]
    @test UnfoldSim.gazevec_from_angle_3d(0,90) ≈ [0; 0; 1]
    @test UnfoldSim.gazevec_from_angle_3d(90,0) ≈ [1; 0; 0]
end