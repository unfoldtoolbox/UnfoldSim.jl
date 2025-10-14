@testset "gazevec_from_angle_3d" begin
    @test UnfoldSim.gazevec_from_angle_3d(0, 0) ≈ [0; 1; 0]
    @test UnfoldSim.gazevec_from_angle_3d(0,90) ≈ [0; 0; 1]
    @test UnfoldSim.gazevec_from_angle_3d(90,0) ≈ [1; 0; 0]
end

@testset "az_simulation" begin
    # test the simulated values of just an eyemovement using the sample data HREF
    # compare a previously simulated and saved set of values with the freshly simulated output
    leadfields_stepbystep_sim = Matrix(CSV.read("src/simulated_2025-07-17.csv",DataFrame));
    data, events = UnfoldSim.az_simulation();
    @test sum(leadfields[:,1:200] .- data[:,1:200]) == 0.0
end