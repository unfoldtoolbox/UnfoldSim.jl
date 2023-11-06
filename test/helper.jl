@testset "helper" begin
    

    @testset "padarray" begin
        @test padarray([2,2],(-3,2),-1) == [-1,-1,-1,2,2,-1,-1]
        @test padarray([2,2],-2,-1) == [-1,-1,2,2]
        @test padarray([2,2],2,-1) == [2,2,-1,-1]
        @test padarray([2,2],2,-1) == [2,2,-1,-1]

    end

    @testset "closest_srcs" begin
        UnfoldSim.closest_srcs([0,0,1],[0 0 0.5; -1 -1 -1; -3 -3 -3])
    
    end
end