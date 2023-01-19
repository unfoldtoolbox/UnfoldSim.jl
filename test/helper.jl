@testset "helper" begin
    

    @testset "padarray" begin
        @test padarray([2,2],(-3,2),-1) == [-1,-1,-1,2,2,-1,-1]
        @test padarray([2,2],-2,-1) == [-1,-1,2,2]
        @test padarray([2,2],2,-1) == [2,2,-1,-1]
        @test padarray([2,2],2,-1) == [2,2,-1,-1]

    end
end