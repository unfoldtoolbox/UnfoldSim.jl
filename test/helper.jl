@testset "helper" begin


    @testset "pad_array" begin
        @test pad_array([2, 2], (-3, 2), -1) == [-1, -1, -1, 2, 2, -1, -1]
        @test pad_array([2, 2], -2, -1) == [-1, -1, 2, 2]
        @test pad_array([2, 2], 2, -1) == [2, 2, -1, -1]
        @test pad_array([2, 2], 2, -1) == [2, 2, -1, -1]

    end

    @testset "closest_src" begin
        @test UnfoldSim.closest_src([0, 0, 1], [0 0 0.5; -1 -1 -1; -3 -3 -3]) == 1
        @test UnfoldSim.closest_src(
            [[0, 0, 1], [-1, -1, 0]],
            [0 0 0.5; -1 -1 -1; -3 -3 -3],
        ) == [1, 2]


    end
end
