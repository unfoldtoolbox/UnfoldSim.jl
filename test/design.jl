@testset "design" begin

    @testset "SingleSubjectDesign" begin
        des = SingleSubjectDesign(;n_trials = 100,conditions= Dict(:A=>nlevels(5),:B=>nlevels(2)))
        @test generate(des) == generate(des)
        @test size(generate(des)) == (100*5*2,2)
        @test names(generate(des)) == ["A","B"]
        @test sum(generate(des).A .== "S2") == 100*2
        @test sum(generate(des).B .== "S2") == 100*5
        @test sum((generate(des).B .== "S2").&&(generate(des).A .== "S2")) == 100
        des = SingleSubjectDesign(;n_trials = 100,conditions= Dict(:A=>nlevels(5),:B=>nlevels(2)),
        tableModifyFun = x->sort(x,order(:B,rev=true)))
        @test generate(des).B[1] == "S2"
    end

    @testset "MultiSubjectDesign" begin
        des = MultiSubjectDesign(;n_subjects=10,
                                n_items = 100,both_within= Dict(:A=>nlevels(5),:B=>nlevels(2)))
        
        
        @test size(generate(des)) == (10*100*5*2,5)
        @test names(generate(des)) == ["subject","item","A","B","dv"]
        @test sum(generate(des).A .== "S2") == 10*100*2
        @test sum(generate(des).B .== "S2") == 10*100*5

        # check that finally everything is sorted by subject
        des = MultiSubjectDesign(;n_subjects=10,
                                n_items = 100,both_within= Dict(:A=>nlevels(5),:B=>nlevels(2)),
                                tableModifyFun = x->sort(x,order(:dv,rev=true)))
        @test generate(des).subject[1] == "S01"
             
        # check tableModifyFun
        des = MultiSubjectDesign(;n_subjects=10,
        n_items = 100,both_within= Dict(:A=>nlevels(5),:B=>nlevels(2)),
        tableModifyFun = x->sort(x,order(:B,rev=true)))
        @test generate(des).B[1] == "S2"

    end
end

