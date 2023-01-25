@testset "design" begin

    @testset "SingleSubjectDesign" begin
        des = SingleSubjectDesign(;conditions= Dict(:A=>nlevels(5),:B=>nlevels(2)))
        @test generate(des) == generate(des)
        @test size(generate(des)) == (5*2,2)
        @test names(generate(des)) == ["A","B"]
        @test sum(generate(des).A .== "S2") == 2
        @test sum(generate(des).B .== "S2") == 5
        @test sum((generate(des).B .== "S2").&&(generate(des).A .== "S2")) == 100
        des = SingleSubjectDesign(;conditions= Dict(:A=>nlevels(5),:B=>nlevels(2)),
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

    @testset "RepeatDesign" begin
            
        designOnce = MultiSubjectDesign(;
        n_items=8,
        n_subjects = 12,
        subjects_between =Dict(:cond=>["levelA","levelB"]),
        items_between =Dict(:cond=>["levelA","levelB"]),
        );

        design = RepeatDesign(designOnce,3); 
        @test size(generate(design)) == (8*12*3,4)
        @test size(design) == (8*3,12)
        #--- single sub

        designOnce = SingleSubjectDesign(;conditions= Dict(:A=>nlevels(5),:B=>nlevels(2)))

        design = RepeatDesign(designOnce,3); 
        @test size(generate(design)) == (5*2*3,2)
        @test size(design) == (3,)
    end
    
end

