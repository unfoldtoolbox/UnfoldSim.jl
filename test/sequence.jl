@testset "Check Sequences" begin
    @test isa(UnfoldSim.check_sequence("bla_"), String)
    @test isa(UnfoldSim.check_sequence("bla"), String)
    @test_throws AssertionError UnfoldSim.check_sequence("b_la_")
    @test_throws AssertionError UnfoldSim.check_sequence("b_la")
    @test_throws AssertionError UnfoldSim.check_sequence("_bla")

    @test length(UnfoldSim.sequencestring(StableRNG(1), "A{10,10}")) == 10
    @test length(UnfoldSim.sequencestring(StableRNG(1), "A{10,10}B")) == 11
    @test length(UnfoldSim.sequencestring(StableRNG(1), "A{10,20}")) >= 10
end

@testset "Simulate Sequences" begin



    design = SingleSubjectDesign(conditions = Dict(:condition => ["one", "two"]))
    design = SequenceDesign(design, "SCR_", StableRNG(1))
    evt = generate_events(design)
    @test size(evt, 1) == 6
    @test evt.event == ['S', 'C', 'R', 'S', 'C', 'R']

    design = RepeatDesign(design, 2)
    evt = generate_events(design)
    @test size(evt, 1) == 12
    @test evt.event == ['S', 'C', 'R', 'S', 'C', 'R', 'S', 'C', 'R', 'S', 'C', 'R']


    # repeat first, then sequence => same sequence
    design = SingleSubjectDesign(conditions = Dict(:condition => ["A", "B"]))
    design = RepeatDesign(design, 2)
    design = SequenceDesign(design, "S[ABCD]", StableRNG(2))
    evt = generate_events(design)

    @test all(evt.event[2:2:end] .== 'B')


    # sequence first, then repeat => different sequence for each repetition
    design = SingleSubjectDesign(conditions = Dict(:condition => ["A", "B"]))
    design = SequenceDesign(design, "S[ABCD]", StableRNG(2))
    design = RepeatDesign(design, 2)
    evt = generate_events(design)
    @test !all(evt.event[2:2:end] .== 'B')
end
