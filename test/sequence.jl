using Automa, Random, Test
@testset "Check Sequences" begin
    @test isa(UnfoldSim.check_sequence("bla_"), String)
    @test isa(UnfoldSim.check_sequence("bla"), String)
    @test_throws AssertionError UnfoldSim.check_sequence("b_la_")
    @test_throws AssertionError UnfoldSim.check_sequence("b_la")
    @test_throws AssertionError UnfoldSim.check_sequence("_bla")

    @test length(UnfoldSim.evaluate_sequencestring(StableRNG(1), "A{10,10}")) == 10
    @test length(UnfoldSim.evaluate_sequencestring(StableRNG(1), "A{10,10}B")) == 11
    @test length(UnfoldSim.evaluate_sequencestring(StableRNG(1), "A{10,20}")) >= 10
end

@testset "Simulate Sequences" begin



    design = SingleSubjectDesign(conditions = Dict(:condition => ["one", "two"]))
    design = SequenceDesign(design, "SCR_")
    evt = generate_events(StableRNG(1), design)
    @test size(evt, 1) == 6
    @test evt.event == ['S', 'C', 'R', 'S', 'C', 'R']

    design = RepeatDesign(design, 2)
    evt = generate_events(StableRNG(1), design)
    @test size(evt, 1) == 12
    @test evt.event == ['S', 'C', 'R', 'S', 'C', 'R', 'S', 'C', 'R', 'S', 'C', 'R']


    # repeat first, then sequence => same sequence
    design = SingleSubjectDesign(conditions = Dict(:condition => ["A", "B"]))
    design = RepeatDesign(design, 2)
    design = SequenceDesign(design, "S[ABCD]")
    evt = generate_events(StableRNG(2), design)

    @test all(evt.event[2:2:end] .== 'B')


    # sequence first, then repeat => different sequence for each repetition
    design = SingleSubjectDesign(conditions = Dict(:condition => ["A", "B"]))
    design = SequenceDesign(design, "S[ABCD]")
    design = RepeatDesign(design, 2)
    evt = generate_events(StableRNG(2), design)
    @test !all(evt.event[2:2:end] .== 'B')
end

@testset "simulate_sequence" begin
    design = SingleSubjectDesign(conditions = Dict(:condition => ["one", "two"]))
    design = SequenceDesign(design, "SCR_")
    c = LinearModelComponent(;
        basis = UnfoldSim.hanning(40),
        formula = @formula(0 ~ 1 + condition),
        β = [1.0, 2.0],
        contrasts = Dict(:cond => EffectsCoding()),
    )
    s, e = simulate(design, c, NoOnset(); return_epoched = true)
    @test size(s) == (40, 6)
    @test s[:, 1] == s[:, 2] # If no component dict is specified, all events have the same component

end

@testset "rand_re" begin

    machine = Automa.compile(Automa.RegExp.RE("b+l+a+"))
    @test UnfoldSim.rand_re(MersenneTwister(2), machine) == "bbbbblllaa"
    # trivial single-character regex
    @test UnfoldSim.rand_re(MersenneTwister(1), Automa.compile(Automa.RegExp.RE("a"))) ==
          "a"

    # different seeds should (very likely) produce different strings
    r1 = UnfoldSim.rand_re(MersenneTwister(1), machine)
    r2 = UnfoldSim.rand_re(MersenneTwister(2), machine)
    @test r1 != r2

    # produced string should match the regex pattern
    for k = 1:10
        s = UnfoldSim.rand_re(MersenneTwister(k), machine)
        @test occursin(r"^b+l+a+$", s)
    end

    # integration: evaluate_sequencestring expands curly braces and delegates to rand_re
    @test UnfoldSim.evaluate_sequencestring(MersenneTwister(1), "bla{3,4}") == "blaaaa"
end