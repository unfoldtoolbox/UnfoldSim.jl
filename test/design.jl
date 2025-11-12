@testset "design" begin

    @testset "SingleSubjectDesign" begin
        des = SingleSubjectDesign(; conditions = Dict(:A => nlevels(5), :B => nlevels(2)))
        @test generate_events(des) == generate_events(des)
        @test size(generate_events(des)) == (5 * 2, 2)
        @test names(generate_events(des)) == ["A", "B"]
        @test sum(generate_events(des).A .== "S2") == 2
        @test sum(generate_events(des).B .== "S2") == 5
        @test sum(
            (generate_events(des).B .== "S2") .&& (generate_events(des).A .== "S2"),
        ) == 1
        des = SingleSubjectDesign(;
            conditions = Dict(:A => nlevels(5), :B => nlevels(2)),
            event_order_function = (rng, x) -> sort(x, order(:B, rev = true)),
        )
        @test generate_events(des).B[1] == "S2"

        des = SingleSubjectDesign(;
            conditions = Dict(:A => nlevels(5), :B => nlevels(2)),
            event_order_function = (rng, x) -> shuffle(rng, x),
        )
        @test generate_events(MersenneTwister(3), des) ==
              generate_events(MersenneTwister(3), des)

        # different sortig seed results in different sorts
        @test generate_events(MersenneTwister(3), des) !=
              generate_events(MersenneTwister(4), des)
    end

    @testset "MultiSubjectDesign" begin
        des = MultiSubjectDesign(;
            n_subjects = 10,
            n_items = 100,
            both_within = Dict(:A => nlevels(5), :B => nlevels(2)),
        )


        @test size(generate_events(des)) == (10 * 100 * 5 * 2, 4)
        @test names(generate_events(des)) == ["subject", "item", "A", "B"]
        @test sum(generate_events(des).A .== "S2") == 10 * 100 * 2
        @test sum(generate_events(des).B .== "S2") == 10 * 100 * 5

        # check that finally everything is sorted by subject
        des = MultiSubjectDesign(;
            n_subjects = 10,
            n_items = 100,
            both_within = Dict(:A => nlevels(5), :B => nlevels(2)),
            event_order_function = (rng, x) -> sort(x, order(:item, rev = true)),
        )
        @test generate_events(des).subject[1] == "S01"

        # check event_order_function
        des = MultiSubjectDesign(;
            n_subjects = 10,
            n_items = 100,
            both_within = Dict(:A => nlevels(5), :B => nlevels(2)),
            event_order_function = (rng, x) -> sort(x, order(:B, rev = true)),
        )
        @test generate_events(des).B[1] == "S2"

        # check event_order_function
        des = MultiSubjectDesign(;
            n_subjects = 10,
            n_items = 100,
            both_within = Dict(:A => nlevels(5), :B => nlevels(2)),
            event_order_function = (rng, x) -> shuffle(rng, x),
        )
        # generating same events with same seed should result in same events
        @test generate_events(MersenneTwister(3), des) ==
              generate_events(MersenneTwister(3), des)

        # different sortig seed results in different sorts
        @test generate_events(MersenneTwister(3), des) !=
              generate_events(MersenneTwister(4), des)

        # check that this throws an error because of `dv` as condition name
        des = MultiSubjectDesign(;
            n_subjects = 10,
            n_items = 100,
            both_within = Dict(:dv => nlevels(5), :B => nlevels(2)),
        )
        @test_throws AssertionError generate_events(des)
    end

    @testset "RepeatDesign" begin

        designOnce = MultiSubjectDesign(;
            n_items = 8,
            n_subjects = 12,
            subjects_between = Dict(:cond => ["levelA", "levelB"]),
            items_between = Dict(:cond => ["levelA", "levelB"]),
        )

        design = RepeatDesign(designOnce, 3)
        # Note: the number of items has to be divided by the number of cond levels because cond is both between item and subject
        @test size(generate_events(design)) == ((8 / 2) * 12 * 3, 3)
        @test size(design) == (8 * 3, 12)
        #--- single sub

        designOnce =
            SingleSubjectDesign(; conditions = Dict(:A => nlevels(5), :B => nlevels(2)))

        design = RepeatDesign(designOnce, 3)
        @test size(generate_events(design)) == (5 * 2 * 3, 2)
        @test size(design) == (3 * 5 * 2,)
    end

    @testset "MultiSubjectDesign without conditions" begin
        # Define number of subjects and items
        n_subjects = 3
        n_items = 5

        # Create a multi-subject design without any conditions
        design_multi = MultiSubjectDesign(; n_subjects = n_subjects, n_items = n_items)

        # Create the events data frame based on the design defined above
        events_multi = generate_events(design_multi)

        # Test that there are only subject and item columns (no condition columns) in the events df
        @test names(events_multi) == ["subject", "item"]

        # Test that the total length of the events df is the number of subjects times
        # the number of items since every item occurs in every subject
        @test size(events_multi, 1) == n_subjects * n_items

        # Test that all subjects and items appear in the events df (and also not more)
        @test length(unique(events_multi.subject)) == n_subjects
        @test length(unique(events_multi.item)) == n_items

        # Compute number of items per subject
        n_items_per_subject = combine(groupby(events_multi, :subject), :item => length)
        # Check that all subjects have the same number of items (and this number equals n_items)
        @test all(.==(n_items, n_items_per_subject.item_length))

        # Simulate some data with a multi-subject design without conditions

        # Define a component that only has a 1 as its basis to facilitate counting the peaks
        component = MixedModelComponent(;
            basis = [1],
            formula = @formula(0 ~ 1 + (1 | subject)),
            β = [5],
            σs = Dict(:subject => [1]),
        )

        # Define an onset without Overlap
        onset = UniformOnset(; width = 10, offset = 3)

        # Simulate data without noise
        data, events = simulate(StableRNG(1), design_multi, component, onset, NoNoise())

        # Check that the simulated data has as many columns as subjects
        @test size(data, 2) == n_subjects

        # Check that (for each subject) the data has as many peaks (i.e. events) as set in n_items
        for s = 1:n_subjects
            @test count(x -> .!iszero(x), data[:, s]) == n_items
        end
    end

    @testset "MultiSubjectDesign with between subject and item factors" begin

        # Define number of subjects and items
        n_subjects = 2
        n_items = 2

        # Create a multi-subject design with the same factor between subject and item
        design = MultiSubjectDesign(;
            n_subjects = n_subjects,
            n_items = n_items,
            subjects_between = Dict(:cond => ["levelA", "levelB"]),
            items_between = Dict(:cond => ["levelA", "levelB"]),
        )

        # Create events data frame based on the design
        events_df = generate_events(design)

        # Extract events for subject 2
        s2 = subset(events_df, :subject => ByRow(==("S2")))

        # Since cond is a between_subject factor, each subject should only be in one condition
        @test length(unique(s2.cond)) == 1

        # Extract events for item 1
        i1 = subset(events_df, :item => ByRow(==("I1")))

        # Since cond is a between_item factor, each item should only be in one condition
        @test length(unique(i1.cond)) == 1
    end

    @testset "Effects Design" begin
        # Begin with simulation design
        design = SingleSubjectDesign(;
            conditions = Dict(
                :condition => ["car", "face"],
                :continuous => range(0, 5, length = 10),
            ),
        )

        # Effects dictionary
        effects_dict_1 = Dict(:condition => ["car", "face"])
        effects_dict_2 = Dict(:condition => ["car", "face"], :continuous => [2, 3, 4])

        # Generate effects design
        ef_design_1 = EffectsDesign(design, effects_dict_1)
        ef_design_2 = EffectsDesign(design, effects_dict_2)

        # Generate events
        ef_events_1 = generate_events(ef_design_1)
        ef_events_2 = generate_events(ef_design_2)

        # SingleSubject tests
        @test size(ef_events_1, 1) == 2 # Test correct length of events df
        @test unique(ef_events_1[!, :continuous])[1] ≈ mean(range(0, 5, length = 10)) # Test that average is calculated correctly and only one value is present in df
        @test size(ef_events_2, 1) == 6 # Test correct length of events df when one inputs values for continuous variable

        # MultiSubjectDesign -> not implemented yet, so should error
        design = MultiSubjectDesign(
            n_subjects = 20,
            n_items = 8,
            items_between = Dict(:condition => ["car", "face"], :continuous => [1, 2]),
        )
        @test_throws ErrorException EffectsDesign(design, effects_dict_1)

    end
end
