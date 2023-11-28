@testset "simulation" begin

    ## Define elements for the simulation

    # Define experimental factors
    conditions = Dict(:cond => ["A", "B"])

    # Create design for one subject with 20 trials (10 in each of the two factor levels)
    repetitions = 10
    design_single_subject = SingleSubjectDesign(;
        conditions = conditions) |> x -> RepeatDesign(x, repetitions)

    # Create design for multiple subjects with conditions as a between-items factor
    n_subjects = 5
    n_items = 4
    design_multiple_subjects = MultiSubjectDesign(;n_subjects=n_subjects, n_items=n_items, items_between = conditions)

    # Linear component for the single-subject simulation
    signal_linear = LinearModelComponent(;
        basis = p100(),
        formula = @formula(0 ~ 1 + cond),
        β = [1, 0.5])

    # Mixed-model component for the multi-subject simulation
    signal_mixed = MixedModelComponent(;
        basis = p100(),
        formula = @formula(0 ~ 1 + cond + (1 + cond|subject)),
        β = [1, 0.5],
        σs = Dict(:subject => [0.2,0.1]))

    # Define headmodel and MultichannelComponent for multi-channel cases
    hartmut_model = headmodel(type="hartmut")
    signal_linear_multichannel = MultichannelComponent(signal_linear, hartmut_model => "Left Central Opercular Cortex")
    signal_mixed_multichannel = MultichannelComponent(signal_mixed, hartmut_model => "Left Central Opercular Cortex")
    
    # Overlap since offset<length(signal.basis)
    onset = UniformOnset(; width=10, offset=5)

    noise = PinkNoise(; noiselevel=0.5)

    @testset "single_subject-single_channel" begin

        ## Simulate data
        simulation = Simulation(design_single_subject, signal_linear, onset, noise)
        data, events = simulate(MersenneTwister(42), simulation)

        ## Tests
        # Check whether the number of events is equal to the number of condition levels * repetitions of the design
        @test nrow(events) == length(conditions[:cond])*repetitions

        @test typeof(data) == Vector{Float64}
        # Check whether the length of the data is equal to the last event onset + the length of the basis
        @test size(data) == (maximum(events.latency)+length(signal_linear.basis),)
    end

    @testset "single_subject-multiple_channels" begin

        ## Simulate data
        data, events = simulate(MersenneTwister(42), design_single_subject, signal_linear_multichannel, onset, noise)

        ## Tests
        @test typeof(data) == Matrix{Float64}

        # Compute expected number of channels and expected length of the eeg signal
        n_channels_exp = length(hartmut_model.electrodes["label"])
        eeg_length_exp = maximum(events.latency) + length(signal_linear.basis)

        # channels x eeg
        @test size(data) == (n_channels_exp, eeg_length_exp)
        
    end

    @testset "multiple_subjects-single_channel" begin

        ## Simulate data
        data, events = simulate(MersenneTwister(42), design_multiple_subjects, signal_mixed, onset, noise)

        ## Tests
        # Compute expected length of the (longest) eeg signal
        eeg_length_exp = maximum(events.latency) + length(signal_mixed.basis)

        # eeg x subjects
        @test size(data) == (eeg_length_exp, n_subjects)
        @test typeof(data) == Matrix{Float64}
    end

    @testset "multiple_subjects-multiple_channels" begin
        
        ## Simulate data
        data, events = simulate(MersenneTwister(42), design_multiple_subjects, signal_mixed_multichannel, onset, noise)

        ## Tests
        # Compute expected number of channels and expected length of the eeg signal
        n_channels_exp = length(hartmut_model.electrodes["label"])
        eeg_length_exp = maximum(events.latency) + length(signal_linear.basis)

        # channels x eeg x subjects
        @test size(data) == (n_channels_exp, eeg_length_exp, n_subjects)
        @test typeof(data) == Array{Float64,3}
    end

end
