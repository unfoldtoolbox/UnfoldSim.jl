using Base: AbstractCartesianIndex
@testset "simulation" begin

    @testset "general_test_simulate" begin
        ## Define elements for the simulation

        # Define experimental factors
        conditions = Dict(:cond => ["A", "B"])

        # Create design for one subject with 20 trials (10 in each of the two factor levels)
        repetitions = 10
        design_single_subject =
            SingleSubjectDesign(; conditions = conditions) |>
            x -> RepeatDesign(x, repetitions)

        # Create design for multiple subjects with conditions as a between-items factor
        n_subjects = 5
        n_items = 4
        design_multiple_subjects = MultiSubjectDesign(;
            n_subjects = n_subjects,
            n_items = n_items,
            items_between = conditions,
        )

        # Linear component for the single-subject simulation
        signal_linear = LinearModelComponent(;
            basis = p100(),
            formula = @formula(0 ~ 1 + cond),
            β = [1, 0.5],
        )

        # Mixed-model component for the multi-subject simulation
        signal_mixed = MixedModelComponent(;
            basis = p100(),
            formula = @formula(0 ~ 1 + cond + (1 + cond | subject)),
            β = [1, 0.5],
            σs = Dict(:subject => [0.2, 0.1]),
        )

        # Define headmodel and MultichannelComponent for multi-channel cases
        hartmut_model = Hartmut()
        signal_linear_multichannel = MultichannelComponent(
            signal_linear,
            hartmut_model => "Left Central Opercular Cortex",
        )
        signal_mixed_multichannel = MultichannelComponent(
            signal_mixed,
            hartmut_model => "Left Central Opercular Cortex",
        )

        # Overlap since offset<length(signal.basis)
        onset = UniformOnset(; width = 10, offset = 5)

        noise = PinkNoise(; noiselevel = 0.5)

        @testset "single_subject-single_channel" begin

            ## Simulate data
            simulation = Simulation(design_single_subject, signal_linear, onset, noise)
            data, events = simulate(MersenneTwister(42), simulation)

            ## Tests
            # Check whether the number of events is equal to the number of condition levels * repetitions of the design
            @test nrow(events) == length(conditions[:cond]) * repetitions

            @test typeof(data) == Vector{Float64}
            # Check whether the length of the data is equal to the last event onset + the length of the basis
            @test size(data) == (maximum(events.latency) + length(signal_linear.basis),)
        end

        @testset "single_subject-multiple_channels" begin

            ## Simulate data
            data, events = simulate(
                MersenneTwister(42),
                design_single_subject,
                signal_linear_multichannel,
                onset,
                noise,
            )

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
            data, events = simulate(
                MersenneTwister(42),
                design_multiple_subjects,
                signal_mixed,
                onset,
                noise,
            )

            ## Tests
            # Compute expected length of the (longest) eeg signal
            eeg_length_exp = maximum(events.latency) + length(signal_mixed.basis)

            # eeg x subjects
            @test size(data) == (eeg_length_exp, n_subjects)
            @test typeof(data) == Matrix{Float64}
        end

        @testset "multiple_subjects-multiple_channels" begin

            ## Simulate data
            data, events = simulate(
                MersenneTwister(42),
                design_multiple_subjects,
                signal_mixed_multichannel,
                onset,
                noise,
            )

            ## Tests
            # Compute expected number of channels and expected length of the eeg signal
            n_channels_exp = length(hartmut_model.electrodes["label"])
            eeg_length_exp = maximum(events.latency) + length(signal_linear.basis)

            # channels x eeg x subjects
            @test size(data) == (n_channels_exp, eeg_length_exp, n_subjects)
            @test typeof(data) == Array{Float64,3}
        end
    end

    # Test the output dimensions for all combinations of subject, channel, onset and return_epoched
    @testset "output dimenions" begin
        for subject ∈ ["single", "multi"]
            if subject == "multi"
                design = MultiSubjectDesign(;
                    n_subjects = 3,
                    n_items = 7,
                    subjects_between = Dict(:cond => nlevels(3, 'C')),
                )

                comp = MixedModelComponent(;
                    basis = [1, 2, 3, 4],
                    formula = @formula(0 ~ 1 + (1 | subject)),
                    β = [1],
                    σs = Dict(:subject => [1]),
                )
            else
                design =
                    SingleSubjectDesign(; conditions = Dict(:a => ["a"])) |>
                    x -> RepeatDesign(x, 7)
                comp = LinearModelComponent(;
                    basis = [1, 2, 3, 4],
                    formula = @formula(0 ~ 1),
                    β = [1],
                )

            end
            for channel ∈ ["single", "multi"]
                if channel == "multi"
                    comp = comp |> x -> MultichannelComponent(x, [1, 2, -3, 4, 5])
                end

                for sim_onset ∈ ["noonset", "yesonset"]
                    if sim_onset == "yesonset"
                        onset = UniformOnset(; width = 20, offset = 4)
                    else
                        onset = NoOnset()
                    end

                    for return_epoched ∈ [false, true]
                        simulation = Simulation(design, comp, onset, NoNoise())

                        local data
                        try
                            data, events = simulate(
                                MersenneTwister(1),
                                simulation;
                                return_epoched = return_epoched,
                            )
                        catch AssertionError
                            # The AssertionError in the simulate function should be elicited only in the case below
                            @test (sim_onset == "noonset") & (return_epoched == false)

                        else
                            sz = size(data)

                            if channel == "multi"
                                @test sz[1] == 5
                                offset = 1
                            else
                                @test all(sz .!= 5)
                                offset = 0
                            end
                            if return_epoched == true
                                @test sz[1+offset] == 4
                                @test sz[2+offset] == 7
                                if subject == "multi"
                                    @test sz[3+offset] == 3
                                end
                            else
                                @test sz[1+offset] > 4
                                if subject == "multi"
                                    @test sz[2+offset] == 3
                                end
                            end
                        end
                    end
                end

            end
        end
    end

    @testset "test data-type" begin
        # Define experimental factors
        conditions = Dict(:cond => ["A", "B"])

        # Create design for one subject with 20 trials (10 in each of the two factor levels)
        repetitions = 10
        design_single_subject =
            SingleSubjectDesign(; conditions = conditions) |>
            x -> RepeatDesign(x, repetitions)

        # Linear component for the single-subject simulation
        signal_linear = LinearModelComponent(;
            basis = p100(),
            formula = @formula(0 ~ 1 + cond),
            β = [1, 0.5],
        )


        # Define headmodel and MultichannelComponent for multi-channel cases
        hartmut_model = Hartmut()
        signal_linear_multichannel = MultichannelComponent(
            signal_linear,
            hartmut_model => "Left Central Opercular Cortex",
        )

        # Overlap since offset<length(signal.basis)
        onset = UniformOnset(; width = 10, offset = 5)

        noise = PinkNoise(; noiselevel = 0.5)


        ## Simulate data
        simulation = Simulation(design_single_subject, signal_linear, onset, noise)
        simulation_c64 =
            Simulation{Complex}(design_single_subject, [signal_linear], onset, noise)
        @test typeof(simulation) == Simulation{Float64}
        @test typeof(simulation_c64) == Simulation{Complex}
        data, events = simulate(MersenneTwister(42), simulation)
        @test eltype(data) == Float64

        data, events = simulate(MersenneTwister(42), simulation_c64)
        @test eltype(data) == Complex



    end
    @testset "multi-component sequence #124" begin
        struct MyLinearModelComponent1 <: AbstractComponent
            comp::Any
        end
        MyLinearModelComponent1(b, f, β) =
            MyLinearModelComponent1(LinearModelComponent(; basis = b, formula = f, β))
        UnfoldSim.simulate_component(
            rng,
            c::MyLinearModelComponent1,
            design::UnfoldSim.SubselectDesign,
        ) = simulate_component(rng, c.comp, design)
        UnfoldSim.length(c::MyLinearModelComponent1) = length(c.comp)
        UnfoldSim.size(c::MyLinearModelComponent1) = size(c.comp)
        sim = Simulation(
            SingleSubjectDesign(conditions = Dict(:event => ['A', 'B'])),
            Dict(
                'A' => [
                    LinearModelComponent(
                        basis = p100(),
                        formula = @formula(0 ~ 1),
                        β = [1],
                    ),
                ],
                'B' => [MyLinearModelComponent1(p100(), @formula(0 ~ 1), [2])],
            ),
            NoOnset(),
            NoNoise(),
        )
        d, e = simulate(UnfoldSim.MersenneTwister(1), sim; return_epoched = true)
        @test d[10, 1] > 0.9 # 1 if the hanning would hit perfectly (currently the peak is between samples)
        @test d[10, 2] > 1.9 # 2 if the hanning would hit perfectly (currently the peak is between samples)
    end

end
