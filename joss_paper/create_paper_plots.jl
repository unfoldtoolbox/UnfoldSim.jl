using UnfoldSim
using Unfold
using UnfoldMakie
using CairoMakie
using DataFrames
using PrettyTables
using StableRNGs
using DSP
using Random
using ProjectRoot

# If the variable has not been set by the Github action then it is set in the script
if !@isdefined save_plots
    # Determines whether the plots should be saved or not
    save_plots = false
end


# Check whether the plot directory already exists otherwise create it
if !isdir(@projectroot("plots"))
    mkdir(@projectroot("plots"))
end

#----
# Simulation example: Custom continuous EEG signal
let
    # Simulating data using UnfoldSim
    design =
        SingleSubjectDesign(;
            conditions = Dict(
                :condition => ["car", "face"],
                :continuous => range(0, 5, length = 10),
            ),
            event_order_function = x -> shuffle(deepcopy(StableRNG(1)), x),
        ) |> x -> RepeatDesign(x, 100)

    n1 = LinearModelComponent(;
        basis = n170(),
        formula = @formula(0 ~ 1 + condition),
        β = [5, 3],
    )

    p3 = LinearModelComponent(;
        basis = p300(),
        formula = @formula(0 ~ 1 + continuous + continuous^2),
        β = [5, 1, 0.2],
    )

    components = [n1, p3]

    onset = UniformOnset(; width = 0, offset = 200)

    noise = PinkNoise(; noiselevel = 2)

    eeg_data, events_df = simulate(StableRNG(1), design, components, onset, noise)

    # Analysing data using Unfold
    m = fit(
        UnfoldModel,
        Dict(
            Any => (
                @formula(0 ~ 1 + condition + spl(continuous, 4)),
                firbasis(τ = [-0.1, 1], sfreq = 100, name = "basis"),
            ),
        ),
        events_df,
        eeg_data,
    )
    #----
    # Visualise event data frame
    pretty_table(
        events_df[1:5, :],
        alignment = :l,
        backend = Val(:markdown),
        header = names(events_df),
    )
    #----
    # Visualise simulated data
    f = Figure(size = (1000, 400))
    ax = Axis(
        f[1, 1],
        title = "Simulated EEG data",
        titlesize = 18,
        xlabel = "Time [samples]",
        ylabel = "Amplitude [µV]",
        xlabelsize = 16,
        ylabelsize = 16,
        xgridvisible = false,
        ygridvisible = false,
    )

    n_samples = 1400
    lines!(eeg_data[1:n_samples]; color = "black")
    v_lines = [
        vlines!(
            [r["latency"]];
            color = ["orange", "teal"][1+(r["condition"]=="car")],
            label = r["condition"],
        ) for r in
        filter(:latency => x -> x < n_samples, events_df)[:, ["latency", "condition"]] |> eachrow
    ]
    xlims!(ax, 0, n_samples)
    axislegend("Event onset"; unique = true)
    #f[1,2] = Legend(f,ax,"Event onset", unique =true)

    if save_plots
        save(@projectroot("plots", "example_simulated_data.svg"), f, pt_per_unit = 1)
    end
    current_figure()

    #----
    # Create a data frame with the model coefficients and extract the coefficient names
    coefs = coeftable(m)
    coefnames = unique(coefs.coefname)

    f = Figure(size = (1000, 400))
    ga = f[1, 1] = GridLayout()
    gb = f[1, 2] = GridLayout()

    # Plot A: Estimated regression parameters
    ax_A = Axis(
        ga[1, 1],
        title = "Estimated regression parameters",
        titlegap = 12,
        xlabel = "Time [s]",
        ylabel = "Amplitude [μV]",
        xlabelsize = 16,
        ylabelsize = 16,
        xgridvisible = false,
        ygridvisible = false,
    )

    for coef in coefnames
        estimate = filter(:coefname => ==(coef), coefs)

        lines!(ax_A, estimate.time, estimate.estimate, label = coef)
    end
    axislegend("Coefficient", framevisible = false)
    hidespines!(ax_A, :t, :r)

    # Plot B: Marginal effects
    plot_B = plot_erp!(
        gb,
        effects(Dict(:condition => ["car", "face"], :continuous => 0:0.5:5), m);
        mapping = (; color = :continuous, linestyle = :condition, group = :continuous),
        legend = (; valign = :top, halign = :right, tellwidth = false),
        categorical_color = false,
        axis = (
            title = "Marginal effects",
            titlegap = 12,
            xlabel = "Time [s]",
            ylabel = "Amplitude [μV]",
            xlabelsize = 16,
            ylabelsize = 16,
            xgridvisible = false,
            ygridvisible = false,
        ),
        layout = (; showlegend = false),
    )

    # Workaround to separate legend and colorbar (will be fixed in a future UnfoldMakie version)
    legend = plot_B.content[2].content
    plot_B[:, 1] = legend


    # Add letter labels to the plots
    # adapted from: https://docs.makie.org/stable/tutorials/layout-tutorial/
    for (label, layout) in zip(["A", "B"], [ga, gb])
        Label(
            layout[1, 1, TopLeft()],
            label,
            fontsize = 26,
            font = :bold,
            padding = (0, 5, 5, 0),
            halign = :right,
        )
    end

    if save_plots
        save(@projectroot("plots", "example_coefficients_effects.svg"), f, pt_per_unit = 1)
    end
    current_figure()
end
#----
# Onset distribution plot
let
    # Create a simple design with many trials
    design =
        SingleSubjectDesign(conditions = Dict(:cond => ["A", "B"])) |>
        x -> RepeatDesign(x, 10000)

    # Define parameter combinations which should be plotted
    parameter_combinations = Dict(
        "Uniform Distribution" => [(70, 0), (60, 20)],
        "Lognormal Distribution" =>
            [(3, 0.35, 0, nothing), (3, 0.35, 20, 25), (3, 0.45, 40, nothing)],
    )

    f = Figure(size = (800, 400))
    axes_list = Array{Any}(undef, length(parameter_combinations))

    for (index, (distribution, parameters)) in enumerate(parameter_combinations)
        ax = Axis(
            f[index, 1],
            title = distribution,
            xgridvisible = false,
            ygridvisible = false,
            xlabelsize = 16,
        )
        axes_list[index] = ax

        for combination in parameters
            if distribution == "Uniform Distribution"
                onset = UniformOnset(combination...)
                global legend_title = "(width,offset)"
            else # Lognormal Distribution
                onset = LogNormalOnset(combination...)
                global legend_title = "(μ,σ,offset,truncate_upper)"
            end

            onsets = UnfoldSim.simulate_interonset_distances(StableRNG(1), onset, design)
            hist!(ax, onsets, bins = range(0, 100, step = 1), label = "$combination")

        end
        hideydecorations!(ax)
        hidespines!(ax, :t, :r)
        axislegend(
            legend_title,
            framevisible = false,
            labelsize = 12,
            markersize = 5,
            patchsize = (10, 10),
            titlesize = 12,
        )
    end
    axes_list[end].xlabel = "Time between events [samples]"

    if save_plots
        save(@projectroot("plots", "onset_distributions.svg"), f, pt_per_unit = 1)
    end

    current_figure()
end
#----
# Noise type plot
let
    f = Figure(size = (1000, 400))
    ga = f[1, 1] = GridLayout()
    gb = f[1, 2] = GridLayout()

    ax_A =
        ga[1, 1] = Axis(
            f;
            title = "Noise samples",
            titlesize = 18,
            xlabel = "Time",
            ylabel = "Amplitude",
            xlabelsize = 16,
            ylabelsize = 16,
            xgridvisible = false,
            ygridvisible = false,
        )

    ax_B =
        gb[1, 1] = Axis(
            f;
            title = "Welch Periodogram",
            titlesize = 18,
            xlabel = "Normalized frequency",
            ylabel = "log(Power)",
            xlabelsize = 16,
            ylabelsize = 16,
            xgridvisible = false,
            ygridvisible = false,
        )

    for n in [PinkNoise RedNoise WhiteNoise NoNoise ExponentialNoise]

        # Generate noise samples
        noisevec = simulate_noise(StableRNG(1), n(), 10000)

        # Plot 1000 samples
        lines!(ax_A, noisevec[1:1000]; label = string(n))

        # Calculate Welch periodogram
        perio = welch_pgram(noisevec)

        # Plot Welch periodogram
        lines!(ax_B, freq(perio), log10.(power(perio)), label = string(n))
    end

    # Add legend
    axislegend(framevisible = false, position = :rt)
    #f[1, 3] = Legend(f, ax_A, "Noise type", framevisible = false)

    # Make the first column (noise samples plot) take up 55% of the available space
    colsize!(f.layout, 1, Relative(0.60))
    #colsize!(f.layout, 2, Relative(0.28))

    # Add letter labels to the plots
    # adapted from: https://docs.makie.org/stable/tutorials/layout-tutorial/
    for (label, layout) in zip(["A", "B"], [ga, gb])
        Label(
            layout[1, 1, TopLeft()],
            label,
            fontsize = 26,
            font = :bold,
            padding = (0, 5, 5, 0),
            halign = :right,
        )
    end

    if save_plots
        save(@projectroot("plots", "noise_types.svg"), f, pt_per_unit = 1)
    end
    current_figure()
end
#----
