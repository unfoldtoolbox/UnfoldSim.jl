# # Overview: Onset types
# The onset types determine the distances between event onsets in the continuous EEG signal. The distances are sampled from a certain probability distribution.
# Currently, there are two types of onset distributions implemented: `UniformOnset` and `LogNormalOnset`.

# ### Setup
# ```@raw html
# <details>
# <summary>Click to expand</summary>
# ```
## Load required packages
using UnfoldSim
using CairoMakie
using Random

## Define a simple design and repeat it 10000.
## This will result in 20000 events i.e. event onsets.
design =
    SingleSubjectDesign(conditions = Dict(:cond => ["A", "B"])) |>
    x -> RepeatDesign(x, 10000);

# ```@raw html
# </details >
# ```

# ## UniformOnset
# The `UniformOnset` is based on a uniform distribution and has two parameters: `width` and `offset`.

# Example:
onset_uniform = UniformOnset(; width = 50, offset = 0);

# The `width` parameter defines the upper bound of the interval of the uniform distribution (its lower bound is 0) i.e. all values between 0 and `width` are equally probable.

# The `offset` parameter determines the minimal distance between two events and its value is added to the value sampled from the uniform distribution i.e. it shifts the distribution.
# Its default value is `0`, i.e. no offset.

# In the figure below, it is illustrated how the onset distribution changes when changing one of its parameters.
let # hide
    f = Figure() # hide

    ## Define parameter combinations # hide
    parameters = [ # hide
        (((50, 0), (80, 0)), "width"), # hide
        (((50, 0), (50, 20)), "offset"), # hide
    ] # hide

    axes_list = Array{Any}(undef, length(parameters)) # hide

    ## Create a subplot for each parameter i.e. one for width and one for offset # hide
    for (index, (combinations, label)) in enumerate(parameters) # hide
        ax = Axis(f[index, 1], title = "Parameter: $label") # hide
        axes_list[index] = ax # hide

        ## Go through all parameter combinations and plot a histogram of the sampled onsets # hide
        for (width, offset) in combinations # hide
            distances = UnfoldSim.simulate_interonset_distances( # hide
                MersenneTwister(42), # hide
                UniformOnset(; width = width, offset = offset), # hide
                design, # hide
            ) # hide

            hist!( # hide
                ax, # hide
                distances, # hide
                bins = range(0, 100, step = 1), # hide
                label = "($width, $offset)", # hide
            ) # hide

            if label == "offset" && offset != 0 # hide 
                vlines!(offset, color = "black") # hide
            end # hide
        end # hide
        hideydecorations!(ax) # hide
        hidespines!(ax, :t, :r) # hide
        axislegend( # hide
            ax, # hide
            framevisible = false, # hide
            labelsize = 12, # hide
            markersize = 5, # hide
            patchsize = (10, 10), # hide
        ) # hide
    end # hide
    axes_list[end].xlabel = "Time between events [samples]" # hide
    linkyaxes!(axes_list...) # hide
    current_figure() # hide
end # hide

## Note: The code is repeated because I did not manage to show the figure but make the code collapsible # hide
# ```@raw html
# <details>
# <summary>Click to show the code for the figure above</summary>
# ```
let
    f = Figure()

    ## Define parameter combinations
    parameters = [(((50, 0), (80, 0)), "width"), (((50, 0), (50, 20)), "offset")]

    axes_list = Array{Any}(undef, length(parameters))

    ## Create a subplot for each parameter i.e. one for width and one for offset
    for (index, (combinations, label)) in enumerate(parameters)
        ax = Axis(f[index, 1], title = "Parameter: $label")
        axes_list[index] = ax

        ## Go through all parameter combinations and plot a histogram of the sampled onsets
        for (width, offset) in combinations
            onsets = UnfoldSim.simulate_interonset_distances(
                MersenneTwister(42),
                UniformOnset(; width = width, offset = offset),
                design,
            )

            hist!(ax, onsets, bins = range(0, 100, step = 1), label = "($width, $offset)")

            if label == "offset" && offset != 0
                vlines!(offset, color = "black")
            end
        end
        hideydecorations!(ax)
        hidespines!(ax, :t, :r)
        axislegend(
            ax,
            framevisible = false,
            labelsize = 12,
            markersize = 5,
            patchsize = (10, 10),
        )
    end
    axes_list[end].xlabel = "Time between events [samples]"
    linkyaxes!(axes_list...)
end
# ```@raw html
# </details >
# ```

# ## LogNormalOnset
# The `LogNormalOnset` is based on a log-normal distribution and has four parameters: `μ`, `σ`, `offset` and `truncate_upper`. 

# Example:
onset_lognormal = LogNormalOnset(; μ = 3, σ = 0.25, offset = 0, truncate_upper = nothing);

# The parameters `μ` and `σ` are the location and scale parameter of the log-normal distribution. However, they are not identical to its mean and standard deviation.
# If a variable $X$ is log-normally distributed then $Y = ln(X)$ is normally distributed with mean `μ` and standard deviation `σ`[^1].

# The `offset` parameter determines the minimal distance between two events and its value is added to the value sampled from the log-normal distribution i.e. it shifts the distribution.
# Its default value is `0`, i.e. no offset.

# The `truncate_upper` parameter allows to truncate the distribution at a certain sample value. Its default value is `nothing`, i.e. no truncation.

# In the figure below, it is illustrated how the onset distribution changes when changing one of its parameters.
let # hide
    f = Figure(size = (600, 800)) # hide

    ## Define parameter combinations # hide
    parameters = [ # hide
        (((3, 0.25, 0, nothing), (2.5, 0.25, 0, nothing)), "μ"), # hide
        (((3, 0.25, 0, nothing), (3, 0.35, 0, nothing)), "σ"), # hide
        (((3, 0.25, 0, nothing), (3, 0.25, 30, nothing)), "offset"), # hide
        (((3, 0.25, 0, nothing), (3, 0.25, 0, 25)), "truncate_upper"), # hide
    ] # hide

    axes_list = Array{Any}(undef, length(parameters)) # hide

    ## Create a subplot for each parameter i.e. one for μ, one for σ etc # hide
    for (index, (combinations, label)) in enumerate(parameters) # hide
        ax = Axis(f[index, 1], title = "Parameter: $label") # hide
        axes_list[index] = ax # hide

        ## Go through all parameter combinations and plot a histogram of the sampled onsets # hide
        for (μ, σ, offset, truncate_upper) in combinations # hide
            onsets = UnfoldSim.simulate_interonset_distances( # hide
                MersenneTwister(42), # hide
                LogNormalOnset(; # hide
                    μ = μ, # hide
                    σ = σ, # hide
                    offset = offset, # hide
                    truncate_upper = truncate_upper, # hide
                ), # hide
                design, # hide
            ) # hide

            hist!( # hide
                ax, # hide
                onsets, # hide
                bins = range(0, 100, step = 1), # hide
                label = "($μ,$σ,$offset,$truncate_upper)", # hide
            ) # hide

            if label == "offset" && offset !== 0 # hide
                vlines!(offset, color = "black") # hide
            elseif label == "truncate_upper" && truncate_upper !== nothing # hide
                vlines!(truncate_upper, color = "black") # hide
            end # hide
        end # hide
        hideydecorations!(ax) # hide
        hidespines!(ax, :t, :r) # hide
        axislegend( # hide 
            ax, # hide
            framevisible = false, # hide
            labelsize = 12, # hide
            markersize = 5, # hide
            patchsize = (10, 10), # hide
        ) # hide
    end # hide
    axes_list[end].xlabel = "Time between events [samples]" # hide
    linkyaxes!(axes_list...) # hide
    current_figure() # hide
end # hide


## Note: The code is repeated because I did not manage to show the figure but make the code collapsible # hide
# ```@raw html
# <details>
# <summary>Click to show the code for the figure above</summary>
# ```
let
    f = Figure(size = (600, 800))

    ## Define parameter combinations
    parameters = [
        (((3, 0.25, 0, nothing), (2.5, 0.25, 0, nothing)), "μ"),
        (((3, 0.25, 0, nothing), (3, 0.35, 0, nothing)), "σ"),
        (((3, 0.25, 0, nothing), (3, 0.25, 30, nothing)), "offset"),
        (((3, 0.25, 0, nothing), (3, 0.25, 0, 25)), "truncate_upper"),
    ]

    axes_list = Array{Any}(undef, length(parameters))

    ## Create a subplot for each parameter i.e. one for μ, one for σ etc
    for (index, (combinations, label)) in enumerate(parameters)
        ax = Axis(f[index, 1], title = "Parameter: $label")
        axes_list[index] = ax

        ## Go through all parameter combinations and plot a histogram of the sampled onsets
        for (μ, σ, offset, truncate_upper) in combinations
            onsets = UnfoldSim.simulate_interonset_distances(
                MersenneTwister(42),
                LogNormalOnset(;
                    μ = μ,
                    σ = σ,
                    offset = offset,
                    truncate_upper = truncate_upper,
                ),
                design,
            )

            hist!(
                ax,
                onsets,
                bins = range(0, 100, step = 1),
                label = "($μ,$σ,$offset,$truncate_upper)",
            )

            if label == "offset" && offset !== 0
                vlines!(offset, color = "black")
            elseif label == "truncate_upper" && truncate_upper !== nothing
                vlines!(truncate_upper, color = "black")
            end
        end
        hideydecorations!(ax)
        hidespines!(ax, :t, :r)
        axislegend(
            ax,
            framevisible = false,
            labelsize = 12,
            markersize = 5,
            patchsize = (10, 10),
        )
    end
    axes_list[end].xlabel = "Time between events [samples]"
    linkyaxes!(axes_list...)
end
# ```@raw html
# </details >
# ```

# # Overlap of subsequent events
# !!! note
#       The overlap of subsequent events can be indirectly controlled by setting the `offset` parameter relative to the length of the component basis.
#       Assuming that `signal` is a component e.g. `LinearModelComponent`,
#        - if `offset` > `length(signal.basis)` -> no overlap
#        - if `offset` < `length(signal.basis)` -> there might be overlap, depending on the other parameters of the onset distribution

# [^1]: Wikipedia contributors. (2023, December 5). Log-normal distribution. In Wikipedia, The Free Encyclopedia. Retrieved 12:27, December 7, 2023, from https://en.wikipedia.org/w/index.php?title=Log-normal_distribution&oldid=1188400077# 
