---
title: 'UnfoldSim.jl: Simulating continuous event-based time series data for EEG and beyond'
tags:
  - Julia
  - EEG
  - ERPs
  - evoked potentials
  - neuroimaging
  - simulation
  - time-series
  - regression ERPs
authors:
  - name: Judith Schepers
    orcid:  0009-0000-9270-730X
    equal-contrib: false
    affiliation: "1"
  - name: Luis Lips
    equal-contrib: false
    affiliation: "1"
  - name: Maanik Marathe
    equal-contrib: false
    affiliation: "1"
  - name: Benedikt V. Ehinger
    orcid:  0000-0002-6276-3332
    equal-contrib: false
    affiliation: "1, 2"
affiliations:
  - name: Institute for Visualisation and Interactive Systems, University of Stuttgart, Germany
    index: 1
  - name: Stuttgart Center for Simulation Science, University of Stuttgart, Germany
    index: 2
date: 06 November 2024
bibliography: paper.bib
---

# Summary

`UnfoldSim.jl` is a Julia package for simulating multivariate time series, with a focus on EEG, especially event-related potentials (ERPs). The user provides four ingredients: 1) an experimental design, with both categorical and continuous variables, 2) event basis functions specified via linear or hierarchical models, 3) an inter-event onset distribution, and 4) a noise specification. `UnfoldSim.jl` then simulates continuous EEG signals with potentially overlapping events. Multi-channel support via EEG-forward models is available as well. `UnfoldSim.jl` is modular, allowing users to implement custom designs, components, onset distributions or noise types to tailor the package to their needs. This allows support even for other modalities, e.g. single-voxel fMRI or pupil dilation signals.

One can find a detailed example of how to use `UnfoldSim.jl` to simulate continuous EEG data in the [documentation](https://unfoldtoolbox.github.io/UnfoldSim.jl/stable/generated/tutorials/simulateERP/).

# Statement of Need
In our work (e.g. @ehinger2019unfold, @dimigen2021regression), we often analyze data containing (temporally) overlapping events (e.g. stimulus onset and button press, or consecutive eye-fixations), non-linear effects, and complex experimental designs. For a multitude of reasons, we often need to simulate such kind of data: Simulated EEG data is useful to test preprocessing and analysis tools, validate statistical methods, illustrate conceptual issues, test toolbox functionalities, and find limitations of traditional analysis workflows. For instance, such simulation tools allow for testing the assumptions of new analysis algorithms and testing their robustness against any violation of these assumptions.

While other EEG simulation toolboxes exist, they each have limitations: they are dominantly based on the proprietary MATLAB software (The MathWorks Inc., Natick, USA), they do not simulate continuous EEG, and they offer little support for designs more complex than two conditions or with non-linear effects. In contrast, UnfoldSim.jl is free and open-source and it allows to simulate continuous EEG signals even for complex designs.

# Functionality
The package provides four abstract types: `AbstractDesign`, `AbstractComponent`, `AbstractOnset` and `AbstractNoise`. In the following, we present the concrete types that are currently implemented.

## Experimental designs
The design contains the levels of all conditions and predictors. Currently, we support a single and a multi-subject design. The multi-subject design uses the `MixedModelsSim.jl` package [@phillip_alday_2024_10669002] and allows a flexible specification of the random-effects structure by indicating which predictors are within- or between-subject (or item). Designs can be encapsulated, for instance, using the `RepeatDesign` type which repeats the generated event table multiple times, thus generating new trials.

## Event basis functions (Components)
`UnfoldSim.jl` provides a `LinearModelComponent` and a `MixedModelComponent` for single- and multi-subject simulation respectively. These components determine the shape of the response to an event. They consist of a basis function which is weighted by the user-defined regression model. Users can specify a basis function by providing a custom vector or selecting from predefined options, such as simplified EEG components like the N170, modelled as temporally shifted Hanning windows. Further, in the components’ model formulae, fixed effects ($\beta s$) and random effects  (`MultiSubjectDesign`s only) need to be specified.

Each component can be nested in a `MultichannelComponent`, which, using a forward head model, projects the simulated source component to the multi-channel electrode space. Using `Artifacts.jl` we provide on-demand access to the HArtMuT [@harmening2022hartmut] model. 

To generate complex activations, it is possible to specify a vector of `<:AbstractComponents`.

## Inter-onset distributions
The inter-onset distribution defines the distance between events in the continuous (EEG) signal. Currently, `UniformOnset` and `LogNormalOnset` are implemented (see \autoref{fig_onset_distributions}). By adjusting the distribution's parameters, one indirectly controls the amount of overlap between the event-related responses.

![Illustration of the inter-onset distributions. The colours indicate different sets of parameter values. Please note that for the lognormal distribution, the parameters are defined on a logarithmic scale, while the distribution is shown on a linear scale. \label{fig_onset_distributions}](plots/onset_distributions.svg){height="180pt"}

## Noise types
UnfoldSim.jl offers different noise types: `WhiteNoise`, `RedNoise`, `PinkNoise` and exponentially decaying autoregressive noise (`ExponentialNoise`) (see \autoref{fig_noise_types}). In the future, we will add simple autoregressive noise and noise based on actual EEG data.

![Illustration of the different noise types (indicated by colour). Panel **A** shows the noise over time. Please note that the noise signals are shifted by 5&nbsp;µV for visualisation purposes. Panel **B** displays its $\text{log}_{\text{10}}\text{(power)}$ at normalized frequencies. \label{fig_noise_types}](plots/noise_types.svg){height="250pt"}

# Related tools
Few toolboxes for simulating EEG data exist, most being proprietary MATLAB tools that have often not received any updates in the past years or have very specific applications (e.g. `EEGg` [@vaziri2023eegg], `SimMEEG` [@herdman2021simmeeg], `SEED-G` [@anzolin2021seed], `EEGSourceSim` [@barzegaran2019eegsourcesim], `simBCI` [@lindgren2018simbci]). 

In the following, we highlight two actively developed MATLAB-based tools: `Brainstorm` [@tadel2011brainstorm] which especially excels at visualizing the forward model and generating ERPs from phase-aligned oscillations, and `SEREEGA` [@krol2018sereega], which offers comprehensive simulation capabilities with a focus on ERP-component simulation, tools for benchmarking like signal-to-noise specification and more realistic noise simulation (e.g. via random sources).

In Python, `MNE-Python` [@GramfortEtAl2013a] provides some tutorials to simulate EEG data, but the functionality is very basic. `HNN-Core` [@Jas2023] can simulate realistic EEG data by parameterising the neuronal activity in cortical columns.

In contrast to these tools, `UnfoldSim.jl` has a higher-level perspective, uniquely focusing on the regression-ERP aspect. It provides functions to simulate multi-condition experiments, uniquely allows for modelling multi-subject EEG datasets, and offers support to model continuous EEG data with overlapping events. Further, the implementation in Julia offers a platform that is free, actively encourages research software engineering methods, makes it easy to add custom expansions via the `AbstractTypes`, and allows easy access from Python and R.

# Acknowledgements
Funded by Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) under Germany's Excellence Strategy – EXC 2075 – 390740016. The authors further thank the International Max Planck Research School for Intelligent Systems (IMPRS-IS) for supporting Judith Schepers. Moreover, the authors would like to thank Tanja Bien for her valuable feedback on the paper manuscript.

# Package references
Please note that we list only the main dependencies here, but the dependencies of the dependencies can be found in the respective `Manifest.toml` files. Furthermore, please note that we only list rather than cite the packages for which we could not find any citation instructions.

**Julia** [@Julia-2017], **DataFrames.jl** [@JSSv107i04], **Distributions.jl** [@Distributions.jl-2019; @JSSv098i16], **Documenter.jl** [@Documenter_jl], **DSP.jl** [@simon_kornblith_2023_8344531], **FileIO.jl**, **Glob.jl**, **HDF5.jl**, **HypothesisTests.jl**, **ImageFiltering.jl**, **Literate.jl**, **LiveServer.jl**, **Makie.jl** [@DanischKrumbiegel2021], **MixedModels.jl** [@douglas_bates_2023_10268806], **MixedModelsSim.jl** [@phillip_alday_2024_10669002], **Parameters.jl**, **PrettyTables.jl** [@ronan_arraes_jardim_chagas_2023_10214175], **ProjectRoot.jl**, **SignalAnalysis.jl**, **Statistics.jl**, **StableRNGs.jl**, **StatsBase.jl**, **StatsModels.jl**, **TimerOutputs.jl**, **ToeplitzMatrices.jl**, **Unfold.jl** [@ehinger2019unfold], **UnfoldMakie.jl** [@mikheev_2023_10235220]  

# References
