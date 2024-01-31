---
title: 'UnfoldSim.jl: A toolbox for simulating continuous event-based time series data for EEG and beyond'

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
date: 31 January 2024
bibliography: paper.bib
---

# Summary

UnfoldSim.jl is a Julia package used to simulate multivariate time series, with a focus on EEG, especially event-related potentials (ERPs). The user provides four ingredients: 1) an experimental design, with both categorical and continuous variables, 2) event basis functions specified via linear or hierarchical models, 3) an inter-event onset distribution, and 4) a noise specification. UnfoldSim.jl then simulates continuous EEG signals with potentially overlapping events. Multi-channel support via EEG-forward models is available as well. UnfoldSim.jl is modular, providing intuitive entrance points for individual customizations. The user can implement custom designs, components, onset distributions or noise types to tailor the toolbox to their needs. This allows support even for other modalities, e.g. single-voxel fMRI or pupil dilation signals.

# Statement of Need
In our work (e.g. @ehinger2019unfold, @dimigen2021regression), we often analyze data containing (temporally) overlapping events (e.g. stimulus onset and button press, or consecutive eye-fixations), non-linear effects, and complex experimental designs. For a multitude of reasons, we need to simulate such kind of data: Simulated EEG data is necessary to test preprocessing and analysis tools, validate statistical methods, illustrate conceptual issues, test toolbox functionalities, and find limitations of traditional analysis workflows.

While other EEG simulation toolboxes exist, they each have limitations: they are dominantly MATLAB-based, they do not simulate continuous EEG, and they offer little support for designs more complex than two conditions or with non-linear effects.

# Functionality
The toolbox provides four abstract types: `AbstractDesign`, `AbstractComponent`, `AbstractOnset` and `AbstractNoise`.

## Concrete Designs
Currently, we support a single and a multi-subject design. They are used to generate an experimental design containing the conditions and levels of all predictors. The multi-subject design uses the MixedModelsSim.jl toolbox [@phillip_alday_2022_7407741] and allows a flexible specification of the random-effects structure by indicating which predictors are within- or between-subject (or item). Tailored randomisation is possible via a user-specified function, which is applied after design generation. Designs can be encapsulated, for instance, the RepeatDesign-type which repeats the generated event tables multiple times, thus generating new trials. Currently, only balanced designs are implemented, i.e. all possible combinations of predictor levels have the same number of trials. However, [a tutorial on how to implement a new design](https://unfoldtoolbox.github.io/UnfoldSim.jl/dev/generated/HowTo/newDesign) for imbalanced datasets is provided.

## Concrete Components
UnfoldSim.jl provides a LinearModelComponent and a MixedModelComponent for multi-subject simulation respectively. These components determine the shape of the response to an event. They consist of a basis function which is multiplied with the user-defined coefficient of a regression model. The user specifies a basis function for the component by either providing a custom vector or choosing one of the prespecified bases. For example, the toolbox provides simplified versions of typical EEG components e.g. N170 which are implemented as temporally shifted hanning windows. Further, in the components’ model formulae, fixed-effects ($\beta s$) and random effects  (MultiSubject designs only) need to be specified.

Each component can be nested in a MultichannelComponent, which, using a forward headmodel, projects the simulated source component to the multi-channel electrode space. Using Artifacts.jl we provide on-demand access to the Hartmut [@harmening2022hartmut] model. 

To generate complex activations, it is possible to specify a vector of `<:AbstractComponents`.

## Concerete Onsets
The inter-onset distribution defines the distance between events in the case of a continuous EEG. Currently, UniformOnset and LogNormalOnset are implemented. By specifying the parameters of the onset distribution, one indirectly controls the amount of overlap between two or more event-related responses.
\ref{fig:onset_distributions} illustrates the parameterization of the two implemented onset distributions. \autoref{fig:onset_distributions2}
![Caption for example figure.\label{fig:onset_distributions}](plots/onset_distributions.pdf){width=90%}
![Caption for example figure.]{label="fig:onset_distributions2"}(plots/onset_distributions.pdf){width=90%}

## Concrete Noise
We provide different noise-types "White","Red" and "Pink", but also an exponentially declining Autoregressive noise type.

# Related tools
Not many toolboxes for simulating EEG data exist. Nearly all toolboxes have been developed in MATLAB (e.g. EEGg - Vaziri, SimMEEG - Herdman 2021, SEED-G-Toolbox Anzolin 2021) but are largely abandoned. MNE-Python provides basic tutorials to simulate EEG data as well, but no dedicated functionality. Thus we highlight in the following two other excellent matlab-based tools: Brainstorm (Tadel 2021) and SEREEGA (Krol 2017). Both toolboxes are based in matlab and provide forward-simulation of EEG signals. Brainstorm especially excells at the visualization of the forward-model, and provides interesting capabilities to generate ERPs based on phase-aligned oscillations. SEREEGA provides the most complete simulation capabilities with a greater focus on ERP-component simulation, tools for benchmarking like signal-to-noise specification, and more realistic noise simulation (e.g. via random sources).

In relation to these tools, Unfold uniquely focuses on the regression-ERP aspect, providing functions to simulate multi-condition experiments, uniquely allows to model hierarchical, that is, multi-subject EEG datasets, and offers unique support to model continuous EEG data with overlap events.

Due to its different focus, UnfoldSim.jl currently lacks advanced visualizations of leadfields and does not provide any tools for simulating oscillations or phase-based simulation of ERPs.



# Other notes
SEEREGA - Matlab, best in class, no continuous data
BRAINSTORM - Matlab, good for forward-model simulation
EEGg - Vaziri - Matlab, very limited functionality
SimMEEG - Herdman 2021 -  Matlab, Discontinued
SEED-G-toolbox - Anzolin 2021 - Matlab, Discontinued, 

MNE-python - matlab, only basic tutorials

## JOSS Stuff
Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

Test
@Documenter_jl

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.

# References
