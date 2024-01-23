
---
title: 'Unfold.jl: Regression models for Cognitive Neuroscience'
tags:
  - Julia
  - EEG
  - neuroimaging
authors:
  - name: Benedikt V. Ehinger
    orcid:  0000-0002-6276-3332 
    equal-contrib: false
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
  
affiliations:
 - name: Stuttgart Center for Simulation Science, University of Stuttgart, Germany
   index: 1
 - name: Institute for Visualization and Interactive Systems, University of Stuttgart, Germany
   index: 2
date: 06 November 2023
bibliography: paper.bib


---
# Summary

UnfoldSim.jl is used to simulate multivariate timeseries, with a focus on EEG, especially event-related potentials (ERPs). The user provides four ingredients: 1) an experimental design, with both categorical and continuous variables, 2) event basis functions specified via linear or hierarchical models, 3) an inter-event onset distribution, and 4) a noise specification. Unfold.jl simulates continuous EEG signals with potentially overlapping events. Multi-channel support via EEG-forward models is available as well. UnfoldSim.jl is modular, providing intuitive entrance points for custom requirements. For instance, support for other modalities, e.g. single-voxel fMRI or pupil dilation signals is easily provided.

# Statement of Need
Simulated EEG data is necessary to test preprocessing and analysis tools, to illustrate issues and to test functions. While other simulation tools exist, they are dominantly based in Matlab (below) and they do not adress our unique challenges and requirements. In our work (e.g. Ehinger & Dimigen 2019), we focus in regression-based deconvolution of ERPs (Smith & Kutas 2015 ab). In short, multiple-regression is used for linear overlap correction, non-linear (hierarchical) effects fitting (similar to GAMMs, Wood et al), often applied to use-cases where events overlap in time (e.g. stimulus and button press or consecutive eye-fixations). The tools used up to now (e.g. mTRF-toolbox, Unfold, Unfold.jl, fitgrid, - add citations!)

# Functionality
The toolbox provides four abstract components: AbstractDesign, AbstractComponent, AbstractOnset and AbstractNoise.

## Concrete Designs
Currently we support a single, and a multi-subject design. They are used to generate an experimental design containing the conditions and levels of all predictors. Randomisation is possible via a user-defined function which is applied after design-generation. Designs could be nested, like our RepeatDesign which simply repeats the generated design multiple times.

## Concrete Components
We provide a LinearModelComponent and a MixedModelComponent for multi-subject simulation respectively. For the components model-formulae, fixed-effects ($\beta s$) and random effects need to be specified. Further the coding-schema can be provided following StatsModels.jl. The component MultichannelComponent can be applied to any component and allows for projecting the simulated source-component to multi-channel electrodes via a headmodel. Using Artifacts.jl we provide on-demand access to the Hartmut (cite) model.

## Concrete Noise
We provide different noise-types "White","Red" and "Pink", but also an exponentially declining Autoregressive noise type.

## Concerete Onsets
The onsets define the distance between events for the continuous EEG. Currently UniformOnset and LogNormalOnset are implemented.



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
