# [![logo_UnfoldSim jl_120px](https://github.com/unfoldtoolbox/UnfoldSim.jl/assets/57703446/139a06c7-55c6-4c2e-8935-627a3c3bf036)](https://github.com/unfoldtoolbox/UnfoldSim.jl/tree/main)

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://unfoldtoolbox.github.io/UnfoldSim.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://unfoldtoolbox.github.io/UnfoldSim.jl/dev/)
[![Build Status](https://github.com/unfoldtoolbox/UnfoldSim.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/unfoldtoolbox/UnfoldSim.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/v/UnfoldSim.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/unfoldtoolbox/UnfoldSim.jl)
[![DOI](https://zenodo.org/badge/413455526.svg)](https://zenodo.org/badge/latestdoi/413455526)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.06641/status.svg)](https://doi.org/10.21105/joss.06641)

|rERP|EEG visualisation|EEG Simulations|BIDS pipeline|Decode EEG data|Statistical testing|
|---|---|---|---|---|---|
| <a href="https://github.com/unfoldtoolbox/Unfold.jl/tree/main"><img src="https://github-production-user-asset-6210df.s3.amazonaws.com/10183650/277623787-757575d0-aeb9-4d94-a5f8-832f13dcd2dd.png"></a> | <a href="https://github.com/unfoldtoolbox/UnfoldMakie.jl"><img  src="https://github-production-user-asset-6210df.s3.amazonaws.com/10183650/277623793-37af35a0-c99c-4374-827b-40fc37de7c2b.png"></a>|<a href="https://github.com/unfoldtoolbox/UnfoldSim.jl"><img src="https://github-production-user-asset-6210df.s3.amazonaws.com/10183650/277623795-328a4ccd-8860-4b13-9fb6-64d3df9e2091.png"></a>|<a href="https://github.com/unfoldtoolbox/UnfoldBIDS.jl"><img src="https://github-production-user-asset-6210df.s3.amazonaws.com/10183650/277622460-2956ca20-9c48-4066-9e50-c5d25c50f0d1.png"></a>|<a href="https://github.com/unfoldtoolbox/UnfoldDecode.jl"><img src="https://github-production-user-asset-6210df.s3.amazonaws.com/10183650/277622487-802002c0-a1f2-4236-9123-562684d39dcf.png"></a>|<a href="https://github.com/unfoldtoolbox/UnfoldStats.jl"><img  src="https://github-production-user-asset-6210df.s3.amazonaws.com/10183650/277623799-4c8f2b5a-ea84-4ee3-82f9-01ef05b4f4c6.png"></a>|

A Julia package to simulate multivariate time series, e.g. model-based ERPs, fMRI activity, pupil dilation etc.
UnfoldSim.jl provides multi-channel support via EEG-forward models. Moreover, it is possible to simulate overlapping event-related activity and to add noise of a certain type e.g. Pink noise.

Many tutorials, guides, how-tos and references are available in the [documentation](https://unfoldtoolbox.github.io/UnfoldSim.jl/)!

![unfoldsim_animation](https://github.com/unfoldtoolbox/UnfoldSim.jl/blob/assets/docs/src/assets/UnfoldSim_features_animation.gif)

## Installation

### Julia
<details>
<summary>Click to expand</summary>

The recommended way to install julia is [juliaup](https://github.com/JuliaLang/juliaup).

TL;DR: If you don't want to read the explicit instructions, just copy the following command:

- Windows: `winget install julia -s msstore`
- Mac/Linux: `curl -fsSL https://install.julialang.org | sh`
</details>

### UnfoldSim.jl

```julia
using Pkg
Pkg.add("UnfoldSim")
```

## Quickstart
We offer some predefined (EEG) signals to get started.

```julia
using UnfoldSim
data, events = UnfoldSim.predef_eeg(; n_repeats = 1, noiselevel = 0.8)
```
Produces continuous "EEG" with PinkNoise and some overlap between 20 events (2 conditions * 10 levels of the continuous variable).

## Slightly longer
All simulation ingredients (design, components, onsets, noise) can be easily modified and you can simply plugin your own!

```julia
using UnfoldSim
using Random

# Start by defining the design / events data frame.
design =
    SingleSubjectDesign(; conditions = Dict(:condA => ["levelA", "levelB"])) |>
    d -> RepeatDesign(d, 10)

# Next define a ground truth signal + relation to events/design with Wilkinson formulas.
signal = LinearModelComponent(;
    basis = [0, 0, 0, 0.5, 1, 1, 0.5, 0, 0],
    formula = @formula(0 ~ 1 + condA),
    β = [1, 0.5],
)
# Finally, define some inter-onset distance distribution and noise, and simulate data!
data, events = simulate(
    Random.MersenneTwister(1),
    design,
    signal,
    UniformOnset(; offset = 5, width = 4),
    PinkNoise(),
)    
```

## Statement of need
EEG researchers often analyze data containing (temporally) overlapping events (e.g. stimulus onset and button press, or consecutive eye-fixations), non-linear effects, and complex experimental designs. For a multitude of reasons, we often need to simulate such kinds of data: Simulated EEG data is useful to test preprocessing and analysis tools, validate statistical methods, illustrate conceptual issues, test toolbox functionalities, and find limitations of traditional analysis workflows. For instance, such simulation tools allow for testing the assumptions of new analysis algorithms and testing their robustness against any violation of these assumptions.

<!--
Note: The statement of need is also used on the documentation landing page (https://unfoldtoolbox.github.io/UnfoldSim.jl/stable/). Make sure that they are synchronized.
-->

## Contributions
Contributions of any kind are very welcome. Please have a look at [CONTRIBUTING.md](https://github.com/unfoldtoolbox/UnfoldSim.jl/blob/main/CONTRIBUTING.md) for guidance on contributing to UnfoldSim.jl.

## Contributors
<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<table>
  <tbody>
    <tr>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/maanikmarathe"><img src="https://avatars.githubusercontent.com/u/66105649?v=4?s=100" width="100px;" alt="Maanik Marathe"/><br /><sub><b>Maanik Marathe</b></sub></a><br /><a href="#doc-maanikmarathe" title="Documentation">📖</a> <a href="#code-maanikmarathe" title="Code">💻</a></td>
      <td align="center" valign="top" width="14.28%"><a href="http://www.benediktehinger.de"><img src="https://avatars.githubusercontent.com/u/10183650?v=4?s=100" width="100px;" alt="Benedikt Ehinger"/><br /><sub><b>Benedikt Ehinger</b></sub></a><br /><a href="#bug-behinger" title="Bug reports">🐛</a> <a href="#code-behinger" title="Code">💻</a> <a href="#doc-behinger" title="Documentation">📖</a> <a href="#ideas-behinger" title="Ideas, Planning, & Feedback">🤔</a> <a href="#infra-behinger" title="Infrastructure (Hosting, Build-Tools, etc)">🚇</a> <a href="#maintenance-behinger" title="Maintenance">🚧</a> <a href="#review-behinger" title="Reviewed Pull Requests">👀</a> <a href="#test-behinger" title="Tests">⚠️</a> <a href="#tutorial-behinger" title="Tutorials">✅</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/llips"><img src="https://avatars.githubusercontent.com/u/38983684?v=4?s=100" width="100px;" alt="Luis"/><br /><sub><b>Luis</b></sub></a><br /><a href="#bug-llips" title="Bug reports">🐛</a> <a href="#code-llips" title="Code">💻</a> <a href="#doc-llips" title="Documentation">📖</a> <a href="#ideas-llips" title="Ideas, Planning, & Feedback">🤔</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/jschepers"><img src="https://avatars.githubusercontent.com/u/22366977?v=4?s=100" width="100px;" alt="Judith Schepers"/><br /><sub><b>Judith Schepers</b></sub></a><br /><a href="#ideas-jschepers" title="Ideas, Planning, & Feedback">🤔</a> <a href="#bug-jschepers" title="Bug reports">🐛</a> <a href="#doc-jschepers" title="Documentation">📖</a> <a href="#tutorial-jschepers" title="Tutorials">✅</a> <a href="#code-jschepers" title="Code">💻</a> <a href="#test-jschepers" title="Tests">⚠️</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/vladdez"><img src="https://avatars.githubusercontent.com/u/33777074?v=4?s=100" width="100px;" alt="Vladimir Mikheev"/><br /><sub><b>Vladimir Mikheev</b></sub></a><br /><a href="#bug-vladdez" title="Bug reports">🐛</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://reboreexplore.github.io/"><img src="https://avatars.githubusercontent.com/u/43548330?v=4?s=100" width="100px;" alt="Manpa Barman"/><br /><sub><b>Manpa Barman</b></sub></a><br /><a href="#infra-ReboreExplore" title="Infrastructure (Hosting, Build-Tools, etc)">🚇</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://reneskukies.de/"><img src="https://avatars.githubusercontent.com/u/57703446?v=4?s=100" width="100px;" alt="René Skukies"/><br /><sub><b>René Skukies</b></sub></a><br /><a href="#doc-ReneSkukies" title="Documentation">📖</a> <a href="#code-ReneSkukies" title="Code">💻</a> <a href="#test-ReneSkukies" title="Tests">⚠️</a> <a href="#ideas-ReneSkukies" title="Ideas, Planning, & Feedback">🤔</a> <a href="#tutorial-ReneSkukies" title="Tutorials">✅</a></td>
    </tr>
  </tbody>
</table>

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->



This project follows the [all-contributors](https://allcontributors.org/docs/en/specification) specification. 
Please reach out, if you have contributed to UnfoldSim.jl but we have not listed you as a contributor yet.

## Citation

If you use UnfoldSim.jl, please acknowledge and support our work by citing
1\. Our [JOSS paper](https://doi.org/10.21105/joss.06641): 

```
Schepers et al., (2025). UnfoldSim.jl: Simulating continuous event-based time series data for EEG and beyond. Journal of Open Source Software, 10(107), 6641, https://doi.org/10.21105/joss.06641
```

<details>
  <summary>BibTeX entry:</summary>

```bib
@article{Schepers2025,
 doi = {10.21105/joss.06641},
 url = {https://doi.org/10.21105/joss.06641},
 year = {2025},
 publisher = {The Open Journal},
 volume = {10},
 number = {107},
 pages = {6641},
 author = {Judith Schepers and Luis Lips and Maanik Marathe and Benedikt V. Ehinger},
 title = {UnfoldSim.jl: Simulating continuous event-based time series data for EEG and beyond},
 journal = {Journal of Open Source Software}
 } 
```
</details>

and

2\. The corresponding [Zenodo DOI](https://doi.org/10.5281/zenodo.7738651) for the specific UnfoldSim.jl version that you are using in your work.

## Acknowledgements

Funded by Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) under Germany´s Excellence Strategy – EXC 2075 – 390740016. Furthermore, the authors thank the International Max Planck Research School for Intelligent Systems (IMPRS-IS) for supporting Judith Schepers.
