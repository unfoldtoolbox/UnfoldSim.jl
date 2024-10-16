# [![logo_UnfoldSim jl_120px](https://github.com/unfoldtoolbox/UnfoldSim.jl/assets/57703446/139a06c7-55c6-4c2e-8935-627a3c3bf036)](https://github.com/unfoldtoolbox/UnfoldSim.jl/tree/main)

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://unfoldtoolbox.github.io/UnfoldSim.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://unfoldtoolbox.github.io/UnfoldSim.jl/dev/)
[![Build Status](https://github.com/unfoldtoolbox/UnfoldSim.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/unfoldtoolbox/UnfoldSim.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/v/UnfoldSim.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/unfoldtoolbox/UnfoldSim.jl)
[![DOI](https://zenodo.org/badge/413455526.svg)](https://zenodo.org/badge/latestdoi/413455526)

|rERP|EEG visualisation|EEG Simulations|BIDS pipeline|Decode EEG data|Statistical testing|
|---|---|---|---|---|---|
| <a href="https://github.com/unfoldtoolbox/Unfold.jl/tree/main"><img src="https://github-production-user-asset-6210df.s3.amazonaws.com/10183650/277623787-757575d0-aeb9-4d94-a5f8-832f13dcd2dd.png"></a> | <a href="https://github.com/unfoldtoolbox/UnfoldMakie.jl"><img  src="https://github-production-user-asset-6210df.s3.amazonaws.com/10183650/277623793-37af35a0-c99c-4374-827b-40fc37de7c2b.png"></a>|<a href="https://github.com/unfoldtoolbox/UnfoldSim.jl"><img src="https://github-production-user-asset-6210df.s3.amazonaws.com/10183650/277623795-328a4ccd-8860-4b13-9fb6-64d3df9e2091.png"></a>|<a href="https://github.com/unfoldtoolbox/UnfoldBIDS.jl"><img src="https://github-production-user-asset-6210df.s3.amazonaws.com/10183650/277622460-2956ca20-9c48-4066-9e50-c5d25c50f0d1.png"></a>|<a href="https://github.com/unfoldtoolbox/UnfoldDecode.jl"><img src="https://github-production-user-asset-6210df.s3.amazonaws.com/10183650/277622487-802002c0-a1f2-4236-9123-562684d39dcf.png"></a>|<a href="https://github.com/unfoldtoolbox/UnfoldStats.jl"><img  src="https://github-production-user-asset-6210df.s3.amazonaws.com/10183650/277623799-4c8f2b5a-ea84-4ee3-82f9-01ef05b4f4c6.png"></a>|

A Julia package to simulate multivariate time series, e.g. model-based ERPs, fMRI activity, pupil dilation etc.
UnfoldSim.jl provides multi-channel support via EEG-forward models. Moreover, it is possible to simulate overlapping event-related activity and to add noise of a certain type e.g. Pink noise.

Many tutorials, guides, how-tos and references are available in the [documentation](https://unfoldtoolbox.github.io/UnfoldSim.jl/dev/)!

![readme_figure](https://github.com/unfoldtoolbox/UnfoldSim.jl/assets/22366977/b69d186c-fd3d-4449-9f2e-03d7e01b8cb3)

## Install

### Julia
<details>
<summary>Click to expand</summary>

The recommended way to install julia is [juliaup](https://github.com/JuliaLang/juliaup).
It allows you to, e.g., easily update Julia at a later point, but also test out alpha/beta versions etc.

TL:DR; If you dont want to read the explicit instructions, just copy the following command

#### Windows

AppStore -> JuliaUp,  or `winget install julia -s msstore` in CMD

#### Mac & Linux

`curl -fsSL https://install.julialang.org | sh` in any shell
</details>

### Unfold.jl

```julia
using Pkg
Pkg.add("UnfoldSim")
```

# Quickstart
```julia
using UnfoldSim
data, events = UnfoldSim.predef_eeg(; n_repeats = 1, noiselevel = 0.8)
```
Produces continuous "EEG" with PinkNoise and some overlap between 20 events (2 conditions * 10 levels of the continuous variable).

## Slightly longer
```julia
using UnfoldSim
using Random

# Start by defining the design / events data frame
design =
    SingleSubjectDesign(; conditions = Dict(:condA => ["levelA", "levelB"])) |>
    d -> RepeatDesign(d, 10);

# Next define a ground truth signal + relation to events/design with Wilkinson formulas
signal = LinearModelComponent(;
    basis = [0, 0, 0, 0.5, 1, 1, 0.5, 0, 0],
    formula = @formula(0 ~ 1 + condA),
    β = [1, 0.5],
);
# finally, define some inter-onset distribution and noise, and simulate data!
data, events = simulate(
    Random.MersenneTwister(1),
    design,
    signal,
    UniformOnset(; offset = 5, width = 4),
    PinkNoise(),
);    
```
All simulation ingredients (design, components, onsets, noise) can be easily modified and you can simply plugin your own!

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
      <td align="center" valign="top" width="14.28%"><a href="https://reneskukies.de/"><img src="https://avatars.githubusercontent.com/u/57703446?v=4?s=100" width="100px;" alt="René Skukies"/><br /><sub><b>René Skukies</b></sub></a><br /><a href="#doc-ReneSkukies" title="Documentation">📖</a></td>
    </tr>
  </tbody>
</table>

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->



This project follows the [all-contributors](https://allcontributors.org/docs/en/specification) specification. 
Please reach out, if you have contributed to UnfoldSim.jl but we have not listed you as a contributor yet.

## Citation

TBA

## Acknowledgements

Funded by Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) under Germany´s Excellence Strategy – EXC 2075 – 390740016. Furthermore, the authors thank the International Max Planck Research School for Intelligent Systems (IMPRS-IS) for supporting Judith Schepers.
