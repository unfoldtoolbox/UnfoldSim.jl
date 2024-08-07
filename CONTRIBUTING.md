# Contribution guide
Contributions are very welcome. These could be typos, bug reports, feature requests, speed optimization, better code, and better documentation.
You are very welcome to raise issues and start pull requests.

## Bug reports
If you notice any bugs, such as crashing code, incorrect results or speed issues, please raise a GitHub issue. The issue should contain a short description of the problem, (optimally) a minimal working example to reproduce the bug and which UnfoldSim.jl version you are using.

## Code contributions (Pull requests)
When opening a pull request, please add a short but meaningful description of the changes/features you implemented. Moreover, please add tests (where appropriate) to ensure that your code is working as expected.

## Adding documentation
1. We recommend to write a Literate.jl document and place it in `docs/literate/FOLDER/FILENAME.jl` with `FOLDER` being `HowTo`, `Explanation`, `Tutorial` or `Reference` ([recommended reading on the 4 categories](https://documentation.divio.com/)).
2. Literate.jl converts the `.jl` file to a `.md` automatically and places it in `docs/src/generated/FOLDER/FILENAME.md`.
3. Edit [make.jl](https://github.com/unfoldtoolbox/Unfold.jl/blob/main/docs/make.jl) with a reference to `docs/src/generated/FOLDER/FILENAME.md`.

## Formatting (Beware of reviewdog :dog:)
We use the [julia-format](https://github.com/julia-actions/julia-format) Github action to ensure that the code follows the formatting rules defined by [JuliaFormatter.jl](https://github.com/domluna/JuliaFormatter.jl).
When opening a pull request [reviewdog](https://github.com/reviewdog/reviewdog) will automatically make formatting suggestions for your code.
