# Developer documentation

Welcome to the developer documentation! We are excited that you are interested in contributing to our package.

Feel free to submit your work in a state you are comfortable with—we genuinely appreciate every contribution! If you are interested in following best practices and learning along the way, keep reading. But don't worry, we welcome your input just as it is 🙂.


!!! note "Contribution guide"
    If you haven't already, please read the [Contribution guide](https://github.com/unfoldtoolbox/UnfoldSim.jl/blob/main/CONTRIBUTING.md) first.

Please note that the following documentation is adapted from the [BestieTemplate.jl developer documentation](https://juliabesties.github.io/BestieTemplate.jl/stable/91-developer/) but has been customized to fit our needs.

## Development/GitHub Workflow

### Before you start coding
1. Check whether there exists a GitHub issue about the topic (e.g. bug, new feature etc). If not create one with a short description of the problem or feature idea.
2. Discuss your approach with the package maintainers either in the issue or via another channel.

### First time clone & development version

If this is the first time you work with this repository, follow the instructions below to clone the repository and create a `dev` version. 

!!! tip "`dev` version of a Julia package"
    Having a `dev` (development) version of a package allows you to import a local version of the package with your changes instead of the registered package version (which is static).

#### a) If you have writing access for the GitHub repository
- Option 1: Clone this repository using `git clone`.
- Option 2 (recommended): Use the Julia `dev` command to create a development version of the package:
    1. Start a Julia session and run `cd("/path/to/your/project")` to navigate to your project folder.
    2. Press `]` to enter `pkg` mode.
    3. Run `dev --local UnfoldSim` to clone the package to `./dev/UnfoldSim` and automatically add it to your Julia project environment.

!!! important
    If you have writing rights, whenever **upstream** is mentioned, use **origin** instead.

#### b) If you don't have writing access for the GitHub repository
1. Fork this repository.
2. Clone your repository (this will create a `git remote` called `origin`).
3. Add this repository as a remote:

   ```bash
   git remote add upstream https://github.com/unfoldtoolbox/UnfoldSim.jl
   ```

This will ensure that you have two remotes in your git: `origin` and `upstream`.
You will create branches and push to `origin`, and you will fetch and update your local `main` branch from `upstream`.

!!! tip "`dev` version without writing rights"
    You can also use the `dev` command on your fork. Run `]dev --local url/of/your/fork` to clone the package to `./dev/UnfoldSim` and automatically add it to your Julia project environment.

### Revise.jl
Further, we recommend to use [`Revise.jl`](https://github.com/timholy/Revise.jl): a Julia package which allows you to track source code changes in a running Julia session without need to restart it and reload the package.

We recommend to install it in the global environment:
```julia-repl
julia> # Press ]
pkg> activate
pkg> add Revise
```

### Working on a new issue

We try to keep a linear history in this repo, so it is important to keep your branches up-to-date.

1. Fetch from the remote and fast-forward your local main

   ```bash
   git fetch upstream
   git switch main
   git merge --ff-only upstream/main
   ```

2. Branch from `main` to address the issue (see below for naming)

   ```bash
   git switch -c 42-add-answer-universe
   ```

3. Push the new local branch to your personal remote repository

   ```bash
   git push -u origin 42-add-answer-universe
   ```

4. Create a pull request to merge your remote branch into the org main.

#### Branch naming

- If there is an associated issue, add the issue number.
- If there is no associated issue, **and the changes are small**, add a prefix such as "typo", "hotfix", "small-refactor", according to the type of update.
If the changes are not small and there is no associated issue, then either create an issue first, or discuss in another channel with the maintainers.
- Use dash separated imperative wording related to the issue (e.g., `14-add-tests`, `15-fix-model`, `16-remove-obsolete-files`).

#### Commit messages

- Use imperative or present tense, for instance: *Add feature* or *Fix bug*.
- Have informative titles.
- When necessary, add a body with details.
- If there are breaking changes, add the information to the commit message.

### Before creating a pull request

- Ideally: Make sure the tests pass.
- Add appropriate documentation (ideally using the [Docstring templates](@ref)).
- Ideally: Follow the formatting rules from `JuliaFormatter.jl` (see [Formatting](@ref)).
- Fetch any `main` updates from upstream and rebase your branch, if necessary:

  ```bash
  git fetch upstream
  git rebase upstream/main BRANCH_NAME
  ```

- Then you can open a pull request and work with the reviewer to address any issues.

!!! important "Best practices not shackles"
    We encourage you to share your contributions in whatever state you are comfortable with. Don’t feel overwhelmed by the number of guides — think of them as helpful resources, not strict requirements. Every contribution is valuable, and we’re happy to refine things together! 

## Formatting

### JuliaFormatter.jl
We use [JuliaFormatter.jl](https://github.com/domluna/JuliaFormatter.jl) for formatting and recommend you to install it in your global environment:
```julia-repl
julia> # Press ]
pkg> activate
pkg> add JuliaFormatter
```

### Beware of reviewdog 🐶
We use the [julia-format](https://github.com/julia-actions/julia-format) Github action to ensure that the code follows the formatting rules defined by [JuliaFormatter.jl](https://github.com/domluna/JuliaFormatter.jl).
When opening a pull request [reviewdog](https://github.com/reviewdog/reviewdog) will automatically make formatting suggestions for your code.

## Testing

As with most Julia packages, you can just open Julia in the repository folder, activate the environment, and run `test`:

```julia-repl
julia> # press ]
pkg> activate .
pkg> test
```

!!! tip "Running single tests"
    Instead fo running all tests, you can also run the `test/setup.jl` to load all required packages, and subsequently run single tests manually either by `include("test/my_test.jl")` or by opening the file and running the specific test block you want to run.

## Documentation

Documentation is key to maintaining a codebase that is easy to understand and extend. Whether it's comments in the code, docstrings, or tutorials, when writing documentation, think about your future self or the next person reading the code or using your functions.

### Building and viewing the documentation locally

We recommend using [LiveServer](https://github.com/JuliaDocs/LiveServer.jl) to build and preview the documentation locally.
To simplify this process we created the `docs/run_liveserver.jl` script.

Please follow these steps:

1. Navigate to the `docs` folder and activate it.
2. Run `using Revise` (in case you decided to install it in your docs environment).
3. If this is the first time building the docs
   1. Press `]` to enter `pkg` mode.
   2. Run `pkg> dev ..` to use the development version of your package.
   3. Press backspace to leave `pkg` mode.
4. Run `include("run_liveserver.jl")`.
5. Click on the provided link or go to `http://0.0.0.0:8000/` in your browser.

!!! tip "Live preview of docstrings"
    Install `Revise.jl` in the docs environment to enable live updating of docstrings in the docs preview.

!!! tip "Separate Julia session for docs preview"
    We recommend using a separate Julia session (in VSCode) to run the `run_liveserver.jl` script, as it continues running. This way, you can avoid "blocking" the REPL and run other code in the meantime.

### Adding a documentation page
1. We recommend to write a Literate.jl document and place it in `docs/literate/FOLDER/FILENAME.jl` with `FOLDER` being `HowTo`, `Explanation`, `Tutorial` or `Reference` ([recommended reading on the 4 categories](https://documentation.divio.com/)).
2. Literate.jl converts the `.jl` file to a `.md` automatically and places it in `docs/src/generated/FOLDER/FILENAME.md`.
3. Edit [make.jl](https://github.com/unfoldtoolbox/Unfold.jl/blob/main/docs/make.jl) with a reference to `docs/src/generated/FOLDER/FILENAME.md`.

### Docstring templates
The following docstring templates are mainly based on the [Julia manual](https://docs.julialang.org/en/v1/manual/documentation/) and [Blue: a Style Guide for Julia](https://github.com/JuliaDiff/BlueStyle?tab=readme-ov-file#documentation).

#### Function docstring template
````julia
"""
    my_function(argument1:Type1; keyword_argument3::Type3 = value3)
    my_function(argument1::Type1, optional_argument2::Type2; keyword_argument3::Type3 = value3)

One-line description using the imperative form ("Do this") instead of the third person and ending with a period.

If the one-line description is not sufficient, one can also write a short paragraph with additional information.

# Arguments (if needed)
- `argument1::Type1`: Description of argument1.
- `optional_argument2::Type2` (optional): Description of optional_argument2.

# Keyword arguments (if needed)
- `keyword_argument3::Type3 = value3`: Description of keyword_argument3.

# Returns
- `result::Type4` : Description of result.

# Examples
```julia-repl
julia> my_function(value1, value2)
result1

julia> my_function(value1; keyword_argument3 = value4)
result2
```

See also [`my_function2`](@ref), [`my_function3`](@ref).
"""
````
**Special cases:**
- If a function accepts many keyword arguments, only include `<keyword arguments>` as a placeholder in the signature and give a keyword list with descriptions in the Keyword arguments section of the docstring.
- If a function returns more than one variable, write the `Returns` section in the following way:
  ```julia
  # Returns
  - (result1, result2)::Tuple{Type1, Type2}:
      - Description of result1
      - Description of result2
  ```

#### Type docstring template
````julia
"""
    MyType <: MyAbstractType

One-line desciption of my type which ends with a period.

If the one-line description is not sufficient, one can also write a short paragraph with additional information.

# Fields
- `field1::Type1`: Description of field1.
- `optional_field2::Type2 = value2` (optional): Description of field2. If not provided, defaults to `value2`.

# Examples
```julia-repl
julia> MyType(field1, field2)
result1
```

See also [`MyType2`](@ref), [`my_function2`](@ref).
"""
````