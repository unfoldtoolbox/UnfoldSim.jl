# Installing Julia & UnfoldSim.jl

## Installing Julia

The recommended way to install julia is [juliaup](https://github.com/JuliaLang/juliaup).

TL;DR: If you don't want to read the explicit instructions, just copy the following command:

- Windows: `winget install julia -s msstore`
- Mac/Linux: `curl -fsSL https://install.julialang.org | sh`

We further recommend to use [VSCode](https://code.visualstudio.com/download) and install the [Julia Extension](https://marketplace.visualstudio.com/items?itemName=julialang.language-julia).

## Installing UnfoldSim.jl

The following instructions are intended for Julia beginners. More advanced Julia users can jump ahead to step 3.

#### 1. Start an interactive Julia session ("REPL")
- Option 1: Type `julia` in the command line.
- Option 2: In VSCode, press `Ctrl + Shift + P` to open the command palette and type in `Julia: Start REPL`.

#### 2. Activate your project environment
Before installing UnfoldSim.jl make sure that you activated your project environment.
- Option 1: `cd("/path/to/your/project")` and `]activate .`
- Option 2: `]activate /path/to/your/project/`

!!! hint
    After activating you should see `(environment) pkg>` (where `environment` is the name of your project folder). If you see `(@v1.11) pkg>` instead, you still have to activate your environment.

Note that by typing `]` you enter the Julia package manager. To get back to the Julia REPL, press backspace.

#### 3. Install the UnfoldSim.jl package
- If you are not in the package manager anymore, type `]`.
- Then type `add UnfoldSim`.

After the installation is finished you can use `using UnfoldSim` in the REPL to import the package.