var documenterSearchIndex = {"docs":
[{"location":"api/","page":"DocStrings","title":"DocStrings","text":"Modules = [UnfoldSim]","category":"page"},{"location":"api/#UnfoldSim.LinearModelComponent","page":"DocStrings","title":"UnfoldSim.LinearModelComponent","text":"A multiple regression component for one subject\n\nbasis: an object, if accessed, provides a 'basis-function', e.g. hanning(40), this defines the response at a single event. It will be weighted by the model-prediction\nformula: StatsModels Formula-Object  @formula 0~1+cond (left side must be 0)\nβ Vector of betas, must fit the formula\ncontrasts: Dict. Default is empty, e.g. Dict(:condA=>EffectsCoding())\n\nAll arguments can be named, in that case contrasts is optional\n\nWorks best with SingleSubjectDesign\n\nLinearModelComponent(;\n    basis=hanning(40),\n    formula=@formula(0~1+cond),\n    β = [1.,2.],\n    contrasts=Dict(:cond=>EffectsCoding())\n)\n\n\n\n\n\n\n","category":"type"},{"location":"api/#UnfoldSim.MixedModelComponent","page":"DocStrings","title":"UnfoldSim.MixedModelComponent","text":"A component that adds a hierarchical relation between parameters according to a LMM defined via MixedModels.jl\n\nbasis: an object, if accessed, provides a 'basis-function', e.g. hanning(40), this defines the response at a single event. It will be weighted by the model-prediction\nformula: Formula-Object in the style of MixedModels.jl e.g. @formula dv~1+cond + (1|subject) - left side must be dv\nβ Vector of betas, must fit the formula\nσs Dict of random effect variances, e.g. Dict(:subject=>[0.5,0.4]) or to specify correlationmatrix Dict(:subject=>[0.5,0.4,I(2,2)],...). Technically, this will be passed to MixedModels.jl create_re function, which creates the θ matrices.\ncontrasts: Dict in the style of MixedModels.jl. Default is empty.\n\nAll arguments can be named, in that case contrasts is optional\n\nWorks best with MultiSubjectDesign\n\nMixedModelComponent(;\n    basis=hanning(40),\n    formula=@formula(dv~1+cond+(1+cond|subject)),\n    β = [1.,2.],\n    σs= Dict(:subject=>[0.5,0.4]),\n    contrasts=Dict(:cond=>EffectsCoding())\n)\n\n\n\n\n\n\n","category":"type"},{"location":"api/#UnfoldSim.MultiSubjectDesign","page":"DocStrings","title":"UnfoldSim.MultiSubjectDesign","text":"n_subjects::Int -> number of subjects\nn_items::Int -> number of items (sometimes ≈trials)\nsubjects_between = nothing -> effects between subjects, e.g. young vs old \nitems_between = nothing -> effects between items, e.g. natural vs artificial images, but shown to all subjects\nboth_within = nothing\t-> effects completly crossed\ntableModifyFun = x->x; # can be used to sort, or x->shuffle(MersenneTwister(42),x) - be sure to fix/update the rng accordingly!!\n\ntipp: check the resulting dataframe using generate(design)\n\n# declaring same condition both sub-between and item-between results in a full between subject/item design\ndesign = MultiSubjectDesignjectDesign(;\n        n_items=10,\n\t\tn_subjects = 30,\n        subjects_between=Dict(:cond=>[\"levelA\",\"levelB\"]),\n\t\titems_between =Dict(:cond=>[\"levelA\",\"levelB\"]),\n        );\n\n\n\n\n\n","category":"type"},{"location":"api/#UnfoldSim.RepeatDesign","page":"DocStrings","title":"UnfoldSim.RepeatDesign","text":"repeat a design DataFrame multiple times to mimick repeatedly recorded trials\n\ndesignOnce = MultiSubjectDesign(;\n        n_items=2,\n\t\tn_subjects = 2,\n        subjects_between =Dict(:cond=>[\"levelA\",\"levelB\"]),\n\t\titems_between =Dict(:cond=>[\"levelA\",\"levelB\"]),\n        );\n\ndesign = RepeatDesign(designOnce,4);\n\n\n\n\n\n","category":"type"},{"location":"api/#UnfoldSim.SingleSubjectDesign","page":"DocStrings","title":"UnfoldSim.SingleSubjectDesign","text":"conditions = Dict of conditions, e.g. Dict(:A=>[\"a_small\",\"a_big\"],:B=>[\"b_tiny\",\"b_large\"])\ntableModifyFun = x->x; # can be used to sort, or x->shuffle(MersenneTwister(42),x) - be sure to fix/update the rng accordingly!!\n\nNumber of trials / rows in generate(design) depend on the full factorial of your conditions.\n\nTo increase the number of repetitions simply use RepeatDesign(SingleSubjectDesign(...),5)\n\ntipp: check the resulting dataframe using generate(design)\n\n\n\n\n\n","category":"type"},{"location":"api/#Base.size-Tuple{MultiSubjectDesign}","page":"DocStrings","title":"Base.size","text":"Returns dimension of experiment design\n\n\n\n\n\n","category":"method"},{"location":"api/#DSP.Windows.hanning-Tuple{Any, Any, Any}","page":"DocStrings","title":"DSP.Windows.hanning","text":"generate a hanning window\n\nduration: in s offset: in s, defines hanning peak sfreq: sampling rate in Hz\n\n\n\n\n\n","category":"method"},{"location":"api/#UnfoldSim.convert-Tuple{Any, Any, Any}","page":"DocStrings","title":"UnfoldSim.convert","text":"Function to convert output similar to unfold (data, evts)\n\n\n\n\n\n","category":"method"},{"location":"api/#UnfoldSim.gen_noise-Tuple{Any, UnfoldSim.RealisticNoise, Int64}","page":"DocStrings","title":"UnfoldSim.gen_noise","text":"gen_noise(t::RealisticNoise, n::Int)\n\nGenerate noise of a given type t and length n\n\n\n\n\n\n","category":"method"},{"location":"api/#UnfoldSim.gen_noise-Tuple{Any, Union{PinkNoise, RedNoise}, Int64}","page":"DocStrings","title":"UnfoldSim.gen_noise","text":"gen_noise(t::Union{PinkNoise, RedNoise}, n::Int)\n\nGenerate noise of a given type t and length n\n\n\n\n\n\n","category":"method"},{"location":"api/#UnfoldSim.gen_noise-Tuple{Any, WhiteNoise, Int64}","page":"DocStrings","title":"UnfoldSim.gen_noise","text":"gen_noise(t::WhiteNoise, n::Int)\n\nGenerate noise of a given type t and length n\n\n\n\n\n\n","category":"method"},{"location":"api/#UnfoldSim.generate-Tuple{MultiSubjectDesign}","page":"DocStrings","title":"UnfoldSim.generate","text":"Generates full factorial Dataframe according to MixedModelsSim.jl 's simdatcrossed function Note: nitems = you can think of it as trials or better, as stimuli\n\nNote: No condition can be named dv which is used internally in MixedModelsSim / MixedModels as a dummy left-side\n\nAfterwards applies expdesign.tableModifyFun.  Could be used to duplicate trials, sort, subselect etc.\n\nFinally it sorts by :subject\n\njulia> d = MultiSubjectDesign(;nsubjects = 10,nitems=20,both_within= Dict(:A=>nlevels(5),:B=>nlevels(2))) julia> generate(d)\n\n\n\n\n\n","category":"method"},{"location":"api/#UnfoldSim.generate-Tuple{SingleSubjectDesign}","page":"DocStrings","title":"UnfoldSim.generate","text":"Generates full-factorial DataFrame of expdesign.conditions\n\nAfterwards applies expdesign.tableModifyFun.\n\njulia> d = SingleSubjectDesign(;conditions= Dict(:A=>nlevels(5),:B=>nlevels(2))) julia> generate(d)\n\n\n\n\n\n","category":"method"},{"location":"api/#UnfoldSim.hrf-Tuple{}","page":"DocStrings","title":"UnfoldSim.hrf","text":"Generate a HRF kernel. \n\nTR = 1/sfreq default parameters taken from SPM\n\nCode adapted from Unfold.jl\n\n\n\n\n\n","category":"method"},{"location":"api/#UnfoldSim.padarray-Tuple{Vector, Tuple, Any}","page":"DocStrings","title":"UnfoldSim.padarray","text":"Pads array with specified value, length padarray(arr, len, val)\n\n\n\n\n\n","category":"method"},{"location":"api/#UnfoldSim.predef_2x2-Tuple{Any}","page":"DocStrings","title":"UnfoldSim.predef_2x2","text":"todo\n\nCareful if you modify nitems with nsubjects = 1, n_items has to be a multiple of 4 (or your equivalent conditions factorial, e.g. all combinations length)\n\n\n\n\n\n","category":"method"},{"location":"api/#UnfoldSim.predef_eeg-Tuple{Any}","page":"DocStrings","title":"UnfoldSim.predef_eeg","text":"Generates a P1/N1/P3 complex.\n\nDefault params:   n_repeats=100   tableModifyFun = x->shuffle(deepcopy(rng),x # random trial order\n\ncomponent / signal\n\nsfreq = 100,   p1 = (p100(;sfreq=sfreq), @formula(0~1),[5],Dict()), # P1 amp 5, no effects   n1 = (n170(;sfreq=sfreq), @formula(0~1+condition),[5,-3],Dict()), # N1 amp 5, dummycoded condition effect (levels \"car\", \"face\") of -3   p3 = (p300(;sfreq=sfreq), @formula(0~1+continuous),[5,1],Dict()), # P3 amp 5, continuous effect range [-5,5] with slope 1\n\nnoise\n\nnoiselevel = 0.2,   noise = PinkNoise(;noiselevel=noiselevel),\n\nonset\n\noverlap = (0.5,0.2), # offset + width/length of Uniform noise. put offset to 1 for no overlap. put width to 0 for no jitter   onset=UniformOnset(;offset=sfreq0.5overlap[1],width=sfreq0.5overlap[2]), \n\n\n\n\n\n","category":"method"},{"location":"api/#UnfoldSim.simulate-NTuple{5, Any}","page":"DocStrings","title":"UnfoldSim.simulate","text":"Simulate eeg data given a simulation design, effect sizes and variances\n\nmake use of return_epoched=true to skip the Onset-calculation + conversion to continuous data and get the epoched data directly\n\n\n\n\n\n","category":"method"},{"location":"api/#UnfoldSim.simulate-Tuple{Any, AbstractComponent, Simulation}","page":"DocStrings","title":"UnfoldSim.simulate","text":"by default call simulate with ::Abstractcomponent,::AbstractDesign`, but allow for custom types\n\nmaking use of other information in simulation\n\n\n\n\n\n","category":"method"},{"location":"api/#UnfoldSim.simulate-Tuple{Any, LinearModelComponent, AbstractDesign}","page":"DocStrings","title":"UnfoldSim.simulate","text":"simulate a linearModel\n\njulia> c = UnfoldSim.LinearModelComponent([0,1,1,0],@formula(0~1+cond),[1,2],Dict()) julia> design = MultiSubjectDesign(;nsubjects=2,nitems=50,item_between=(;:cond=>[\"A\",\"B\"])) julia> simulate(StableRNG(1),c,design)\n\n\n\n\n\n","category":"method"},{"location":"api/#UnfoldSim.simulate-Tuple{Any, MixedModelComponent, AbstractDesign}","page":"DocStrings","title":"UnfoldSim.simulate","text":"simulate MixedModelComponent\n\njulia> design = MultiSubjectDesign(;nsubjects=2,nitems=50,item_between=(;:cond=>[\"A\",\"B\"])) julia> c = UnfoldSim.MixedModelComponent([0.,1,1,0],@formula(dv~1+cond+(1|subject)),[1,2],Dict(:subject=>[2],),Dict()) julia> simulate(StableRNG(1),c,design)\n\n\n\n\n\n","category":"method"},{"location":"api/#UnfoldSim.simulate-Tuple{Any, Vector{<:AbstractComponent}, Simulation}","page":"DocStrings","title":"UnfoldSim.simulate","text":"Simulates erp data given the specified parameters \n\n\n\n\n\n","category":"method"},{"location":"api/#UnfoldSim.weight_σs-Tuple{Dict, Float64, Float64}","page":"DocStrings","title":"UnfoldSim.weight_σs","text":"Weights a σs Dict for MixedModels.jl by a Float64\n\nFinally sales it by σ_lmm, as a trick to simulate noise-free LMMs\n\nI anticipate a function     function weight_σs(σs::Dict,b_σs::Dict,σ_lmm::Float64) where each σs entry can be weighted individually\n\n\n\n\n\n","category":"method"},{"location":"literate/reference/basistypes/","page":"ComponentBasisTypes","title":"ComponentBasisTypes","text":"EditURL = \"https://github.com/unfoldtoolbox/UnfoldSim.jl/blob/main/docs/src/literate/reference/basistypes.jl\"","category":"page"},{"location":"literate/reference/basistypes/","page":"ComponentBasisTypes","title":"ComponentBasisTypes","text":"using UnfoldSim\nusing CairoMakie\nusing DSP\nusing StableRNGs\n#","category":"page"},{"location":"literate/reference/basistypes/#Basistypes","page":"ComponentBasisTypes","title":"Basistypes","text":"","category":"section"},{"location":"literate/reference/basistypes/","page":"ComponentBasisTypes","title":"ComponentBasisTypes","text":"There are several bases types directly implemented. They can be easily used for the components.","category":"page"},{"location":"literate/reference/basistypes/","page":"ComponentBasisTypes","title":"ComponentBasisTypes","text":"note: Note\nYou can use any arbitrary shape defined by yourself! We often make use of hanning(50) from the DSP.jl package","category":"page"},{"location":"literate/reference/basistypes/#EEG","page":"ComponentBasisTypes","title":"EEG","text":"","category":"section"},{"location":"literate/reference/basistypes/","page":"ComponentBasisTypes","title":"ComponentBasisTypes","text":"by default, the EEG bases assume a sampling rate of 100, which can easily be changed by e.g. p100(;sfreq=300)","category":"page"},{"location":"literate/reference/basistypes/","page":"ComponentBasisTypes","title":"ComponentBasisTypes","text":"f = Figure()\nax = f[1,1] = Axis(f)\nfor b in [p100,n170,p300,n400]\n    lines!(ax,b(),label=string(b))\n    scatter!(ax,b(),label=string(b))\nend\naxislegend(ax,merge=true)\nf","category":"page"},{"location":"literate/reference/basistypes/#fMRI","page":"ComponentBasisTypes","title":"fMRI","text":"","category":"section"},{"location":"literate/reference/basistypes/","page":"ComponentBasisTypes","title":"ComponentBasisTypes","text":"default hrf TR is 1. Get to know all your favourite shapes!","category":"page"},{"location":"literate/reference/basistypes/","page":"ComponentBasisTypes","title":"ComponentBasisTypes","text":"##--\nf = Figure()\nplotConfig = (:peak=>1:3:10,\n             :psUnder=>10:5:30,\n             :amplitude=>2:5,\n             :shift=>0:3:10,\n             :peak_width => 0.1:0.5:1.5,\n             :psUnder_width => 0.1:0.5:1.5,\n             )\n\nfor (ix,pl) = enumerate(plotConfig)\n    col = (ix-1)%3 +1\n    row = Int(ceil(ix/3))\n\n    ax = f[row,col] = Axis(f)\n    cfg = collect(pl)\n    for k = cfg[2]\n        lines!(ax,UnfoldSim.hrf(;TR=0.1,(cfg[1]=>k,)...),label=string(k))\n    end\n\n    axislegend(string(cfg[1]);merge=true,)\nend\nf\n\n# pupil","category":"page"},{"location":"literate/reference/basistypes/","page":"ComponentBasisTypes","title":"ComponentBasisTypes","text":"we use the simplified PuRF from Hoeks & Levelt, 1993. Note that https://www.science.org/doi/10.1126/sciadv.abi9979 show some evidence in their supplementary material, that the convolution model is not fully applicable.","category":"page"},{"location":"literate/reference/basistypes/","page":"ComponentBasisTypes","title":"ComponentBasisTypes","text":"f = Figure()\nplotConfig = (:n=>5:3:15,\n             :tmax=>0.5:0.2:1.1,\n             )\n\nfor (ix,pl) = enumerate(plotConfig)\n    ax = f[1,ix] = Axis(f)\n    cfg = collect(pl)\n    for k = cfg[2]\n        lines!(ax,UnfoldSim.PuRF(;(cfg[1]=>k,)...),label=string(k))\n    end\n\n    axislegend(string(cfg[1]);merge=true,)\nend\nf","category":"page"},{"location":"literate/reference/basistypes/","page":"ComponentBasisTypes","title":"ComponentBasisTypes","text":"","category":"page"},{"location":"literate/reference/basistypes/","page":"ComponentBasisTypes","title":"ComponentBasisTypes","text":"This page was generated using Literate.jl.","category":"page"},{"location":"literate/reference/noisetypes/","page":"NoiseTypes","title":"NoiseTypes","text":"EditURL = \"https://github.com/unfoldtoolbox/UnfoldSim.jl/blob/main/docs/src/literate/reference/noisetypes.jl\"","category":"page"},{"location":"literate/reference/noisetypes/","page":"NoiseTypes","title":"NoiseTypes","text":"using UnfoldSim\nusing CairoMakie\nusing DSP\nusing StableRNGs\nimport StatsBase.autocor","category":"page"},{"location":"literate/reference/noisetypes/#What's-the-noise?","page":"NoiseTypes","title":"What's the noise?","text":"","category":"section"},{"location":"literate/reference/noisetypes/","page":"NoiseTypes","title":"NoiseTypes","text":"There are several noise-types directly implemented. Here is a comparison","category":"page"},{"location":"literate/reference/noisetypes/","page":"NoiseTypes","title":"NoiseTypes","text":"f = Figure()\nax_sig = f[1,1:2] = Axis(f;title=\"1.000 samples of noise\")\nax_spec = f[2,1] = Axis(f;title=\"Welch Periodigram\")\nax_auto = f[2,2] = Axis(f;title=\"Autocorrelogram (every 10th lag)\")\nfor n = [PinkNoise RedNoise WhiteNoise NoNoise ExponentialNoise]\n\n    # generate\n    noisevec = gen_noise(StableRNG(1),n(),10000)\n\n    # plot 1000 samples\n    lines!(ax_sig,noisevec[1:1000];label=string(n))\n\n    # calc spectrum\n    perio = welch_pgram(noisevec)\n\n    # plot spectrum\n    lines!(ax_spec,freq(perio),log10.(power(perio)))\n\n    lags = 0:10:500\n    autocor_vec = autocor(noisevec,lags)\n    lines!(ax_auto,lags,autocor_vec)\n\nend\nf[1:2,3] = Legend(f,ax_sig,\"NoiseType\")\nf","category":"page"},{"location":"literate/reference/noisetypes/","page":"NoiseTypes","title":"NoiseTypes","text":"!!! Recommendation    We recommed for smaller signals the ExponentialNoise, maybe with a removed DC offset or a HighPass filter. For long signals, this Noise requires lot's of memory though. maybe Pinknoise is a better choice","category":"page"},{"location":"literate/reference/noisetypes/","page":"NoiseTypes","title":"NoiseTypes","text":"","category":"page"},{"location":"literate/reference/noisetypes/","page":"NoiseTypes","title":"NoiseTypes","text":"This page was generated using Literate.jl.","category":"page"},{"location":"literate/tutorials/quickstart/","page":"Quickstart","title":"Quickstart","text":"EditURL = \"https://github.com/unfoldtoolbox/UnfoldSim.jl/blob/main/docs/src/literate/tutorials/quickstart.jl\"","category":"page"},{"location":"literate/tutorials/quickstart/","page":"Quickstart","title":"Quickstart","text":"using UnfoldSim\nusing Random\nusing CairoMakie","category":"page"},{"location":"literate/tutorials/quickstart/","page":"Quickstart","title":"Quickstart","text":"tipp: Tipp\nUse subtypes(AbstractNoise) (or subtypes(AbstractComponent) etc.) to find already implemented building blocks","category":"page"},{"location":"literate/tutorials/quickstart/#\"Experimental\"-Design","page":"Quickstart","title":"\"Experimental\" Design","text":"","category":"section"},{"location":"literate/tutorials/quickstart/","page":"Quickstart","title":"Quickstart","text":"Define a 1 x 2 design with 20 trials. That is, one condition (condaA) with two levels.","category":"page"},{"location":"literate/tutorials/quickstart/","page":"Quickstart","title":"Quickstart","text":"design = SingleSubjectDesign(;\n        conditions=Dict(:condA=>[\"levelA\",\"levelB\"])\n        ) |> x->RepeatDesign(x,10);\nnothing #hide","category":"page"},{"location":"literate/tutorials/quickstart/#Component-/-Signal","page":"Quickstart","title":"Component / Signal","text":"","category":"section"},{"location":"literate/tutorials/quickstart/","page":"Quickstart","title":"Quickstart","text":"Define a simple component and ground truth simulation formula. Akin to ERP components, we call one simulation signal a component.","category":"page"},{"location":"literate/tutorials/quickstart/","page":"Quickstart","title":"Quickstart","text":"highlight: Highlight\nYou could easily specify multiple components by providing a vector of components, which are automatically added at the same onsets. This procedure simplifies to generate some response that is independent of simulated condition, whereas other depends on it.","category":"page"},{"location":"literate/tutorials/quickstart/","page":"Quickstart","title":"Quickstart","text":"signal = LinearModelComponent(;\n        basis=[0,0,0,0.5,1,1,0.5,0,0],\n        formula = @formula(0~1+condA),\n        β = [1,0.5]\n        );\nnothing #hide","category":"page"},{"location":"literate/tutorials/quickstart/#Onsets-and-Noise","page":"Quickstart","title":"Onsets and Noise","text":"","category":"section"},{"location":"literate/tutorials/quickstart/","page":"Quickstart","title":"Quickstart","text":"We will start with a uniform (but overlapping, offset < length(signal.basis)) onset-distribution","category":"page"},{"location":"literate/tutorials/quickstart/","page":"Quickstart","title":"Quickstart","text":"onset = UniformOnset(;width=20,offset=4);\nnothing #hide","category":"page"},{"location":"literate/tutorials/quickstart/","page":"Quickstart","title":"Quickstart","text":"And we will use some noise","category":"page"},{"location":"literate/tutorials/quickstart/","page":"Quickstart","title":"Quickstart","text":"noise = PinkNoise(;noiselevel=0.2);\nnothing #hide","category":"page"},{"location":"literate/tutorials/quickstart/#Combine-and-Generate","page":"Quickstart","title":"Combine & Generate","text":"","category":"section"},{"location":"literate/tutorials/quickstart/","page":"Quickstart","title":"Quickstart","text":"We will put it all together in one Simulation type","category":"page"},{"location":"literate/tutorials/quickstart/","page":"Quickstart","title":"Quickstart","text":"simulation = Simulation(design, signal,  onset, noise);\nnothing #hide","category":"page"},{"location":"literate/tutorials/quickstart/","page":"Quickstart","title":"Quickstart","text":"finally, we will simulate some data","category":"page"},{"location":"literate/tutorials/quickstart/","page":"Quickstart","title":"Quickstart","text":"data,events = simulate(MersenneTwister(1),simulation);\nnothing #hide","category":"page"},{"location":"literate/tutorials/quickstart/","page":"Quickstart","title":"Quickstart","text":"Data is a n-sample Vector (but could be a Matrix for e.g. MultiSubjectDesign).","category":"page"},{"location":"literate/tutorials/quickstart/","page":"Quickstart","title":"Quickstart","text":"events is a DataFrame that contains a column latency with the onsets of events.","category":"page"},{"location":"literate/tutorials/quickstart/#Plot-them!","page":"Quickstart","title":"Plot them!","text":"","category":"section"},{"location":"literate/tutorials/quickstart/","page":"Quickstart","title":"Quickstart","text":"lines(data;color=\"black\")\nvlines!(events.latency;color=[\"orange\",\"teal\"][1 .+ (events.condA.==\"levelB\")])\ncurrent_figure()","category":"page"},{"location":"literate/tutorials/quickstart/","page":"Quickstart","title":"Quickstart","text":"","category":"page"},{"location":"literate/tutorials/quickstart/","page":"Quickstart","title":"Quickstart","text":"This page was generated using Literate.jl.","category":"page"},{"location":"literate/HowTo/newDesign/","page":"New Experimental Design","title":"New Experimental Design","text":"EditURL = \"https://github.com/unfoldtoolbox/UnfoldSim.jl/blob/main/docs/src/literate/HowTo/newDesign.jl\"","category":"page"},{"location":"literate/HowTo/newDesign/","page":"New Experimental Design","title":"New Experimental Design","text":"using UnfoldSim\nusing StableRNGs\nusing DataFrames\nusing Parameters","category":"page"},{"location":"literate/HowTo/newDesign/#Define-a-new-Design","page":"New Experimental Design","title":"Define a new Design","text":"","category":"section"},{"location":"literate/HowTo/newDesign/","page":"New Experimental Design","title":"New Experimental Design","text":"A design specifies how much data is generated, and how the event-table(s) should be generated. Already implemented examples are MultiSubjectDesign and SingleSubjectdesign","category":"page"},{"location":"literate/HowTo/newDesign/","page":"New Experimental Design","title":"New Experimental Design","text":"We need 3 things for a new design: a struct<:AbstractDesign, a size and a generate function","category":"page"},{"location":"literate/HowTo/newDesign/#)-type","page":"New Experimental Design","title":"1) type","text":"","category":"section"},{"location":"literate/HowTo/newDesign/","page":"New Experimental Design","title":"New Experimental Design","text":"We need a ImbalanceSubjectDesign struct. You are free to implement it as you wish, as long as the other two functions are implemented","category":"page"},{"location":"literate/HowTo/newDesign/","page":"New Experimental Design","title":"New Experimental Design","text":"@with_kw struct ImbalanceSubjectDesign <: UnfoldSim.AbstractDesign\n    nTrials::Int\n    balance::Float64 = 0.5 # default balanced\nend;\nnothing #hide","category":"page"},{"location":"literate/HowTo/newDesign/#)-size","page":"New Experimental Design","title":"2) size","text":"","category":"section"},{"location":"literate/HowTo/newDesign/","page":"New Experimental Design","title":"New Experimental Design","text":"we need a size(design::ImbalanceSubjectDesign) function to tell how many events we will have. This is used at different places, e.g. in the Default onset implementation","category":"page"},{"location":"literate/HowTo/newDesign/","page":"New Experimental Design","title":"New Experimental Design","text":"# note the trailling , to make it a Tuple\nsize(design::ImbalanceSubjectDesign) = (design.nTrials,);\nnothing #hide","category":"page"},{"location":"literate/HowTo/newDesign/#)-generate","page":"New Experimental Design","title":"3) generate","text":"","category":"section"},{"location":"literate/HowTo/newDesign/","page":"New Experimental Design","title":"New Experimental Design","text":"We need a type generate(design::ImbalanceSubjectDesign) function. This function should return the actual table as a DataFrame","category":"page"},{"location":"literate/HowTo/newDesign/","page":"New Experimental Design","title":"New Experimental Design","text":"function generate(design::ImbalanceSubjectDesign)\n    nA = Int(round.(design.nTrials .* design.balance))\n    nB = Int(round.(design.nTrials .* (1-design.balance)))\n    @assert nA + nB  ≈ design.nTrials\n    levels = vcat(repeat([\"levelA\"],nA),repeat([\"levelB\"],nB))\n    return DataFrame(Dict(:condition=>levels))\nend;\nnothing #hide","category":"page"},{"location":"literate/HowTo/newDesign/","page":"New Experimental Design","title":"New Experimental Design","text":"Finally, we can test the function and see whether it returns a Design-DataFrame as we requested","category":"page"},{"location":"literate/HowTo/newDesign/","page":"New Experimental Design","title":"New Experimental Design","text":"design = ImbalanceSubjectDesign(;nTrials=6,balance=0.2)\ngenerate(design)","category":"page"},{"location":"literate/HowTo/newDesign/","page":"New Experimental Design","title":"New Experimental Design","text":"important: Important\nit is the users task to ensure that each run is reproducible. So if you have a random process (e.g. shuffling), be sure to   safe a RNG object in your struct and use it in your generate function.","category":"page"},{"location":"literate/HowTo/newDesign/","page":"New Experimental Design","title":"New Experimental Design","text":"","category":"page"},{"location":"literate/HowTo/newDesign/","page":"New Experimental Design","title":"New Experimental Design","text":"This page was generated using Literate.jl.","category":"page"},{"location":"literate/tutorials/poweranalysis/","page":"Poweranalysis","title":"Poweranalysis","text":"EditURL = \"https://github.com/unfoldtoolbox/UnfoldSim.jl/blob/main/docs/src/literate/tutorials/poweranalysis.jl\"","category":"page"},{"location":"literate/tutorials/poweranalysis/","page":"Poweranalysis","title":"Poweranalysis","text":"using UnfoldSim\nusing Unfold\nusing Statistics\nusing HypothesisTests\nusing DataFrames\nusing Random","category":"page"},{"location":"literate/tutorials/poweranalysis/#Simple-Poweranalysis-Script","page":"Poweranalysis","title":"Simple Poweranalysis Script","text":"","category":"section"},{"location":"literate/tutorials/poweranalysis/","page":"Poweranalysis","title":"Poweranalysis","text":"For a power analysis, we will repeatedly simulate data, and check whether we can find a significant effect.","category":"page"},{"location":"literate/tutorials/poweranalysis/","page":"Poweranalysis","title":"Poweranalysis","text":"We perform the power analysis on epoched data.","category":"page"},{"location":"literate/tutorials/poweranalysis/","page":"Poweranalysis","title":"Poweranalysis","text":"pvals = fill(NaN,100)\n@time for seed = eachindex(pvals)\n    # Simulate data of 30 subjects\n    data,evts = UnfoldSim.predef_2x2(MersenneTwister(seed);\n                n_subjects=20, ## 30 subjects\n                overlap=(1,0), ## deactivate overlap\n                noiselevel=10,  ## add more noise to make it more challenging\n                return_epoched=true, ## saves us the epoching step\n                )\n\n\n    # take the mean over a pre-specified timewindow\n    evts.y = dropdims(mean(data[40:60,:],dims=1),dims=(1))\n\n    # extract the two levels of condition A\n    evts_reduced = combine(groupby(evts,[:subject,:A]),:y=>mean)\n    y_big = evts_reduced[evts_reduced.A .==\"a_big\",:y_mean]\n    y_small = evts_reduced[evts_reduced.A .==\"a_small\",:y_mean]\n\n    # calculate a one-sided t-test\n    pvals[seed] = pvalue(OneSampleTTest(y_big,y_small))\nend","category":"page"},{"location":"literate/tutorials/poweranalysis/","page":"Poweranalysis","title":"Poweranalysis","text":"let's calculate the power","category":"page"},{"location":"literate/tutorials/poweranalysis/","page":"Poweranalysis","title":"Poweranalysis","text":"power = mean(pvals .<0.05)*100","category":"page"},{"location":"literate/tutorials/poweranalysis/","page":"Poweranalysis","title":"Poweranalysis","text":"","category":"page"},{"location":"literate/tutorials/poweranalysis/","page":"Poweranalysis","title":"Poweranalysis","text":"This page was generated using Literate.jl.","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = UnfoldSim","category":"page"},{"location":"#UnfoldSim","page":"Home","title":"UnfoldSim","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for UnfoldSim.","category":"page"},{"location":"#Start-simulating-timeseries","page":"Home","title":"Start simulating timeseries","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"We offer some predefined signals, check them out!","category":"page"},{"location":"","page":"Home","title":"Home","text":"For instance an P1/N170/P300 complex.","category":"page"},{"location":"","page":"Home","title":"Home","text":"using UnfoldSim\nusing CairoMakie\ndata,evts = UnfoldSim.predef_eeg(;n_repeats=1,noiselevel=0.8)\n\nlines(data;color=\"black\")\nvlines!(evts.latency;color=[\"orange\",\"teal\"][1 .+ (evts.condition .==\"car\")])\n\ncurrent_figure()","category":"page"},{"location":"#Or-simulate-epoched-data-directly","page":"Home","title":"Or simulate epoched data directly","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"\ndata,evts = UnfoldSim.predef_eeg(;n_repeats=20,noiselevel=0.8,return_epoched=true)\nheatmap(data[:,sortperm(evts,[:condition,:continuous])])\n","category":"page"},{"location":"literate/HowTo/repeatTrials/","page":"Repeating Trials within a Design","title":"Repeating Trials within a Design","text":"EditURL = \"https://github.com/unfoldtoolbox/UnfoldSim.jl/blob/main/docs/src/literate/HowTo/repeatTrials.jl\"","category":"page"},{"location":"literate/HowTo/repeatTrials/","page":"Repeating Trials within a Design","title":"Repeating Trials within a Design","text":"using UnfoldSim","category":"page"},{"location":"literate/HowTo/repeatTrials/#Repeating-Design-entries","page":"Repeating Trials within a Design","title":"Repeating Design entries","text":"","category":"section"},{"location":"literate/HowTo/repeatTrials/","page":"Repeating Trials within a Design","title":"Repeating Trials within a Design","text":"Sometimes we want to repeat a design, that is, have multiple trials with identical values, but it is not always straight forward to implement For instance, there is no way to easily modify MultiSubjectDesign to have multiple identical subject/item combinations, without doing awkward repetitions of condition-levels or something.","category":"page"},{"location":"literate/HowTo/repeatTrials/","page":"Repeating Trials within a Design","title":"Repeating Trials within a Design","text":"If you struggle with this problem RepeatDesign is an easy tool for you","category":"page"},{"location":"literate/HowTo/repeatTrials/","page":"Repeating Trials within a Design","title":"Repeating Trials within a Design","text":"designOnce = MultiSubjectDesign(;\n        n_items=2,\n\t\tn_subjects = 2,\n        subjects_between =Dict(:cond=>[\"levelA\",\"levelB\"]),\n\t\titems_between =Dict(:cond=>[\"levelA\",\"levelB\"]),\n        );\n\ndesign = RepeatDesign(designOnce,4);\ngenerate(design)","category":"page"},{"location":"literate/HowTo/repeatTrials/","page":"Repeating Trials within a Design","title":"Repeating Trials within a Design","text":"As you can see, the design was simply repeated. As always, you can ignore the dv column, it is for internal consistency with MixedModelsSim.jl","category":"page"},{"location":"literate/HowTo/repeatTrials/","page":"Repeating Trials within a Design","title":"Repeating Trials within a Design","text":"note: Note\nif you implemented your own AbstractDesign, you need to define the size function accordingly. E.g.:   Base.size(design::RepeatDesign{SingleSubjectDesign}) = size(design.design).*design.repeat","category":"page"},{"location":"literate/HowTo/repeatTrials/","page":"Repeating Trials within a Design","title":"Repeating Trials within a Design","text":"","category":"page"},{"location":"literate/HowTo/repeatTrials/","page":"Repeating Trials within a Design","title":"Repeating Trials within a Design","text":"This page was generated using Literate.jl.","category":"page"},{"location":"literate/tutorials/simulateERP/","page":"Simulate ERPs","title":"Simulate ERPs","text":"EditURL = \"https://github.com/unfoldtoolbox/UnfoldSim.jl/blob/main/docs/src/literate/tutorials/simulateERP.jl\"","category":"page"},{"location":"literate/tutorials/simulateERP/","page":"Simulate ERPs","title":"Simulate ERPs","text":"using UnfoldSim\nusing CairoMakie\nusing Random\nusing Unfold\nusing UnfoldMakie","category":"page"},{"location":"literate/tutorials/simulateERP/#ERP-Complex","page":"Simulate ERPs","title":"ERP Complex","text":"","category":"section"},{"location":"literate/tutorials/simulateERP/","page":"Simulate ERPs","title":"Simulate ERPs","text":"here we will learn how to simulate a typical ERP complex with P100, N170, P300","category":"page"},{"location":"literate/tutorials/simulateERP/","page":"Simulate ERPs","title":"Simulate ERPs","text":"let's grab a SingleSubjectDesign and add a continuous predictor","category":"page"},{"location":"literate/tutorials/simulateERP/","page":"Simulate ERPs","title":"Simulate ERPs","text":"design = SingleSubjectDesign(;\n        conditions=Dict(:condition=>[\"car\",\"face\"],:continuous=>range(-5,5,length=10))\n        ) |> x->RepeatDesign(x,100);\nnothing #hide","category":"page"},{"location":"literate/tutorials/simulateERP/","page":"Simulate ERPs","title":"Simulate ERPs","text":"let's make use of the prespecified basis functions, but use different formulas + parameters for each!","category":"page"},{"location":"literate/tutorials/simulateERP/","page":"Simulate ERPs","title":"Simulate ERPs","text":"p100 is unaffected by our design and has amplitude of 5","category":"page"},{"location":"literate/tutorials/simulateERP/","page":"Simulate ERPs","title":"Simulate ERPs","text":"p1 =  LinearModelComponent(;\n        basis = p100(),\n        formula = @formula(0~1),\n        β = [5]\n        );\nnothing #hide","category":"page"},{"location":"literate/tutorials/simulateERP/","page":"Simulate ERPs","title":"Simulate ERPs","text":"n170 has a condition effect, faces are more negative than cars","category":"page"},{"location":"literate/tutorials/simulateERP/","page":"Simulate ERPs","title":"Simulate ERPs","text":"n1 =  LinearModelComponent(;\n        basis = n170(),\n        formula = @formula(0~1+condition),\n        β = [5,-3]\n        );\nnothing #hide","category":"page"},{"location":"literate/tutorials/simulateERP/","page":"Simulate ERPs","title":"Simulate ERPs","text":"p300 has a continuous effect, higher continuous values will result in larger P300's","category":"page"},{"location":"literate/tutorials/simulateERP/","page":"Simulate ERPs","title":"Simulate ERPs","text":"p3 =  LinearModelComponent(;\n        basis = p300(),\n        formula = @formula(0~1+continuous),\n        β = [5,1]\n        );\nnothing #hide","category":"page"},{"location":"literate/tutorials/simulateERP/","page":"Simulate ERPs","title":"Simulate ERPs","text":"now we can simply combine the components and simulate","category":"page"},{"location":"literate/tutorials/simulateERP/","page":"Simulate ERPs","title":"Simulate ERPs","text":"components = [p1,n1,p3]\ndata,evts = simulate(MersenneTwister(1),design,[p1,n1,p3],UniformOnset(;width=0,offset=1000),PinkNoise());\n\n# Analysis","category":"page"},{"location":"literate/tutorials/simulateERP/","page":"Simulate ERPs","title":"Simulate ERPs","text":"Let's check that everything worked out well, by using Unfold","category":"page"},{"location":"literate/tutorials/simulateERP/","page":"Simulate ERPs","title":"Simulate ERPs","text":"m = fit(UnfoldModel,Dict(Any=>(@formula(0~1+condition+continuous),firbasis(τ=[-0.1,1],sfreq=100,name=\"basis\"))),evts,data);\nnothing #hide","category":"page"},{"location":"literate/tutorials/simulateERP/","page":"Simulate ERPs","title":"Simulate ERPs","text":"first the \"pure\" beta/linear regression parameters","category":"page"},{"location":"literate/tutorials/simulateERP/","page":"Simulate ERPs","title":"Simulate ERPs","text":"plot_erp(coeftable(m))","category":"page"},{"location":"literate/tutorials/simulateERP/","page":"Simulate ERPs","title":"Simulate ERPs","text":"and now beautifully visualized as marginal betas / predicted ERPs","category":"page"},{"location":"literate/tutorials/simulateERP/","page":"Simulate ERPs","title":"Simulate ERPs","text":"plot_erp(effects(Dict(:condition=>[\"car\",\"face\"],:continuous=>-5:5),m);\n        setMappingValues=(:color=>:continuous,linestyle=:condition,group=:continuous),\n        setExtraValues=(;categoricalColor=false))","category":"page"},{"location":"literate/tutorials/simulateERP/","page":"Simulate ERPs","title":"Simulate ERPs","text":"","category":"page"},{"location":"literate/tutorials/simulateERP/","page":"Simulate ERPs","title":"Simulate ERPs","text":"This page was generated using Literate.jl.","category":"page"}]
}
