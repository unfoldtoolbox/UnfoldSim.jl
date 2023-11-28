# here we can define some predefined Simulations. Convenient if you just want to have a quick simulation :)



predef_2x2(;kwargs...) = predef_2x2(MersenneTwister(1);kwargs...) # without rng always call same one
predef_eeg(;kwargs...) = predef_eeg(MersenneTwister(1);kwargs...) # without rng always call same one
predef_eeg(nsubjects::Int;kwargs...) = predef_eeg(MersenneTwister(1),nsubjects;kwargs...) # without rng always call same one

"""
predef_eeg(;kwargs...)
predef_eeg(rng;kwargs...)
predef_eeg(rng,n_subjects;kwargs...)

Gene
rates a P1/N1/P3 complex.
In case `n_subjects` is defined - `MixedModelComponents`` are generated, else `LinearModelComponents`

Most used kwargs is: `return_epoched=true` to ignore the overlap/onset bits and return already epoched data

## Default params:

. n_repeats=100
- tableModifyFun = x->shuffle(deepcopy(rng),x # random trial order
- conditions = Dict(...),

#### component / signal
- sfreq = 100,
- p1 = (p100(;sfreq=sfreq), @formula(0~1),[5],Dict()), # P1 amp 5, no effects
- n1 = (n170(;sfreq=sfreq), @formula(0~1+condition),[5,-3],Dict()), # N1 amp 5, dummycoded condition effect (levels "car", "face") of -3
- p3 = (p300(;sfreq=sfreq), @formula(0~1+continuous),[5,1],Dict()), # P3 amp 5, continuous effect range [-5,5] with slope 1

#### noise
- noiselevel = 0.2,
- noise = PinkNoise(;noiselevel=noiselevel),
                
#### onset
- overlap = (0.5,0.2), # offset + width/length of Uniform noise. put offset to 1 for no overlap. put width to 0 for no jitter
- onset=UniformOnset(;offset=sfreq*0.5*overlap[1],width=sfreq*0.5*overlap[2]), 
"""
function predef_eeg(rng;
                    # design
                    n_repeats=100,
                    tableModifyFun = x->shuffle(deepcopy(rng),x),
                    
                    # component / signal
                    sfreq = 100,
                    p1 = (p100(;sfreq=sfreq), @formula(0~1),[5],Dict()),
                    n1 = (n170(;sfreq=sfreq), @formula(0~1+condition),[5,-3],Dict()),
                    p3 = (p300(;sfreq=sfreq), @formula(0~1+continuous),[5,1],Dict()),

                    
                    kwargs...
                    )

    design = SingleSubjectDesign(;
        conditions=Dict(:condition=>["car","face"],
                        :continuous=>range(-5,5,length=10)),
        tableModifyFun = tableModifyFun
        ) |> x->RepeatDesign(x,n_repeats);
        return predef_eeg(rng,design,LinearModelComponent,[p1,n1,p3];sfreq,kwargs...)
end

function predef_eeg(rng::AbstractRNG,design::AbstractDesign,T::Type{<:AbstractComponent},comps;
                    sfreq=100,
                    # noise
                    noiselevel = 0.2,
                    noise = PinkNoise(;noiselevel=noiselevel),
                    
                    # onset
                    overlap = (0.5,0.2),
                    onset=UniformOnset(;offset=sfreq*overlap[1],width=sfreq*overlap[2]), #put offset to 1 for no overlap. put width to 0 for no jitter
                    kwargs...
                    )

    components = []
    for c = comps
        append!(components,[T(c...)])
    end
    return simulate(rng,design, components,  onset, noise;kwargs...);
end

function predef_eeg(rng::AbstractRNG,n_subjects;
                    # design
                    n_items=100,
                    tableModifyFun = x->shuffle(deepcopy(rng),x),
                    
                    conditions = Dict(:condition=>["car","face"],
                    :continuous=>range(-5,5,length=10)),
                    # component / signal
                    sfreq = 100,
                    p1 = (;s=p100(;sfreq=sfreq), f=@formula(0~1+(1|subject)+(1|item)),β=[5],σs=Dict(:subject=>[1],:item=>[1]),c=Dict()),
                    n1 = (;s=n170(;sfreq=sfreq), f=@formula(0~1+condition+ (1+condition|subject)+ (1+condition|item)), β=[5,-3],σs=Dict(:subject=>[1,1],:item=>[0.5,0.5]),c=Dict()),
                    p3 = (;s=p300(;sfreq=sfreq), f=@formula(0~1+continuous+(1+continuous|subject)+(1+continuous|item)),β=[5,1], σs=Dict(:subject=>[1,1],:item=>[0.5,0.5]),c=Dict()),
                    kwargs...)


    design = MultiSubjectDesign(;
        n_subjects = n_subjects,
        n_items = n_items,
        items_between = conditions,
        tableModifyFun = tableModifyFun)

    return predef_eeg(rng,design,MixedModelComponent,[p1,n1,p3];sfreq,kwargs...)
    
end
"""

    predef_2x2(rng::AbstractRNG;kwargs...)

Most used kwargs is: `return_epoched=true` to ignore the overlap/onset bits and return already epoched data

#### design
- `n_items`=100,
- `n_subjects`=1,
- `conditions` = Dict(:A=>["a_small","a_big"],:B=>["b_tiny","b_large"]),
- `tableModifyFun` = x->shuffle(deepcopy(rng),x),

#### component / signal
- `signalsize` = 100, length of simulated hanning window
- `basis`` = hanning(signalsize), the actual "function", `signalsize` is only used here
- `β` = [1,-0.5,.5,+1], the parameters
- `σs` = Dict(:subject=>[1,0.5,0.5,0.5],:item=>[1]), - only in n_subjects>=2 case, specifies the random effects
- `contrasts` = Dict(:A=>EffectsCoding(),:B=>EffectsCoding()) - effect coding by default
- `formula` = n_subjects==1 ? @formula(0~1+A*B) : @formula(dv~1+A*B+(A*B|subject)+(1|item)),

#### noise
- `noiselevel` = 0.2,
- `noise` = PinkNoise(;noiselevel=noiselevel),

#### onset
- `overlap` = (0.5,0.2),
- `onset`=UniformOnset(;offset=signalsize*overlap[1],width=signalsize*overlap[2]), #put offset to 1 for no overlap. put width to 0 for no jitter


Careful if you modify n_items with n_subjects = 1, n_items has to be a multiple of 4 (or your equivalent conditions factorial, e.g. all combinations length)
"""
function predef_2x2(rng::AbstractRNG;
                    # design
                    n_items=100,
                    n_subjects=1,
                    conditions = Dict(:A=>["a_small","a_big"],:B=>["b_tiny","b_large"]),
                    tableModifyFun = x->shuffle(deepcopy(rng),x),
                    
                    # component / signal
                    signalsize = 100,
                    basis = hanning(signalsize), # the component "function"
                    β = [1,-0.5,.5,+1],
                    σs = Dict(:subject=>[1,0.5,0.5,0.5],:item=>[1]), # only in n_subjects>2 case
                    contrasts = Dict(:A=>EffectsCoding(),:B=>EffectsCoding()),
                    
                    formula = n_subjects==1 ? @formula(0~1+A*B) : @formula(dv~1+A*B+(A*B|subject)+(1|item)),
                    
                    # noise
                    noiselevel = 0.2,
                    noise = PinkNoise(;noiselevel=noiselevel),
                    
                    # onset
                    overlap = (0.5,0.2),
                    onset=UniformOnset(;offset=signalsize*overlap[1],width=signalsize*overlap[2]), #put offset to 1 for no overlap. put width to 0 for no jitter
                    kwargs...
                    )
                    
    if n_subjects == 1
        design = SingleSubjectDesign(;
            conditions = conditions,
            tableModifyFun = tableModifyFun) |> x->RepeatDesign(x,n_items./length(x))

        signal = LinearModelComponent(;
            basis=basis,
            formula = formula,
            β = β,
            contrasts = contrasts,
    );
    else
        design = MultiSubjectDesign(;           
                n_subjects = n_subjects,
                n_items = n_items  ,              
                items_between = conditions,
                tableModifyFun = tableModifyFun)
        signal = MixedModelComponent(;
                basis=basis,
                formula = formula,
                β = β,
                σs = σs,

                contrasts = contrasts,
        );
        
    end
 
    data,events = simulate(rng,design, signal,  onset, noise;kwargs...);
    return data,events
end
