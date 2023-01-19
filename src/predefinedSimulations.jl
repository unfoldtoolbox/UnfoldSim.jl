# here we can define some predefined Simulations. Convenient if you just want to have a quick simulation :)

predef_2x2(;kwargs...) = predef_2x2(MersenneTwister(1);kwargs...) # without rng always call same one

function predef_2x2(rng;
                    # design
                    n_trials=100,
                    n_subjects=1,
                    conditions = Dict(:A=>["a_small","a_big"],:B=>["b_tiny","b_large"]),
                    tableModifyFun = x->shuffle(rng,x),
                    
                    # component / signal
                    signalsize = 100,
                    basis = hanning(signalsize), # the component "function"
                    β = [1,-0.5,.5,+1],
                    σs = Dict(:subject=>[1.,0.5,0.5,0.5],:item=>[1.]), # only in n_subjects>2 case
                    contrasts = Dict(:A=>EffectsCoding(),:B=>EffectsCoding()),
                    
                    formula = n_subjects==1 ? @formula(0~1+A*B) : @formula(dv~1+A*B+(A*B|subject)+(1|item)),
                    
                    # noise
                    noiselevel = 0.2,
                    noise = PinkNoise(;noiselevel=noiselevel),
                    
                    # onset
                    overlap = (0.5,0.2),
                    onset=UniformOnset(;offset=signalsize*overlap[1],width=signalsize*overlap[2]), #put offset to 1 for no overlap. put width to 0 for no jitter
                    
                    )
                    
    if n_subjects == 1
        design = SingleSubjectDesign(;
            n_trials=n_trials,
            conditions = conditions,
            tableModifyFun = tableModifyFun)

        signal = LinearModelComponent(;
            basis=basis,
            formula = formula,
            β = β,
            contrasts = contrasts,
    );
    else
        design = MultiSubjectDesign(;           
                n_subj = n_subjects,
                n_item = n_trials  ,              
                item_btwn = conditions,
                tableModifyFun = tableModifyFun)
        signal = MixedModelComponent(;
                basis=basis,
                formula = formula,
                β = β,
                σs = σs,

                contrasts = contrasts,
        );
        
    end
 
    data,events = simulate(rng,design, signal,  onset, noise);
    return data,events
end