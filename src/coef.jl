
get_df(;kwargs...) = get_df(MersenneTwister();kwargs...) # without seed
function get_df(rng;n_subjects=20,n_items=10)
    # simulate a 2x2 dataset with 4x n_items per subject
    df = simdat_crossed(rng,n_subjects,n_items;both_win=Dict(:catA => ["1","2"],:catB=>["-1","-2"]))
    return df
end

function signal_coef(df)
    β  = [1,0.5,-0.2]
    σs = [0.5, 0.5 ,0.5]
    ρ  = [1.   0.8 0.0;
          0.8  1  -0.5;
          0.0 -0.5 1]
σ = 1
θ_subjects = create_re(σs...,corrmat = ρ)
θ_items = create_re(1.)
f1 = @formula(dv ~ 1 + catA+catB + (1|item) + (1+catA+catB|subj));
m1 = fit(MixedModel, f1, df)
MixedModelsSim.update!(m1;item=θ_items,subj=θ_subjects)
refit!(simulate!(m1,β=β,σ=σ))
return ranef(m1) # or is coef the right call here?
end