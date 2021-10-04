using MixedModelsSim
using MixedModels
using Random

function get_df(seed)
n_subjects = 20
n_items = 10
# 3 predictors, 1 + catA + catB


coef = simdat_crossed(MersenneTwister(seed),n_subjects,n_items;both_win=Dict(:catA => ["1","2"],:catB=>["-1","-2"]))
return coef
end

dfA = get_df(1)
dfB = get_df(2)


function signal_coef(df)
β = [1,0.5,-0.2]
σs = [0.5, 0.5 ,0.5]
ρ = [1.  0.8 0.0;
     0.8  1  -0.5;
     0.0 -0.5 1]
σ = 1
θ_subjects = create_re(σs...,corrmat = ρ)
θ_items = create_re(1.)
f1 = @formula(dv ~ 1 + catA+catB + (1|item) + (1+catA+catB|subj));
m1 = fit(MixedModel, f1, df)
MixedModelsSim.update!(m1;item=θ_items,subj=θ_subjects)
refit!(simulate!(m1,β=β,σ=σ))


end


