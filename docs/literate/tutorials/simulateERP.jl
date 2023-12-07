using UnfoldSim
using CairoMakie
using Random
using Unfold
using UnfoldMakie
# ## ERP Complex
# Here we will learn how to simulate a typical ERP complex with P100, N170, P300.

# Let's grab a SingleSubjectDesign and add a continuous predictor
design = SingleSubjectDesign(;
        conditions=Dict(:condition=>["car","face"],:continuous=>range(-5,5,length=10))
        ) |> x->RepeatDesign(x,100);

# Let's make use of the prespecified basis functions, but use different formulas + parameters for each!

# **p100** is unaffected by our design and has amplitude of 5
p1 =  LinearModelComponent(;
        basis = p100(),
        formula = @formula(0~1),
        β = [5]
        );

# **n170** has a condition effect, faces are more negative than cars
n1 =  LinearModelComponent(;
        basis = n170(),
        formula = @formula(0~1+condition),
        β = [5,-3]
        );
# **p300** has a continuous effect, higher continuous values will result in larger P300's
p3 =  LinearModelComponent(;
        basis = p300(),
        formula = @formula(0~1+continuous),
        β = [5,1]
        );

# Now we can simply combine the components and simulate 
components = [p1,n1,p3] 
data,evts = simulate(MersenneTwister(1),design,[p1,n1,p3],UniformOnset(;width=0,offset=1000),PinkNoise());


# ## Analysis
# Let's check that everything worked out well, by using Unfold

m = fit(UnfoldModel,Dict(Any=>(@formula(0~1+condition+continuous),firbasis(τ=[-0.1,1],sfreq=100,name="basis"))),evts,data);

# first the "pure" beta/linear regression parameters
plot_erp(coeftable(m))

# and now beautifully visualized as marginal betas / predicted ERPs
plot_erp(effects(Dict(:condition=>["car","face"],:continuous=>-5:5),m);
        mapping=(:color=>:continuous,linestyle=:condition,group=:continuous),
        categorical_color=false)
