# # Overview of functionality
# UnfoldSim has many modules, here we try to collect them to provide you with an overview
using UnfoldSim
using InteractiveUtils

# ## Design
# Designs define the experimental design. They can be nested, e.g. `RepeatDesign(SingleSubjectDesign,10)` would repeat the generated design-dataframe 10x.
subtypes(AbstractDesign)

# ## Component
# components define a signal. Some components can be nested, e.g. `LinearModelComponent|>MultichannelComponent`, see the multi-channel tutorial for more information
subtypes(AbstractComponent)

# ## Onsets
# Onsets define the distance between events in the continuous signal
subtypes(AbstractOnset)

# ## Noise
# Choose the noise you need!
subtypes(AbstractNoise)


