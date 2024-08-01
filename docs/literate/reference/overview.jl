# # Overview of functionality
# A UnfoldSim simulation has four ingredients: Design, Component, Onset and Noise. Here we provide a short overview of the implemented types.

# ### Setup
# ```@raw html
# <details>
# <summary>Click to expand</summary>
# ```
## Load required packages
using UnfoldSim
using InteractiveUtils
# ```@raw html
# </details>
# ```

# ## Design
# Designs define the experimental design. They can be nested, e.g. `RepeatDesign(SingleSubjectDesign,10)` would repeat the generated design-dataframe 10x.
subtypes(AbstractDesign)

# ## Component
# Components define a signal. Some components can be nested, e.g. `LinearModelComponent|>MultichannelComponent`, see the multi-channel tutorial for more information.
subtypes(AbstractComponent)

# ## Onsets
# Onsets define the distance between events in the continuous signal.
subtypes(AbstractOnset)

# ## Noise
# Choose the noise you need!
subtypes(AbstractNoise)
