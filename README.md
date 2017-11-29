# dmc
## Distinguishing among modes of convergent adaptation using population genomic data

This repository contains code associated with Lee and Coop (2017).
### 1. Scripts to run DMC
#### A fully worked example for how to use DMC can be found at [dmc_example.md](https://github.com/kristinmlee/dmc/blob/master/dmc_example.md)


#### Short descriptions of the scripts needed to run DMC found in this repository.

See individual files for more details and example for their usage in context.


+ [dmc/calcNeutralF.R](https://github.com/kristinmlee/dmc/blob/master/calcNeutralF.R) contains code to generate and save **F** as an R object.

+ [dmc/genSelMatrices_individualModes.R](https://github.com/kristinmlee/dmc/blob/master/genSelMatrices_individualModes.R) contains functions to generate **F^(S)^** for single modes of convergent adaptation.

+ [dmc/genSelMatrices_multipleModes.R](https://github.com/kristinmlee/dmc/blob/master/genSelMatrices_multipleModes.R) contains functions to generate **F^(S)^** for multiple modes of convergence.

+ [dmc/calcInvDetSelMatrices_all.R](https://github.com/kristinmlee/dmc/blob/master/calcInvDetSelMatrices_all.R) contains scripts to generate all the possible inverses and determinants.

+ [dmc/calcCompositeLike.R](https://github.com/kristinmlee/dmc/blob/master/calcCompositeLike.R) contains functions for calculating the composite log-likelihoods for all models.

+ [dmc/getMCLE.R](https://github.com/kristinmlee/dmc/blob/master/getMCLE.R) contains functions to get the maximum composite likelihood estimates of parameters for all convergence models.


### 2. Scripts associated with the [*Mimulus*](https://github.com/kristinmlee/dmc/blob/mimulusAnalysis) and [killifish](https://github.com/kristinmlee/dmc/blob/killifishAnalysis) analyses in Lee and Coop (2017)

### 3. [Modified version of mssel](https://github.com/kristinmlee/dmc/blob/mssel_modified), a version of ms [1] that allows for the incorporation of selection at a single site, modified by KL to allow for multiple independent origins of beneficial lineages.

[1] Hudson, R. R. (2002). Generating samples under a Wright–Fisher neutral model of genetic variation. Bioinformatics 18, 337–338.

