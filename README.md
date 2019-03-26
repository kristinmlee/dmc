# dmc
## Distinguishing among modes of convergent adaptation using population genomic data

This repository contains code associated with Lee and Coop (2017).

### MYAdIDAS extension
Code for extension MYAdIDAS by Sivan Yair and Kristin Lee, associated with Oziolor et al. (2019) and modified for cases of adaptive introgression, can be found in [myAdidas](https://github.com/kristinmlee/dmc/tree/master/myAdidas). This also includes modifications for deletions, strong selection coefficients, and differences in effective population sizes among populations. See folder for complete details.

### 1. Scripts to run DMC
+ A fully worked example for how to use DMC can be found at [dmc_example.md](https://github.com/kristinmlee/dmc/blob/master/dmc_example.md)


+ Short descriptions of the scripts needed to run DMC found in this repository.

	See individual files for more details and example for their usage in context.


	+ [calcNeutralF.R](https://github.com/kristinmlee/dmc/blob/master/calcNeutralF.R) contains code to generate and save **F** as an R object.

	+ [genSelMatrices_individualModes.R](https://github.com/kristinmlee/dmc/blob/master/genSelMatrices_individualModes.R) contains functions to generate **F<sup>(S)</sup>** for single modes of convergent adaptation.

	+ [genSelMatrices_multipleModes.R](https://github.com/kristinmlee/dmc/blob/master/genSelMatrices_multipleModes.R) contains functions to generate **F<sup>(S)</sup>** for multiple modes of convergence.

	+ [calcInvDetSelMatrices_all.R](https://github.com/kristinmlee/dmc/blob/master/calcInvDetSelMatrices_all.R) contains scripts to generate all the possible inverses and determinants.

	+ [calcCompositeLike.R](https://github.com/kristinmlee/dmc/blob/master/calcCompositeLike.R) contains functions for calculating the composite log-likelihoods for all models.

	+ [getMCLE.R](https://github.com/kristinmlee/dmc/blob/master/getMCLE.R) contains functions to get the maximum composite likelihood estimates of parameters for all convergence models.


### 2. Scripts associated with the [*Mimulus*](https://github.com/kristinmlee/dmc/tree/master/mimulusAnalysis) and [killifish](https://github.com/kristinmlee/dmc/tree/master/killifishAnalysis) analyses in Lee and Coop (2017)

### 3. [Modified version of mssel](https://github.com/kristinmlee/dmc/tree/master/mssel_modified), a version of ms [1] that allows for the incorporation of selection at a single site, modified by KL to allow for multiple independent origins of beneficial lineages.

[1] Hudson, R. R. (2002). Generating samples under a Wright–Fisher neutral model of genetic variation. Bioinformatics 18, 337–338.

