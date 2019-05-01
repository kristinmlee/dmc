# **DMC-MYAdIDAS**

## **DMC** **M**odified for **Y**oung **Ad**aptive **I**ntrogression **D**istinguished from **A**llele **S**tanding 

**Sivan Yair and Kristin M. Lee**

This folder contains code associated with Oziolor et al. (2019). In this publication, we extended DMC (introduced in Lee and Coop (2017); see parent directory) to confirm and characterize adaptive introgression from *Fundulus heteroclitus* into *Fundulus grandis*. We ran the method on the AHR deletion on chromosome 1 and the ARNT region on chromosome 10. Scripts (with descriptions) and results for both regions are available in the subdirectories: [Chr1_AHR_Analysis](https://github.com/kristinmlee/dmc/blob/master/myAdidas/Chr1_AHR_Analysis) and [Chr10_ARNT_Analysis](https://github.com/kristinmlee/dmc/blob/master/myAdidas/Chr10_ARNT_Analysis). 

### Extensions to DMC

* Addition of the "Staggered Sweeps" Model: Additional migration model that allows the sweep to finish in the source population before the selected variant migrates to the recipient population. This is meant to be contrasted to the original migration model from Lee and Coop (2017), which assumes the sweeps occur at the same time. Here we call that original model "Concurrent Sweeps".
* Allowing for extremely strong selection regimes
* Allowing for differences in effective population size among populations
* For analysis of AHR on chromosome 1, we assigned the selected site to be the deletion midpoint instead of making it a free parameter
* For analysis of ARNT on chromosome 10, we modified the models to exclude selection in the source population of the putatively introgressed haplotype

### Steps to run DMC-MYAdIDAS

While there are some differences between our analyses of the AHR (on chromosome 1) and ARNT (on chromosome 10) regions, they follow the same general framework to calculate the composite likelihoods of each model. Below, we list scripts that must be run in the order they are listed; in each script, R objects are written that will be used as input for subsequent scripts. 

1. **calcNeutralF.R** contains code to generate and save **F** as an R object. This script is almost identical between the two analyses, so we leave it in this directory. 
2. **genSelMatrices_fxns.R** contains functions for each model to generate **F<sup>(S)</sup>** at a neutral site that is a certain recombination distance away from the selected site.
3. **genSelMatrices_exec.R** uses those functions in step 2 to generate **F<sup>(S)</sup>** for each recombination distance to the selected site and combination of free parameters. We define 1000 bins of recombination distances within which a neutral site can be categorized.
4. **calcInvDetSelMatrices.R** computes the inverses and determinants of all **F<sup>(S)</sup>** generated in step 3. These are used in the likelihood functions (see following steps).
5. **calcCompLike_fxns.R** contains functions to calculate the composite log-likelihood of each model. This script is identical between the two analyses, so we leave it in this directory.
6. **calcCompLike_exec.R** implements the functions defined in step 5. 
7. **genFigs.R** identifies and plots the profile composite log-likelihood surfaces of each model for each parameter. These plots were generated for the supplementary figures.

### Necessary input to run this method

The beginning of each script lists the R objects you must provide (within the code) to run it. Below is a list of all input that is required at some point when running the method. Estimates from each population must be listed in the same order for every R object that contains this information.

* matrix of allele frequencies at putatively neutral sites with	dimensions numberOfPopulations x numberOfSites.
* matrix of allele frequencies for region you are analyzing with dimension numberOfPopulations x numberOfSites. Exclude sites that are monomorphic across all populations.
* vector of genomic positions for region you are analyzing. Order of positions should match order that allele frequencies are listed (see above bullet)
* vector of sample sizes of length numberOfPopulations (i.e. number of chromosomes sampled in each population). 
* vector of effective population sizes of length numberOfPopulations
* per base pair recombination rate estimate for the region

* Parameter spaces for likelihood calculations. Not all models include all parameters.
    * selSite: vector of positions of proposed selected sites
    * sels: vector of proposed selection coefficients
    * intro_times: vector of proposed time in generations the variant is standing in populations after migration and before selection occurs
    * stdVar_times: vector of proposed time in generations the variant is standing in populations after selected populations split and before selection occurs
    * gs: vector of proposed frequencies of the standing variant
    * migs: vector of proposed migration rates (proportion of individuals from source each generation)
        * Note: cannot include 0

### Acknowledgments
We would like to thank Doc Edge for helping us come up with the name for this method, and [the music that inspired us](https://www.youtube.com/watch?v=JNua1lFDuDI).

