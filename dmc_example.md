
# DMC Example
**Kristin M Lee (10/23/2017)** 


This is an example of how to use the code found in
<https://github.com/kristinmlee/dmc> to calculate composite
log-likelihoods for various modes of convergent adaptation.

We simulate under the sceanario depicted in the picture below. The red
stars represent the three populations with selection. Ten alleles are
sampled from each of the six populations.

![Example demographic
sceanario](https://github.com/kristinmlee/dmc/blob/master/example/example_figures/example_tree.png)

Population allele frequencies and the positions of the SNPs are located in [example/selectedRegionAlleleFreqs_example.RDS](https://github.com/kristinmlee/dmc/blob/master/example/selectedRegionAlleleFreqs_example.RDS) and [example/selectedRegionPositions_example.RDS](https://github.com/kristinmlee/dmc/blob/master/example/selectedRegionPositions_example.RDS), respectively.

Calculate Neutral **F** Matrix
------------------------------

First, we calculate the neutral variance/covariance matrix (**F**).
[dmc/calcNeutralF.R](https://github.com/kristinmlee/dmc/blob/master/calcNeutralF.R)
contains code to generate and save **F** as an R object.

We need to specify (1) a matrix of allele frequencies with dimension
numberOfPopulations x numberOfSites, (2) a vector of sample sizes for
each population, (3) a string of filename/path we want our output to be
saved as. See
[dmc/calcNeutralF.R](https://github.com/kristinmlee/dmc/blob/master/calcNeutralF.R)
for more information about arguments.

10,000 neutral SNPs generated from simulating under the same demography and population structure with no selection are found here: [example/neutralAlleleFreqs_example.RDS](https://github.com/kristinmlee/dmc/blob/master/example/neutralAlleleFreqs_example.RDS).

    allFreqs = readRDS("example/neutralAlleleFreqs_example.RDS")
    sampleSizes = rep(10, 6)
    neutralF_filename = "example_output/neutralF_example"

    source("calcNeutralF.R")

Calculate **F<sup>(S)</sup>** Matrix
------------------------------------

Here, we generate matrices for the five following modes: 1) all selected
populations have independent mutations of beneficial allele 2) all
selected populations share beneficial allele via migration 3) the
beneficial allele was standing in the ancestor of all selected
populations 4) first two selected populations (1,3) share beneficial
allele via migration, third selected population (5) has an independent
mutation at the same locus  
5) the beneficial allele was standing in the ancestor of the first two
selected populations (1,3), third selected population (5) has an
independent mutation at the same locus

We will generate the mean-centered covariance matrices with selection,
**F<sup>(S)</sup>**, over our specified parameter spaces.
**F<sup>(S)</sup>** is a function of these parameters and the
recombination distance a site is away from a proposed selected site. We
generate **F<sup>(S)</sup>** for given bins of distance away from our
proposed putative selected sites. We first divide the distances in our
window into **numBins** bins and take the midpoint of the distances in
these bins to calculate **numBins** corresponding **F<sup>(S)</sup>**s
for a given set of parameters. See Appendix A.2.6 of Lee and Coop (2017)
for more information.

We'll start with models 1-3, where all selected populations share the
same mode of convergence.
[dmc/genSelMatrices\_individualModes.R](https://github.com/kristinmlee/dmc/blob/master/genSelMatrices_individualModes.R)
contains functions to generate **F<sup>(S)</sup>** for single modes of
convergent adaptation. See the R script for more information about
arguments.

    rec = 0.005 #per base pair recombination rate estimate for the region
    Ne = 10000
    numPops = 6
    sampleSizes = rep(10, 6)
    selPops = c(1, 3, 5)

    F_estimate = readRDS("example_output/neutralF_example.RDS")

    positions = readRDS("example/selectedRegionPositions_example.RDS")

    numBins = 1000

    selSite = seq(min(positions), max(positions), length.out = 10)
    sels = c(1e-4, 1e-3, 0.01, seq(0.02, 0.14, by = 0.01), seq(0.15, 0.3, by = 0.05), 
             seq(0.4, 0.6, by = 0.1))
    times = c(0, 5, 25, 50, 100, 500, 1000, 1e4, 1e6)
    gs = c(1/(2*Ne), 10^-(4:1))
    migs = c(10^-(seq(5, 1, by = -2)), 0.5, 1)
    sources = selPops

    source("genSelMatrices_individualModes.R")

In the following code, we generate and save nested lists of these
matrices over our specified parameter spaces. For a given parameter set,
we have a list of length **numBins** where each element contains the
corresponding **F<sup>(S)</sup>** for that bin.

For model 1 (all selected populations have independent mutations of
beneficial allele):

    FOmegas_ind = lapply(sels, function(sel) {
      calcFOmegas_indSweeps(sel)
    })

    saveRDS(FOmegas_ind, "example_output/FOmegas_ind_example.RDS")

For model 2 (all selected populations share beneficial allele via
migration):

    FOmegas_mig = lapply(sels ,function(sel) {
      lapply(migs, function(mig) {
        lapply(sources, function(my.source) {
          calcFOmegas_mig(sel, mig, my.source)
        })
      })
    })

    saveRDS(FOmegas_mig, "example_output/FOmegas_mig_example.RDS")

For model 3 (the beneficial allele was standing in the ancestor of all
selected populations):

    FOmegas_sv = lapply(sels, function(sel) {
      lapply(gs, function(g) {
        lapply(times, function(time) {
          lapply(sources, function(my.source) {
            calcFOmegas_stdVar.source(sel, g, time, my.source)
          })
        })
      })
    })

    saveRDS(FOmegas_sv, "example_output/FOmegas_sv_example.RDS")

Now we'll work with models 4 and 5, where there are multiple modes of
convergent adaptation (i.e. different sets of populations are aquiring
the beneficial allele at the same site via different modes).
[dmc/genSelMatrices\_multipleModes.R](https://github.com/kristinmlee/dmc/blob/master/genSelMatrices_multipleModes.R)
contains functions to generate **F<sup>(S)</sup>** for multiple modes of
convergence. See the R script for more information about arguments.

    rec = 0.005 #per base pair recombination rate estimate for the region
    Ne = 10000
    numPops = 6
    sampleSizes = rep(10, 6)
    selPops = c(1, 3, 5)

    sets = list(c(1, 3), 5)
    #populations 1 and 3 will share a mode, population 5 has separate mode

    F_estimate = readRDS("example_output/neutralF_example.RDS")

    positions = readRDS("example/selectedRegionPositions_example.RDS")

    numBins = 1000

    selSite = seq(min(positions), max(positions), length.out = 10)
    sels = c(1e-4, 1e-3, 0.01, seq(0.02, 0.14, by = 0.01), seq(0.15, 0.3, by = 0.05), 
             seq(0.4, 0.6, by = 0.1))
    times = c(0, 5, 25, 50, 100, 500, 1000, 1e4, 1e6)
    gs = c(1/(2*Ne), 10^-(4:1))
    migs = c(10^-(seq(5, 1, by = -2)), 0.5, 1)
    sources = selPops

    source("genSelMatrices_multipleModes.R")

For model 4 (first two selected populations (1,3) share beneficial
allele via migration, third selected population (5) has an independent
mutation at the same locus):

    #we first define vector of modes corresponding the populations specified in sets
    my.modes_migInd = c("mig", "ind")

    #the parameters time and g are not involved in the migration model so we only loop over
    ## the first element of these vectors
    FOmegas_mixed_migInd = lapply(sels ,function(sel) {
        lapply(gs[1], function(g) {
            lapply(times[1], function(time) {
                lapply(migs, function(mig) {
                    lapply(sources, function(my.source) {
                        calcFOmegas_mixed(sel, g, time, mig, my.source, my.modes_migInd)
                    })
                })
            })
        })
    })

    saveRDS(FOmegas_mixed_migInd, "example_output/FOmegas_mixed_migInd_example.RDS")

For model 5 (the beneficial allele was standing in the ancestor of the
first two selected populations (1,3), third selected population (5) has
an independent mutation at the same locus):

    #we first define vector of modes corresponding the populations specified in sets
    my.modes_svInd = c("sv", "ind")

    #the parameter mig is not involved in the standing variant model so we only loop over
    ## the first element of this vector
    FOmegas_mixed_svInd = lapply(sels ,function(sel) {
        lapply(gs, function(g) {
            lapply(times, function(time) {
                lapply(migs[1], function(mig) {
                    lapply(sources, function(my.source) {
                        calcFOmegas_mixed(sel, g, time, mig, my.source, my.modes_svInd)
                    })
                })
            })
        })
    })

    saveRDS(FOmegas_mixed_svInd, "example_output/FOmegas_mixed_svInd_example.RDS")

Generate and save inverses and determinants for **F<sup>(S)</sup>** matrices
----------------------------------------------------------------------------

To avoid the costly step of recomputing the corresponding inverses and
determinants of the covariance matrices needed for selection used in
likelihood calculations, we do this step first and use these values for
all SNPs in a given bin, and store them and reuse them over all
locations of the selected site. See Appendix A.2.6 of Lee and Coop
(2017) for more information.

[dmc/calcInvDetSelMatrices\_all.R](https://github.com/kristinmlee/dmc/blob/master/calcInvDetSelMatrices_all.R)
contains scripts to generate all the possible inverses and determinants.
We repeat the relevant code for our specified models here to both
compute and save these values.

    library("MASS")

    ## Neutral model
    sampleSizes = rep(10, 6)
    F_estimate = readRDS("example_output/neutralF_example.RDS")

    numPops = 6
    M = numPops
    Tmatrix = matrix(data = rep(-1 / M, (M - 1) * M), nrow = M - 1, ncol = M)
    diag(Tmatrix) = (M - 1) / M 

    sampleErrorMatrix = diag(1/sampleSizes, nrow = numPops, ncol = numPops)

    det_FOmegas_neutral = det(Tmatrix %*% (F_estimate + sampleErrorMatrix) %*% t(Tmatrix))
    saveRDS(det_FOmegas_neutral, "example_output/det_FOmegas_neutral_example.RDS")

    inv_FOmegas_neutral = ginv(Tmatrix %*% (F_estimate + sampleErrorMatrix) %*% t(Tmatrix))
    saveRDS(inv_FOmegas_neutral, "example_output/inv_FOmegas_neutral_example.RDS")



    ## Model 1
    FOmegas_ind = readRDS("example_output/FOmegas_ind_example.RDS")

    det_FOmegas_ind = lapply(FOmegas_ind, function(sel) {
        lapply(sel, function(dist) {
            det(dist)
        })
    })
    saveRDS(det_FOmegas_ind, "example_output/det_FOmegas_ind_example.RDS")

    inv_FOmegas_ind = lapply(FOmegas_ind, function(sel) {
        lapply(sel, function(dist) {
            ginv(dist)
        })
    })
    saveRDS(inv_FOmegas_ind, "example_output/inv_FOmegas_ind_example.RDS")



    ## Model 2
    FOmegas_mig = readRDS("example_output/FOmegas_mig_example.RDS")

    det_FOmegas_mig = lapply(FOmegas_mig, function(sel) {
        lapply(sel, function(mig) {
            lapply(mig, function(source) {
                lapply(source, function(dist) {
                    det(dist)
                })
            })
        })
    })
    saveRDS(det_FOmegas_mig_output, "example_output/det_FOmegas_mig_example.RDS")

    inv_FOmegas_mig = lapply(FOmegas_mig, function(sel) {
        lapply(sel, function(mig) {
            lapply(mig, function(source) {
                lapply(source, function(dist) {
                    ginv(dist)
                })
            })
        })
    })
    saveRDS(inv_FOmegas_mig_output, "example_output/inv_FOmegas_mig_example.RDS")



    ## Model 3
    FOmegas_sv = readRDS("example_output/FOmegas_sv_example.RDS")

    det_FOmegas_sv = lapply(FOmegas_sv, function(sel) {
        lapply(sel, function(g) {
            lapply(g, function(time) {
                lapply(time, function(my.source) {
                    lapply(my.source, function(dist) {
                        det(dist)
                    })
                })
            })
        })
    })
    saveRDS(det_FOmegas_sv, "example_output/det_FOmegas_sv_example.RDS")

    inv_FOmegas_sv = lapply(FOmegas_sv, function(sel) {
        lapply(sel, function(g) {
            lapply(g, function(time) {
                lapply(time, function(my.source) {
                    lapply(my.source, function(dist) {
                        ginv(dist)
                    })
                })
            })
        })
    })
    saveRDS(inv_FOmegas_sv, "example_output/inv_FOmegas_sv_example.RDS")

    ## Model 4
    FOmegas_mixed_migInd = readRDS("example_output/FOmegas_mixed_migInd_example.RDS")

    detFOmegas_mixed_migInd = lapply(FOmegas_mixed_migInd, function(sel) {
        lapply(sel, function(g) {
            lapply(g, function(time) {
                lapply(time, function(mig) {
                    lapply(mig, function(source) {
                        lapply(source, function(dist) {
                            det(dist)
                        })
                    })  
                })
            })
        })
    })
    saveRDS(detFOmegas_mixed_migInd, "example_output/det_FOmegas_mixed_migInd_example.RDS")

    invFOmegas_mixed_migInd = lapply(FOmegas_mixed_migInd, function(sel) {
        lapply(sel, function(g) {
            lapply(g, function(time) {
                lapply(time, function(mig) {
                    lapply(mig, function(source) {
                        lapply(source, function(dist) {
                            ginv(dist)
                        })
                    })  
                })
            })
        })
    })
    saveRDS(invFOmegas_mixed_migInd, "example_output/inv_FOmegas_mixed_migInd_example.RDS")

    ## Model 5
    FOmegas_mixed_svInd = readRDS("example_output/FOmegas_mixed_svInd_example.RDS")

    detFOmegas_mixed_svInd = lapply(FOmegas_mixed_svInd, function(sel) {
        lapply(sel, function(g) {
            lapply(g, function(time) {
                lapply(time, function(mig) {
                    lapply(mig, function(source) {
                        lapply(source, function(dist) {
                            det(dist)
                        })
                    })  
                })
            })
        })
    })
    saveRDS(detFOmegas_mixed_svInd, "example_output/det_FOmegas_mixed_svInd_example.RDS")

    invFOmegas_mixed_svInd = lapply(FOmegas_mixed_svInd, function(sel) {
        lapply(sel, function(g) {
            lapply(g, function(time) {
                lapply(time, function(mig) {
                    lapply(mig, function(source) {
                        lapply(source, function(dist) {
                            ginv(dist)
                        })
                    })  
                })
            })
        })
    })
    saveRDS(invFOmegas_mixed_svInd, "example_output/inv_FOmegas_mixed_svInd_example.RDS")

Calculate composite likelihoods
-------------------------------

We use the inverses and determinants of the covariance matrices
generated above to calculate the log-likelihood of a site a given
distance away from a proposed selected site. We sum the log-likelihoods
across all SNPs in the window to obtain a composite log-likelihood under
a given convergence model with a set of parameters for a proposed
selected site.
[dmc/calcCompositeLike.R](https://github.com/kristinmlee/dmc/blob/master/calcCompositeLike.R)
contains functions for calculating the composite log-likelihoods for all
models.

Below we calculate and save the composite log-likelihoods for all five
models (and the neutral model). Composite likelihoods are stored as
nested lists for the parameters specified (see
[dmc/calcCompositeLike.R](https://github.com/kristinmlee/dmc/blob/master/calcCompositeLike.R)
and below for more information).

First, we need to randomize with respect the reference allele. We will
save these allele frequencies since they should be used for all
likelihood calculations.

    freqs_notRand = readRDS("example/selectedRegionAlleleFreqs_example.RDS")

    randFreqs = apply(freqs_notRand, 2, function(my.freqs) {
        if(runif(1) < 0.5) {
            my.freqs = 1 - my.freqs
        }
        my.freqs
    })

    saveRDS(randFreqs, "example_output/selectedRegionAlleleFreqsRand_example.RDS")

    numPops = 6
    positions = readRDS("example/selectedRegionPositions_example.RDS")

    freqs = readRDS("example_output/selectedRegionAlleleFreqsRand_example.RDS")

    #these values must be same as used to calculate Calculate F^(S) matrices above
    numBins = 1000
    selSite = seq(min(positions), max(positions), length.out = 10)
    sels = c(1e-4, 1e-3, 0.01, seq(0.02, 0.14, by = 0.01), seq(0.15, 0.3, by = 0.05), 
             seq(0.4, 0.6, by = 0.1))
    times = c(0, 5, 25, 50, 100, 500, 1000, 1e4, 1e6)
    Ne = 10000 #only necessary to specify since we define lowest value of g by Ne
    gs = c(1/(2*Ne), 10^-(4:1))
    migs = c(10^-(seq(5, 1, by = -2)), 0.5, 1)
    selPops = c(1, 3, 5) 
    #only necessary to specify since we define possible source populations as all
    ## selected populations
    sources = selPops

    source("calcCompositeLike.R")

    ## Neutral model
    det_FOmegas_neutral = readRDS("example_output/det_FOmegas_neutral_example.RDS")
    inv_FOmegas_neutral = readRDS("example_output/inv_FOmegas_neutral_example.RDS")
    compLikelihood_neutral = lapply(1 : length(selSite), function(j) {
      calcCompLikelihood_neutral(j, det_FOmegas_neutral, inv_FOmegas_neutral))
    }
    saveRDS(compLikelihood_neutral, "example_output/compLikelihood_neutral_example.RDS")

    ## Model 1
    det_FOmegas_ind = readRDS("example_output/det_FOmegas_ind_example.RDS")
    inv_FOmegas_ind = readRDS("example_output/inv_FOmegas_ind_example.RDS")
    compLikelihood_ind = lapply(1 : length(selSite), function(j) {
      lapply(1 : length(sels), function(sel) calcCompLikelihood_1par(j, det_FOmegas_ind,
                                                                     inv_FOmegas_ind, sel))
    })
    saveRDS(compLikelihood_ind, "example_output/compLikelihood_ind_example.RDS")

    ## Model 2
    det_FOmegas_mig = readRDS("example_output/det_FOmegas_mig_example.RDS")
    inv_FOmegas_mig = readRDS("example_output/inv_FOmegas_mig_example.RDS")
    compLikelihood_mig = lapply(1 : length(selSite), function(j) {
        lapply(1 : length(sels), function(sel) {
            lapply(1 : length(migs), function(mig) {
                lapply(1 : length(sources), function(my.source) {
                    calcCompLikelihood_3par(j, det_FOmegas_mig, inv_FOmegas_mig, sel, mig,
                                            my.source)
                })
            })
        })
    })
    saveRDS(compLikelihood_mig, "example_output/compLikelihood_mig_example.RDS")

    ## Model 3
    det_FOmegas_sv = readRDS("example_output/det_FOmegas_sv_example.RDS")
    inv_FOmegas_sv = readRDS("example_output/inv_FOmegas_sv_example.RDS")
    compLikelihood_sv = lapply(1 : length(selSite), function(j) {
        lapply(1 : length(sels), function(sel) {
            lapply(1 : length(gs), function(g) {
                lapply(1 : length(times), function(t) {
                    lapply(1: length(sources), function(my.source) {
                    calcCompLikelihood_4par(j, det_FOmegas_sv, inv_FOmegas_sv, sel, g, t,
                                            my.source)
                    })
                })
            })
        })
    })
    saveRDS(compLikelihood_sv, "example_output/compLikelihood_sv_example.RDS")

    ## Model 4
    det_FOmegas_mixed_migInd = readRDS("example_output/det_FOmegas_mixed_migInd_example.RDS")
    inv_FOmegas_mixed_migInd = readRDS("example_output/inv_FOmegas_mixed_migInd_example.RDS")

    # same trick as above (the parameters time and g are not involved in the migration
    ## model so we only loop over the first element of these vectors)
    # now save lists for each proposed selected site (may want to do this for other 
    ## models/more elegantly depending on density of parameter space)
    for(j in 1 : length(selSite)) {
        compLikelihood_mixed_migInd = lapply(1 : length(sels), function(sel) {
            lapply(1 : length(gs[1]), function(g) {
                lapply(1 : length(times[1]), function(t) {
                  lapply(1 : length(migs), function (mig) {
                      lapply(1: length(sources), function(my.source) {
                        calcCompLikelihood_5par(j, det_FOmegas_mixed_migInd,
                                                inv_FOmegas_mixed_migInd, sel, g, t, mig,
                                                my.source)
                      })
                    })
                })
            })
        })
      saveRDS(compLikelihood_mixed_migInd,
              paste("example_output/compLikelihood_mixed_migInd_example_selSite", j, ".RDS",
                    sep = ""))
    }

    ## Model 5
    det_FOmegas_mixed_svInd = readRDS("example_output/det_FOmegas_mixed_svInd_example.RDS")
    inv_FOmegas_mixed_svInd = readRDS("example_output/inv_FOmegas_mixed_svInd_example.RDS")

    #same trick as above (the parameter mig is not involved in the migration model so we
    ##only loop over the first element of this vector)
    # now save lists for each proposed selected site (may want to do this for other
    ## models/more elegantly depending on density of parameter space)
    for(j in 1 : length(selSite)) {
        compLikelihood_mixed_svInd = lapply(1 : length(sels), function(sel) {
            lapply(1 : length(gs), function(g) {
                lapply(1 : length(times), function(t) {
                  lapply(1 : length(migs[1]), function (mig) {
                      lapply(1: length(sources), function(my.source) {
                        calcCompLikelihood_5par(j, det_FOmegas_mixed_svInd,
                                                inv_FOmegas_mixed_svInd, sel, g, t, mig,
                                                my.source)
                      })
                    })
                })
            })
        })
      saveRDS(compLikelihood_mixed_svInd,
              paste("example_output/compLikelihood_mixed_svInd_example_selSite", j, ".RDS",
                    sep = ""))
    }

Plotting and obtaining maximum composite likelihood parameters and interpreting results
---------------------------------------------------------------------------------------

First, let's combine files for models 4 and 5 that are separated by
proposed selected site.

    positions = readRDS("example/selectedRegionPositions_example.RDS")
    selSite = seq(min(positions), max(positions), length.out = 10)

    ## Model 4
    compLikelihood_mixed_migInd_all = lapply(1: length(selSite), function(i) {
      readRDS(paste("example_output/compLikelihood_mixed_migInd_example_selSite", i, ".RDS",
                    sep = ""))
    })

    saveRDS(compLikelihood_mixed_migInd_all, "example_output/compLikelihood_mixed_migInd_example.RDS")

    ## Model 5
    compLikelihood_mixed_svInd_all = lapply(1: length(selSite), function(i) {
      readRDS(paste("example_output/compLikelihood_mixed_svInd_example_selSite", i, ".RDS",
                    sep = ""))
    })

    saveRDS(compLikelihood_mixed_svInd_all, "example_output/compLikelihood_mixed_svInd_example.RDS")

### Plot maximum composite likelihood ratio (model - neutral) for all models over proposed selected sites

    positions = readRDS("example_output/selectedRegionPositions_example.RDS")
    selSite = seq(min(positions), max(positions), length.out = 10)

    #read in composite likelihood files and calculate max for all proposed selected sites
    compLikelihood_neutral = readRDS("example_output/compLikelihood_neutral_example.RDS")
    compLikelihood_neutral_site = sapply(1 : length(selSite), function(i) {
      max(unlist(compLikelihood_neutral[[i]]))
    })

    compLikelihood_ind = readRDS("example_output/compLikelihood_ind_example.RDS")
    compLikelihood_ind_site = sapply(1 : length(selSite), function(i) {
      max(unlist(compLikelihood_ind[[i]]))
    })

    compLikelihood_mig = readRDS("example_output/compLikelihood_mig_example.RDS")
    compLikelihood_mig_site = sapply(1 : length(selSite), function(i) {
      max(unlist(compLikelihood_mig[[i]]))
    })

    compLikelihood_sv = readRDS("example_output/compLikelihood_sv_example.RDS")
    compLikelihood_sv_site = sapply(1 : length(selSite), function(i) {  
      max(unlist(compLikelihood_sv[[i]]))
    })

    compLikelihood_mixed_migInd = readRDS("example_output/compLikelihood_mixed_migInd_example.RDS")
    compLikelihood_mixed_migInd_site = sapply(1 : length(selSite), function(i) {
      max(unlist(compLikelihood_mixed_migInd[[i]]))
    })

    compLikelihood_mixed_svInd = readRDS("example_output/compLikelihood_mixed_svInd_example.RDS")
    compLikelihood_mixed_svInd_site = sapply(1 : length(selSite), function(i) {
      max(unlist(compLikelihood_mixed_svInd[[i]]))
    })

    plot_range = range(c((compLikelihood_ind_site - compLikelihood_neutral_site),
                         (compLikelihood_mig_site - compLikelihood_neutral_site), 
                         (compLikelihood_sv_site - compLikelihood_neutral_site), 
                         (compLikelihood_mixed_migInd_site - compLikelihood_neutral_site),
                         (compLikelihood_mixed_svInd_site - compLikelihood_neutral_site)))

    plot(selSite, compLikelihood_ind_site - compLikelihood_neutral_site, type = "b",
         ylim = c(plot_range[1] - 50, plot_range[2] + 50),
         xlab = "Proposed position selected site",
         ylab = "Composite log-likelihood (model - neutral)")
    lines(selSite, compLikelihood_mig_site - compLikelihood_neutral_site, col = "red",
          type = "b")
    lines(selSite, compLikelihood_sv_site - compLikelihood_neutral_site, col = "blue",
          type = "b")
    lines(selSite, compLikelihood_mixed_migInd_site - compLikelihood_neutral_site,
          col = "orange", lty = 2, type = "b")
    lines(selSite, compLikelihood_mixed_svInd_site - compLikelihood_neutral_site,
          col = "green", lty = 2, type = "b")
    legend("topright", col = c("black", "red", "blue", "orange", "green"),
           lty = c(rep(1, 3), rep(2, 2)), sapply(1 : 5, function(i) paste("Model", i)),
           cex = 0.5)

![Figure 1](https://github.com/kristinmlee/dmc/blob/master/example/example_figures/figure1.png)

### Get maximum composite likelihood estimates (MCLEs)

[dmc/getMCLE.R](https://github.com/kristinmlee/dmc/blob/master/getMCLE.R)
contains functions to get the maximum composite likelihood estimates of
parameters for all convergence models.

    source("getMCLE.R")

    sels = c(1e-4, 1e-3, 0.01, seq(0.02, 0.14, by = 0.01), seq(0.15, 0.3, by = 0.05),
             seq(0.4, 0.6, by = 0.1))
    times = c(0, 5, 25, 50, 100, 500, 1000, 1e4, 1e6)
    Ne = 10000 #only necessary to specify since we define lowest value of g by Ne
    gs = c(1/(2*Ne), 10^-(4:1))
    migs = c(10^-(seq(5, 1, by = -2)), 0.5, 1)
    selPops = c(1, 3, 5)
    #only necessary to specify since we define possible source populations as all
    ## selected populations
    sources = selPops

    ## Model 1
    getMCLEind(compLikelihood_ind, selSite, sels)

    ##      maxLoc maxSel
    ## [1,] 0.0017   0.03

    ## Model 2
    getMCLEmig(compLikelihood_mig, selSite, sels, migs, sources)

    ##      maxLoc maxSel maxMig maxSource
    ## [1,] 0.0017   0.03  1e-05         3

    ## Model 3
    getMCLEsv_source(compLikelihood_sv, selSite, sels, gs, times, sources)

    ##      maxLoc maxSel maxG maxTime maxSource
    ## [1,] 0.0017   0.05 0.01   10000         5

    ## Model 4
    getMCLEmixed(compLikelihood_mixed_migInd, selSite, sels, gs[1], times[1], migs, sources)

    ##      maxLoc maxSel  maxG maxTime maxSource
    ## [1,] 0.0017   0.03 5e-05       0         1

    ## Model 5
    getMCLEmixed(compLikelihood_mixed_svInd, selSite, sels, gs, times, migs[1], sources)

    ##      maxLoc maxSel maxG maxTime maxSource
    ## [1,] 0.0017   0.03 0.01   10000         5

### Plotting some profile likelihood surfaces for parameters

In Figure 1, we see that the models 3 and 5 (where selection is on
standing variation present in all or a subset of populations) have the
highest composite log-likelihoods. For these models, we plot the profile
likelihood surfaces for the parameter of time the beneficial allele has
been standing independently in the selected populations, at the most
likely selected site. See Lee and Coop (2017) for more information.

    ## Model 3
    mcle_sv = getMCLEsv_source(compLikelihood_sv, selSite, sels, gs, times, sources)

    compLike_sv_byTime = lapply(1 : length(times), function(time) {
      sapply(1: length(sels), function(sel) {
        sapply(1 : length(gs), function(g) {
          compLikelihood_sv[[mcle_sv[1]]][[sel]][[g]][[time]]
        })
      })
    })

    profileLike_time_sv = sapply(1: length(compLike_sv_byTime), function(i) {
      max(unlist(compLike_sv_byTime[[i]]))
    })

    plot(times, profileLike_time_sv, type = "b", xlab = "Time",
         ylab = "Profile composite log-likelihood", main = "Figure 2: Model 3")
    abline(v = mcle_sv[4], lty = 2, col = "red")

![Figure 2](https://github.com/kristinmlee/dmc/blob/master/example/example_figures/figure2.png)

    ## Model 5
    mcle_mixed_svInd = getMCLEmixed(compLikelihood_mixed_svInd, selSite, sels,
                                    gs, times, migs[1], sources)

    compLike_mixed_svInd_byTime = lapply(1 : length(times), function(time) {
      sapply(1: length(sels), function(sel) {
        sapply(1 : length(gs), function(g) {
          compLikelihood_mixed_svInd[[mcle_mixed_svInd[1]]][[sel]][[g]][[time]]
        })
      })
    })

    profileLike_time_mixed_svInd = sapply(1: length(compLike_mixed_svInd_byTime), function(i) {
      max(unlist(compLike_mixed_svInd_byTime[[i]]))
    })

    plot(times, profileLike_time_mixed_svInd, type = "b", xlab = "Time",
         ylab = "Profile composite log-likelihood", main = "Figure 3: Model 5")
    abline(v = mcle_mixed_svInd[4], lty = 2, col = "red")

![Figure 3](https://github.com/kristinmlee/dmc/blob/master/example/example_figures/figure3.png)

In both cases, we see the MCLE of t tending towards infinity. In this
case, we know our standing variant model overlaps with our independent
mutations model (model 1). The slight increase in composite
log-likelihood values between these models and model 1 most likely
reflects the fact that models 3 and 5 are more parameterized.
Additionally, we see little difference in the composite log-likelihoods
of models 3 and 5 where three and two populations, respectively, are
sharing the sweep via selection and standing variation and the third
population has an independent mutation at the same locus. This provides
further evidence that independent mutations (or selection very old
standing variation) generated the patterns observed in the data. The
point estimates of MCLE of t are marked by red lines in Figures 2 and 3.

##### The true data was simulated with independent mutations in all three selected populations.

Additionally, we are interested in the MCLE for the strength of
selection. We see in our best fitting models (models 1, 3, and 5) that
the MCLE of s is at 0.05. We plot profile likelihood surfaces of s at
the most likely selected site for each of these models.

    ## Model 1
    mcle_ind = getMCLEind(compLikelihood_ind, selSite, sels)

    profileLike_sel_ind = sapply(1: length(sels), function(i) {
      max(unlist(compLikelihood_ind[[mcle_ind[1]]][[i]]))
    })

    ## Model 3
    profileLike_sel_sv = sapply(1: length(sels), function(i) {
      max(unlist(compLikelihood_sv[[mcle_sv[1]]][[i]]))
    })

    ## Model 5
    profileLike_sel_mixed_svInd = sapply(1: length(sels), function(i) {
      max(unlist(compLikelihood_mixed_svInd[[mcle_mixed_svInd[1]]][[i]]))
    })

    par(mfrow = c(1, 3))
    plot(sels, profileLike_sel_sv, type = "b", xlab = "Sel",
         ylab = "Profile composite log-likelihood", main = "Model 3")
    abline(v = mcle_sv[2], lty = 2, col = "blue")

    plot(sels, profileLike_sel_mixed_svInd, type = "b", xlab = "Sel",
         ylab = "Profile composite log-likelihood", main = "Model 5")
    abline(v = mcle_mixed_svInd[2], lty = 2, col = "blue")

    plot(sels, profileLike_sel_ind, type = "b", xlab = "Sel",
         ylab = "Profile composite log-likelihood", main = "Model 1")
    abline(v = mcle_ind[2], lty = 2, col = "blue")

![Figure 4](https://github.com/kristinmlee/dmc/blob/master/example/example_figures/figure4.png)

Models 1 and 5 show a peak at s = 0.03 and model 3 at s = 0.05.

##### The true data was simulated with s = 0.05.

Lastly, we look at the location of the selected site. In Figure 1 and in
the MCLE values, we see the highest likelihoods are at the leftmost end
of the scaffold (position 0.0017, the leftmost proposed selected site).

##### The true data was simulated with the selection site at position 0.
