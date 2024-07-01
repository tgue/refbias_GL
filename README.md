# Estimating ancestry proportions and genotype likelihoods in the presence of mapping bias

This repository contains scripts used to conduct simulations of realistic aDNA sequence data from several populations under a demographic model. The data is then used as input for different ancestry estimation tools. Furthermore, we provide scripts to calculate a modified genoype likelihood for ascertained biallelic SNPs accounting for mapping bias.

## Calculating a corrected genotype likelihoods


## Simulations

The folder `simulations` contains scripts for running the simulations for the study, using `msprime` and `gargammel` to generate the data and then `ADMIXTURE`, `qpAdm`, `ngsadmix` and `fastngsadmix` to estimate ancestry proportions. The simulations are started from the script `msp_gg_loop_fngsa_correctedGL.sh`.
