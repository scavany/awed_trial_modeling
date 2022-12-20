# wolbachia trial modeling

This repository contains code to reproduce results from the following pre-print:
SM Cavany<sup>&#42;</sup>, JH Huber<sup>&#42;</sup>, A Wieler, QM Tran, M Alkuzweny, M Elliott, G Espa&#241;a, SM Moore, TA Perkins (2022) **Ignoring transmission dynamics leads to underestimation of the impact of a novel intervention against mosquito-borne disease**. *medRxiv* doi:[10.1101/2021.11.19.21266602 ](https://doi.org/10.1101/2021.11.19.21266602)

### License

This code is released under the GNU General Public License v3.0

### Overview of R scripts
Scripts should be run in the following order:
* foi_calc.R estimates the force of infection in Yogyakarta based on seropositivity data from [Indriani et al](doi.org/10.4269/ajtmh.18-0315)
* fit_tempSIR.R fits the temperature dependent SIR model to the epidemiological data from Indriani et al
* hum_mvmt_data.R estimates the value of b based on the per-protocol analysis in [Utarini et al](doi.org/10.1056/NEJMoa2030243)
* run_Fig1_tempSIRvarC.R does the calculations underpinning figure 1
* fig_Fig1_tempSIRvarC.R creates figure 1
* run_Fig2_tempSIRvarC.R does the calculations underpinning figure 2
* fig_Fig2_tempSIRvarC.R creates figure 2
* fig_FigS7.R creates figure S7
* fig_FigS8.R creates figure S8
* fig_FigS9.R creates figure S9

The following two files contain necessary functions for the above scripts:
* functions_immune_history.R contains functions for calculating seropositivity from force of infection
* functions_trial_sim.R contains the ODE models and the functions for calculating the human movement parameters from the trial layout

### Overview of data
* AWED_Data.csv contains data from Utarini et al on the proportion of sampled mosquitoes that were wMel positive by cluster
* pop_by_age_Indonesia.csv contains the population age structure of Indonesia
* wei_movement.csv contains estimates of the efficacy from the per-protocol analysis of Utarini et al
* yogyakarta_dengue_incidence.csv and yogyakarta_seropositivity.csv contain data on dengue incidence and seropositivity in Yogyakarta respectively from Indriani et al