Scripts for analysis of experimental data from PdR1 preference and transmission experiment.

Directory "R_functions" contains scripts of functions used repeatedly for statistics and epidemic models
Directory "R_analyses" contains scripts for statistical analysis and running epidemic models
Directory "output" contains cleaned and organized data files (in .rds format)
Directory "data" contains some raw data files (in .xlsx format)

Within R_analyses directory:
1) "preference_consumer_movement_model.R" contains code for estimating experimental vector insect preference using Consumer Movement Model.
2) "pdr1_transmission_analysis.R" contains code for analyses of PD symptom severity, Xylella fastidiosa population size, vector acquisition, and vector transmission data. Includes non-linear model fitting and assessment for acquisition and transmission data.
3) "pdr1_sem_transmission_analysis.R" contains code for structural equation modeling (SEM) of 2017 experimental data.
4) "pdr1_reinfection_analysis.R" contains code for analysis of reinfection experiment data reported in Wallis et al. paper on phytochemistry
5) "cmm_assumptions_tests.R" contains code for testing of assumtions of the Consumer Movement Model.
6) "pdr1_synthesis_analyses.R" contains code for elastic net algorithm fitting of plant phenolic data.
