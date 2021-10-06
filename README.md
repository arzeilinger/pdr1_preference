Scripts for analyses for publication Zeilinger et al. (2021) “Plant defense against a pathogen drives nonlinear transmission....” Ecosphere.

Data are from an experiment testing the effect of a disease-resistant grape cultivar (PdR1) on feeding preference by insect vectors, Blue-green sharpshooter (or BGSS), and on transmission of the pathogen Xylella fastidiosa (Xf). Experiment was repeated over two years, 2016 and 2017.

Vector preference between two choices plants (Xf-infected vs. Xf-free) was assessed by recording the number of insects on each plant choice over time in trials. These data were then analyzed in the Consumer Movement Model (CMM) of Zeilinger et al. (2014) to estimate attraction and leaving rates of insects from each plant. The model fitting can be found in the script “R_analyses/preference_consumer_movement_model.R”. Additionally, analyses testing the assumptions of the CMM can be found in “R_analyses/cmm_assumptions_tests.R”.

Analyses of plant infection and transmission data can be found in the script “R_analyses/pdr1_transmission_analysis.R”. Specifically, the analyses in this script include:

1. Data on symptom severity of Xf infection were collected as a 1-5 scale (see Zeilinger et al. (2021) for details) and differences between resistant and susceptible grapevines were tested with a partial-odds logistic model because of the ordinal scale of the data. 

2. Data on pathogen population size in infected plants were collected using bacterial culturing and expressed as colony forming units (CFU)/g plant tissue. Population size data were analyzed using a quasi-Poisson generalized linear model because the population size data were Poisson distributed and over-dispersed. Preliminary analyses suggested that the quasi-Poisson model was superior to log10 transformation of the response and a negative Binomial model (not shown).

3. Data on infection of insect vectors (vector acquisition data) were collected using qPCR and are expressed as a binary variable of infectious (1) or non-infectious (0). These data were to a series of non-linear models relating how infectiousness of vectors changed over time. The models included the Logistic Growth, Holling Type II, Ricker, and Holling Type IV models. The first two are saturating models while the last two are uni-modal or “humped”. The Logistic Growth and Ricker belong to the exponential family of models while the Holling Type II and Type IV belong to the rational family of models. See Appendix S3 of Zeilinger et al. (2021) for more details on the models. I fit the models using maximum likelihood estimation and selected the best model using AIC.

4. Data on infection of test plants (initially uninfected or “Xf-free”) was collected using bacterial culturing and again expressed as a binary variable. For analysis, I used the same scheme of fitting non-linear models to the data using MLE as described above for vector acquisition data.

Finally, in "R_analyses/pdr1_sem_analysis.R", I tested which variables described above best explains variation in transmission to the test plants using structural equation modeling (SEM). I chose SEM because it allowed me to explore multiple causal pathways that lead to transmission. I tested multiple models based on our hypotheses on the transmission process and selected the best model using AIC and Fisher’s C statistic. 


[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4547775.svg)](https://doi.org/10.5281/zenodo.4547775)
