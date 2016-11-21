#### ANALYSIS OF PDR1 TRANSMISSION DATA

my.packages <- c("tidyr", "dplyr", "data.table", "lme4", "lmerTest", "xlsx")
lapply(my.packages, require, character.only = TRUE)


##############################################################################################################
#### Combining data sets

# Getting symptoms index
prefdata <- read.xlsx("data/pdr1_preference_data.xlsx", sheetName = "data")
str(prefdata)

# Separate leaf and symptom data from preference count data
leafdata <- prefdata[,c("week", "trt", "rep", "n_leaves_test", "n_leaves_source_start", 
                        "n_leaves_source_end", "n_pd_leaves", "n_ms_petioles", "pd_index")] %>% dplyr::filter(., !is.na(n_pd_leaves))
leafdata$n_leaves_test <- factor2numeric(leafdata$n_leaves_test)
leafdata$n_leaves_source_start <- factor2numeric(leafdata$n_leaves_source_start)
leafdata$n_leaves_source_end <- factor2numeric(leafdata$n_leaves_source_end)
# Some trials I only counted leaves at beginning, other trials only at end; combine total leaf data
leafdata$n_leaves_source <- leafdata %>% with(., cbind(n_leaves_source_start, n_leaves_source_end)) %>% rowMeans(., na.rm = TRUE)
leafdata$n_leaves_source[is.nan(leafdata$n_leaves_source)] <- NA
str(leafdata)

#### NOTE: Could convert my number of PD leaves and matchstick petioles to the 0-5 scale of Rashed et al. 2013, I think. Look into it.

################################################################################################################
#### Constructing PD indices

#### PD index of Rashed et al. 2013:
# 0 = asymptomatic
# 1 = 1-2 scorched leaves
# 2 = 3-4 scorched leaves
# 3 = all leaves scorched; a few matchstick petioles
# 4 = all leaves heavily scorched and many matchstick petioles
# 5 = only a few leaves at the end of the cane

table(leafdata$n_pd_leaves, leafdata$n_ms_petioles)

#### Alternative PD index: 
leafdata$pd_index2 <- leafdata %>% with(., (n_pd_leaves/n_leaves_source) + n_ms_petioles)
plot(x = leafdata$pd_index, y = leafdata$pd_index2)


################################################################################################################
#### Constructing vector infection index
# In each trial, I have Xf pops for each of the vectors. I could include in transmission model:
# total Xf population among all vectors
# proportion of vectors infected (b/c some cages have fewer than 8 vectors)
# an evenness index of Xf pops amongs vectors
# Matt's transmission parameters paper might have something to say about this
# multiple possibilities can be evaluated using AIC


# draft model
# test.infection.status ~ genotype*week + xf.pop.source + pd_index(1 or 2) + "vector.infection.index" + p1 + p2 + mu1 + mu2