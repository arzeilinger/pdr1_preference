#### Analysis of PdR1 re-infection experiment 2017

rm(list = ls())
# Load packages
my.packages <- c("tidyr", "dplyr", "data.table", "openxlsx", "ggplot2", "DHARMa",
                 "MASS", "logistf", "multcomp", "bbmle", "lme4", "googlesheets", "broom")
lapply(my.packages, require, character.only = TRUE)

source("R_functions/factor2numeric.R")
source("R_functions/standardize.R")

#### Import culturing data
# reinfdata <- read.xlsx("data/2017_data/PdR1_2017_preference-transmission_experiment_data.xlsx", sheet = "Re-infection_culturing_data")
# Import from Googlesheets
pdr1DataURL <- gs_url("https://docs.google.com/spreadsheets/d/14uJLfRL6mPrdf4qABeGeip5ZkryXmMKkan3mJHeK13k/edit?usp=sharing",
                      visibility = "private")
reinfdataGS <- gs_read(pdr1DataURL, ws = "Re-infection_culturing_data")


reinfdata <- reinfdataGS
str(reinfdata)
summary(reinfdata)
# Make xf_cfu_per_g a numeric variable and genotype and trt factors
reinfdata$xf_cfu_per_g <- as.numeric(reinfdata$xf_cfu_per_g)
reinfdata$genotype <- factor(reinfdata$genotype)
reinfdata$trt <- factor(reinfdata$trt)
summary(reinfdata)

#### Where are the NAs?
reinfdata[is.na(reinfdata$xf_cfu_per_g),]
# NAs for xf_cfu_per_g are mostly controls and dead plants
# For sample 9-007S-4, xf_cfu_d0 was too contaminated but d1 and d2 were 0 could be imputed, but I'll remove it
# Remove NAs and controls
reinfdata <- reinfdata %>% dplyr::filter(., trt != "CAB-") %>% dplyr::filter(!is.na(xf_cfu_per_g)) 
# Remove samples cultured on 2017-07-21 because all were negative, something went wrong, not sure what
reinfdata <- reinfdata %>% dplyr::filter(xf_plant_date_cultured != "2017-07-21")

#### Label infection treatments according to reps
reinfdata$inoc.time <- factor(with(reinfdata, ifelse(rep <= 4, "inoc1",
                                                    ifelse(rep >= 9, "inoc2",
                                                           "inoc1-2"))))
summary(reinfdata)

#### Do I have positive culture from all of the plants?
sumcfu <- reinfdata %>% group_by(genotype, rep) %>% summarise(sum = sum(xf_cfu_per_g, na.rm = TRUE))
print.data.frame(sumcfu)

# Check how many times plants were negative that were positive previously
reinflist <- reinfdata %>% dplyr::select(week, genotype, rep, xf_cfu_per_g, inoc.time) %>% data.table() %>% split(., by = c("rep", "genotype"))
reinflist
# Hard to say when zeros were false negatives, can't justify altering them



########################################################################################################
#### Analysis
#### Analysis of re-infection data using Poisson GLMM
# Random effect on plantID because of repeated measures -- is very hard to fit though -- need to re-scale
reinfdata$plantID <- factor(with(reinfdata, paste(genotype, trt, rep, sep = "-")))
# Re-scale week variable
reinfdata$std.week <- standardize(reinfdata$week)
# Remove week == 5 points because I'm lacking data for 092 and 094
reinfdata <- reinfdata %>% dplyr::filter(week != 5)

reinfMod1 <- glmer(xf_cfu_per_g ~ week*genotype + inoc.time + inoc.time:genotype + (1|plantID), data = reinfdata, family = "poisson",
                   control = glmerControl(optimizer = "bobyqa"))
reinfMod2 <- glmer(xf_cfu_per_g ~ week*trt + inoc.time + inoc.time:trt + (1|plantID), data = reinfdata, family = "poisson",
                   control = glmerControl(optimizer = "bobyqa"))
plot(simulateResiduals(reinfMod1)) ### Not a good model
summary(reinfMod1)



########################################################################################################
#### Plotting
#### Summarize and plot reinfection data
reinfSummary <- reinfdata %>% mutate(xfpoplog = log10(xf_cfu_per_g + 1)) %>% 
  group_by(week, trt, inoc.time) %>% 
  summarise(mean = mean(xfpoplog, na.rm = TRUE),
            n = sum(!is.na(xfpoplog)),
            se = sd(xfpoplog)/sqrt(n))
print.data.frame(reinfSummary)


#### Plotting mean xf pops
meanreinfplot <- ggplot(data=reinfSummary, aes(x=week, y=mean)) +
  geom_line(aes(linetype=inoc.time, colour = trt), size=1.25) +
  geom_point(aes(shape=inoc.time, colour = trt), size=3.5) +
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), width=0.2) +
  scale_x_continuous(name = "Weeks post inoculation", 
                     breaks = unique(reinfSummary$week)) + 
  scale_y_continuous(name = "Xylella populations in plant (log10)",
                     limits = c(0,10)) +
  # ylab("% insects on source plant") + 
  # ylim(c(0,100)) +
  # xlab("Weeks post inoculation") +
  theme_bw(base_size=18) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank()) 

meanreinfplot

ggsave("results/figures/2017_figures/reinfection_line_plot_2017.jpg", plot = meanreinfplot,
       width = 7, height = 7, units = "in")


#### Plotting genotypes separately using facetting
reinfSummary2 <- reinfdata %>% mutate(xfpoplog = log10(xf_cfu_per_g + 1)) %>% 
  group_by(week, genotype, trt, inoc.time) %>% 
  summarise(mean = mean(xfpoplog, na.rm = TRUE),
            n = sum(!is.na(xfpoplog)),
            se = sd(xfpoplog)/sqrt(n))
reinfSummary2$se[is.na(reinfSummary2$se)] <- 0 # 26wk-092S only has one observation; set SE = 0 since xf_cfu = 0
print.data.frame(reinfSummary2)

# Plot of reinfection dynamis by genotype and inoculation time (inoc.time)
reinfplot <- ggplot(data=reinfSummary2, aes(x=week, y=mean)) +
  geom_line(aes(linetype=inoc.time, colour = trt), size=1.25) +
  geom_point(aes(shape=inoc.time, colour = trt), size=3.5) +
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), width=0.2) +
  facet_wrap(~genotype) +
  scale_x_continuous(name = "Weeks post inoculation", 
                     breaks = unique(reinfSummary2$week)) + 
  scale_y_continuous(name = "Xylella populations in plant (log10)",
                     limits = c(0,10)) +
  theme_bw(base_size=18) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank()) 

reinfplot


## BW plot for paper
reinfplotBW <- ggplot(data=reinfSummary2, aes(x=week, y=mean)) +
  geom_line(aes(linetype=inoc.time), colour = "black", size=1.25) +
  geom_point(aes(shape=inoc.time), colour = "black", size=3) +
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), width=0.2) +
  facet_wrap(~genotype) +
  scale_x_continuous(name = "Weeks post inoculation", 
                     breaks = unique(reinfSummary2$week)) + 
  scale_y_continuous(name = "Xylella population size (log10, CFU/mg)",
                     limits = c(0,10)) +
  scale_shape_manual(values = c(0,1,2)) +
  labs(linetype = "Treatment", shape = "Treatment") +
  theme_bw(base_size=14) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank(),
        strip.background = element_blank()) 

reinfplotBW

ggsave("results/figures/2017_figures/reinfection_genotype_line_plot_2017_BW.jpg", plot = reinfplotBW,
       width = 7, height = 7, units = "in")




#### Plotting the dyanmics of individual plants
reinfplantplot <- ggplot(data=reinfdata, aes(x = week, y = log10(xf_cfu_per_g+1))) +
  geom_point(aes(shape=inoc.time, colour = trt), size=3.5) +
  facet_wrap(~genotype) +
  scale_x_continuous(name = "Weeks post inoculation", 
                     breaks = unique(reinfdata$week)) + 
  scale_y_continuous(name = "Xylella CFU per g plant tissue (log 10)") +
  # ylab("% insects on source plant") + 
  # ylim(c(0,100)) +
  # xlab("Weeks post inoculation") +
  theme_bw(base_size=18) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank()) 
reinfplantplot
# ggsave("results/figures/vector_xf_line_plot.jpg", plot = vectorxfplot,
#        width = 7, height = 7, units = "in")




########################################################################################################
########################################################################################################
#### Analysis and plotting of just 2nd inoculation (weeks 21 and 26) using ANCOVA
inoc2dat <- reinfdata %>% dplyr::filter(week >= 21 & inoc.time != "inoc1")
summary(inoc2dat)


## GLM Poisson and quasiPoisson don't work well
## Try negative binomial
inoc2NB <- glm.nb((xf_cfu_per_g+1) ~ week*trt*inoc.time, data = inoc2dat,
                  control = glm.control(epsilon = 1e-16, maxit = 100))
summary(inoc2NB)
## Plotting options not available for glm.nb
## Check constant error variance
residNB <- resid(inoc2NB)
predNB <- predict(inoc2NB)
plot(x = predNB, y = residNB)
## Not constant

## Linear model
boxcox((xf_cfu_per_g+1) ~ week*trt*inoc.time, data = inoc2dat, lambda = seq(-2, 2, by=0.5))
inoc2mod <- lm(sqrt(xf_cfu_per_g+1) ~ week*trt*inoc.time, data = inoc2dat)
plot(simulateResiduals(inoc2mod)) 
summary(inoc2mod)

## Get 95% CI
tidy(inoc2NB, conf.int = TRUE)


## RESULTS: Both NB and linear model return similar results: week main effect is significant, other terms are not significant


#### Plotting
inoc2Summary <- inoc2dat %>% mutate(xfpoplog = log10(xf_cfu_per_g + 1)) %>% 
  group_by(week, trt, inoc.time) %>% 
  summarise(mean = mean(xfpoplog, na.rm = TRUE),
            n = sum(!is.na(xfpoplog)),
            se = sd(xfpoplog)/sqrt(n))
print.data.frame(inoc2Summary)


#### Plotting mean xf pops
meaninoc2plot <- ggplot(data=inoc2Summary, aes(x=week, y=mean)) +
  ## inoc.time: inoc1-2 = open circles; inoc1 = open triangles
  ## trt: R = solid line; S = dashed line
  geom_line(aes(linetype=trt, colour = inoc.time), size=1.25) +
  geom_point(aes(shape=inoc.time), colour = "black", size=3) +
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), width=0.2) +
  scale_x_continuous(name = "Weeks post inoculation", 
                     breaks = unique(inoc2Summary$week)) + 
  scale_y_continuous(name = "Xylella population size (log10, CFU/g)",
                     limits = c(0,10)) +
  scale_color_manual(values = c("black", "black")) + 
  scale_shape_manual(values = c(1,2)) +
  theme_bw(base_size=14) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank(),
        legend.position = "none")

meaninoc2plot

ggsave("results/figures/2017_figures/reinfection_inoc2_line_plot_2017.jpg", plot = meaninoc2plot,
       width = 7, height = 7, units = "in")



####################################################################################################
#### Cleaning up data set for Chris Wallis

pdr1DataURL <- gs_url("https://docs.google.com/spreadsheets/d/14uJLfRL6mPrdf4qABeGeip5ZkryXmMKkan3mJHeK13k/edit?usp=sharing",
                      visibility = "private")
reinfdataGS <- gs_read(pdr1DataURL, ws = "Re-infection_culturing_data")


reinfdata <- reinfdataGS
str(reinfdata)
summary(reinfdata)
# Make xf_cfu_per_g a numeric variable and genotype and trt factors
reinfdata$xf_cfu_per_g <- as.numeric(reinfdata$xf_cfu_per_g)
reinfdata$genotype <- factor(reinfdata$genotype)
reinfdata$trt <- factor(reinfdata$trt)
summary(reinfdata)

#### Where are the NAs?
reinfdata[is.na(reinfdata$xf_cfu_per_g),]
# NAs for xf_cfu_per_g are mostly controls and dead plants
# For sample 9-007S-4, xf_cfu_d0 was too contaminated but d1 and d2 were 0 could be imputed, but I'll remove it

## For Chris' analysis, he only needs genotypes 092 and 094, as these were the only to get assayed for phytochemistry 
## Also remove unnecessary columns
reinfdataCW <- reinfdata %>% 
  dplyr::filter(., genotype == "092" | genotype == "094") %>%
  dplyr::select(week, genotype, trt, rep, xf_plant_date_cultured, xf_cfu_per_g, notes)

#### Label infection treatments according to reps
reinfdataCW$inoc.time <- factor(with(reinfdataCW, ifelse(rep <= 4, "inoc1",
                                                         ifelse(rep >= 9, "inoc2",
                                                                "inoc1-2"))))
summary(reinfdataCW)

write.csv(reinfdataCW, "output/xf_populations_for_Chris.csv", row.names = FALSE)
