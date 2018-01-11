#### Analysis of PdR1 re-infection experiment 2017

rm(list = ls())
# Load packages
my.packages <- c("tidyr", "dplyr", "data.table", "openxlsx", "ggplot2",
                 "MASS", "logistf", "multcomp", "bbmle", "lme4")
lapply(my.packages, require, character.only = TRUE)

source("R_functions/factor2numeric.R")

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
# Remove NAs
reinfdata <- reinfdata %>% dplyr::filter(!is.na(xf_cfu_per_g))
# Remove samples cultured on 2017-07-21 because all were negative, something went wrong, not sure what
reinfdata <- reinfdata %>% dplyr::filter(xf_plant_date_cultured != "2017-07-21")

#### Label infection treatments according to reps
reinfdata$inoc.trt <- factor(with(reinfdata, ifelse(rep <= 4, "inoc1",
                                                    ifelse(rep >= 9, "inoc2",
                                                           "inoc1-2"))))
summary(reinfdata)

#### Do I have positive culture from all of the plants?
reinfdata <- reinfdata %>% dplyr::filter(., trt != "CAB-") %>% dplyr::filter(!is.na(xf_cfu_per_g)) 
sumcfu <- reinfdata %>% group_by(genotype, rep) %>% summarise(sum = sum(xf_cfu_per_g, na.rm = TRUE))
print.data.frame(sumcfu)

# Check how many times plants were negative that were positive previously
reinflist <- reinfdata %>% dplyr::select(week, genotype, rep, xf_cfu_per_g, inoc.trt) %>% data.table() %>% split(., by = c("rep", "genotype"))
reinflist
# Hard to say when zeros were false negatives, can't justify altering them

#### Analysis of re-infection data using Poisson GLMM
# Random effect on plantID because of repeated measures -- is very hard to fit though -- need to re-scale
reinfdata$plantID <- factor(with(reinfdata, paste(genotype, trt, rep, sep = "-")))
reinfMod1 <- glmer(xf_cfu_per_g ~ week + genotype + inoc.trt + (1|plantID), data = reinfdata, family = "poisson")
plot(reinfMod1)
summary(reinfMod1)

#### Summarize and plot reinfection data
reinfSummary <- reinfdata %>% mutate(xfpoplog = log10(xf_cfu_per_g + 1)) %>% 
  dplyr::filter(week != 5) %>% # Remove week 5 because it's a mess
  group_by(week, trt, inoc.trt) %>% 
  summarise(mean = mean(xfpoplog, na.rm = TRUE),
            n = sum(!is.na(xfpoplog)),
            se = sd(xfpoplog)/sqrt(n))
print.data.frame(reinfSummary)


#### Plotting mean xf pops
meanreinfplot <- ggplot(data=reinfSummary, aes(x=week, y=mean)) +
  geom_line(aes(linetype=inoc.trt, colour = trt), size=1.25) +
  geom_point(aes(shape=inoc.trt, colour = trt), size=3.5) +
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



#### Plotting the dyanmics of individual plants
reinfplantplot <- ggplot(data=reinfdata, aes(x = week, y = log10(xf_cfu_per_g+1))) +
  #geom_line(aes(colour = trt), size=1.25) +
  geom_point(aes(shape = genotype, colour = trt), size=2.5) +
  #geom_errorbar(aes(ymax=meancfu+secfu, ymin=meancfu-secfu), width=0.2) +
  scale_x_continuous(name = "Weeks post inoculation", 
                     breaks = c(2,5,9)) + 
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
