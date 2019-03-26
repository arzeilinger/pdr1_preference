#### Analysis of PdR1 preference-transmission experiment 2017

rm(list = ls())
# Load packages
my.packages <- c("tidyr", "dplyr", "data.table", "openxlsx", "ggplot2",
                 "MASS", "logistf", "multcomp", "bbmle", "lme4", "googlesheets", "DHARMa")
lapply(my.packages, require, character.only = TRUE)

source("R_functions/factor2numeric.R")


#### Munging source plant culturing data moved to "munging_pdr1_preference_transmission_data.R"
## Load complete data set
transVCPdata <- readRDS("output/complete_2017_transmission-preference_dataset.rds")


#### Analysis of source plant xf populations
## Linear model with transformation
boxcox(xfpop + 1 ~ block + week*genotype, data = transVCPdata, lambda = seq(-2, 2, by=0.5))
# Best transformation is either sqrt or log; go with square root because the residuals look better
sourcepopMod1 <- lm(sqrt(xfpop) ~ block + week*genotype, data = transVCPdata)
plot(simulateResiduals(sourcepopMod1))
## Residuals don't look very good
summary(sourcepopMod1)

## Generalized linear model with quasipoisson distribution
poispopMod1 <- glm(xfpop ~ block + week*genotype, data = transVCPdata, family = "quasipoisson")
plot(poispopMod1)
# Residuals look about the same as lm() with sqrt()
summary(poispopMod1)

## Negative binomial GLM
sourcexfNB <- glm.nb(xfpop ~ block + week*genotype, data = transVCPdata)
summary(sourcexfNB)
## Results are qualitatively similar to quasipoisson but returns a warning. Quasipoisson also works better for 2016 data; go with that.

#### NOTE ON ALTERNATIVE MODEL: Removing false negatives and all negatives had negligible effects on results.
#### These models had reduced significance of week main effect but were otherwise unchanged.
#### Log10 transformation works better if all negatives are removed


#### Plotting source plant xf populations
sourceSummary <- sourcedata2 %>% mutate(logxfpop = log10(xfpop+1), sqrtxfpop = sqrt(xfpop)) %>% 
  group_by(week, genotype, trt) %>% 
  summarise_at(c("logxfpop", "sqrtxfpop"), funs(mean = mean(.), n = sum(!is.na(.)), se = sd(.)/sqrt(sum(!is.na(.)))))


# Xf pops in source plant plot
sourcexfplot <- ggplot(data=sourceSummary, aes(x=week, y=sqrtxfpop_mean)) +
  geom_line(aes(linetype=genotype, colour = trt), size=1.25) +
  geom_point(aes(shape=genotype, colour = trt), size=3.5) +
  geom_errorbar(aes(ymax=sqrtxfpop_mean+sqrtxfpop_se, ymin=sqrtxfpop_mean-sqrtxfpop_se), width=0.2) +
  scale_x_continuous(name = "Weeks post inoculation", 
                     breaks = c(2,5,8,14)) + 
  scale_y_continuous(name = "Xylella populations in source plant (square root)") +
                    # limits = c(0,10)) +
  # ylab("% insects on source plant") + 
  # ylim(c(0,100)) +
  # xlab("Weeks post inoculation") +
  theme_bw(base_size=18) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank()) 

sourcexfplot

ggsave("results/figures/2017_figures/source_xf_sqrt_line_plot_2017.jpg", plot = sourcexfplot,
       width = 7, height = 7, units = "in")

