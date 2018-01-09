#### Analysis of PdR1 preference-transmission experiment 2017

rm(list = ls())
# Load packages
my.packages <- c("tidyr", "dplyr", "data.table", "openxlsx", "ggplot2",
                 "MASS", "logistf", "multcomp", "bbmle", "lme4", "googlesheets")
lapply(my.packages, require, character.only = TRUE)

source("R_functions/factor2numeric.R")


#### Import culturing data
# Import from local .xlsx file
#sourcedata <- read.xlsx("data/2017_data/PdR1_2017_preference-transmission_experiment_data.xlsx", sheet = "source_plant_culturing", detectDates = TRUE)

# Import from Googlesheets
pdr1DataURL <- gs_url("https://docs.google.com/spreadsheets/d/14uJLfRL6mPrdf4qABeGeip5ZkryXmMKkan3mJHeK13k/edit?usp=sharing",
                     visibility = "private")
sourcedataGS <- gs_read(pdr1DataURL, ws = "source_plant_culturing")
sourcedata <- sourcedataGS
str(sourcedata)
summary(sourcedata)

# Convert classes for columns
sourcedata$block <- factor(sourcedata$block)
sourcedata$genotype <- factor(sourcedata$genotype)
sourcedata$trt <- factor(sourcedata$trt)
sourcedata$xf_cfu_per_g <- as.numeric(sourcedata$xf_cfu_per_g)

# Remove NAs and Control samples
print.data.frame(sourcedata[is.na(sourcedata$xf_cfu_per_g), ])
sourcedata <- sourcedata %>% dplyr::filter(!is.na(xf_cfu_per_g))


# #### Which plants were negative at least once?
# # make a plant ID
# sourcedata$plantID <- with(sourcedata, paste(week, block, genotype, trt, rep, sep = "-"))
# 
# negativeSamples <- sourcedata %>% filter(., xf_cfu_d0 == 0 & xf_cfu_d1 == 0 & xf_cfu_d2 == 0 & genotype != "CAB-") %>%
#   dplyr::select(., plantID, times_cultured, starts_with("xf_cfu"), notes) %>% arrange(., plantID, times_cultured)
# 
# # get all times cultured for each plant that was negative once
# negativePlants <- sourcedata[sourcedata$plantID %in% negativeSamples$plantID,] %>% 
#   dplyr::select(., plantID, times_cultured, starts_with("xf_cfu"), notes) %>% arrange(., plantID, times_cultured)


#### For plants that were cultured multiple times, get the best estimate of Xf pop
# Remove unnecessary columns and use spread() to get each plant on a single row
sourcedata <- sourcedata %>% mutate(plantID = factor(paste(week, block, genotype, rep, sep = "-"))) %>%
  mutate(times_cultured = paste("culture", times_cultured, sep = "_"))
sourcedata2 <- sourcedata %>% dplyr::select(-xf_plant_sample_mass, -xf_cfu_d0, -xf_cfu_d1, -xf_cfu_d2, -xf_plant_date_cultured, -notes, -notes2)
sourcedata2 <- sourcedata2 %>% spread(., key = times_cultured, value = xf_cfu_per_g)
print.data.frame(sourcedata2)
length(unique(sourcedata2$plantID)) == nrow(sourcedata2) # Check that spread() correctly prdouced a unique row for each plant ID; if so should be "TRUE"
# Make sure that the latest re-culture was always the best
recultures <- sourcedata2 %>% dplyr::filter(!is.na(culture_2))
# Things are complicated, latest re-culture wasn't always the best. Use complicated ifelse statements
# Make the final xf pop estimates
sourcedata2$xfpop <- with(sourcedata2, 
                          ifelse(culture_1 == 0 & (culture_2 == 0 | is.na(culture_2)) & (culture_3 == 0 | is.na(culture_3)), 0,
                                 ifelse((culture_1 == 0 | is.na(culture_1)) & (culture_2 > 0 | culture_3 > 0), 100,
                                        culture_1)))
print.data.frame(sourcedata2)


#### Analysis of source plant xf populations
# Linear model with transformation
boxcox(xfpop + 1 ~ block + week*genotype, data = sourcedata2, lambda = seq(-2, 2, by=0.5))
sourcepopMod1 <- lm(sqrt(xfpop) ~ block + week*genotype, data = sourcedata2)
plot(sourcepopMod1)
summary(sourcepopMod1)
# Generalized linear model with quasipoisson distribution
poispopMod <- glm(xfpop ~ block + week*genotype, data = sourcedata2, family = "quasipoisson")
plot(poispopMod)
summary(poispopMod)

#### Plotting source plant xf populations
sourceSummary <- sourcedata2 %>% mutate(logxfpop = log10(xfpop + 1)) %>% 
  group_by(week, genotype, trt) %>% 
  summarise(mean = mean(logxfpop),
            se = sd(logxfpop)/sqrt(length(logxfpop)))

# Xf pops in source plant plot
sourcexfplot <- ggplot(data=sourceSummary, aes(x=week, y=mean)) +
  geom_line(aes(linetype=genotype, colour = trt), size=1.25) +
  geom_point(aes(shape=genotype, colour = trt), size=3.5) +
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), width=0.2) +
  scale_x_continuous(name = "Weeks post inoculation", 
                     breaks = c(2,5,8,14)) + 
  scale_y_continuous(name = "Xylella populations in source plant (log10)",
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

sourcexfplot
ggsave("results/figures/2017_figures/source_xf_line_plot_2017.jpg", plot = PDplot,
       width = 7, height = 7, units = "in")

