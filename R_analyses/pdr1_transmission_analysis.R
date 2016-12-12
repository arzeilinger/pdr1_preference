#### ANALYSIS OF PDR1 TRANSMISSION DATA

my.packages <- c("tidyr", "dplyr", "data.table", "xlsx", "ggplot2", "MASS")
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
prop_pd_leaves <- leafdata %>% with(., ifelse(n_leaves_source == 0, 1, n_pd_leaves/n_leaves_source))
leafdata$pd_index2 <- leafdata %>% with(., prop_pd_leaves + n_ms_petioles)
plot(x = leafdata$pd_index, y = leafdata$pd_index2)


################################################################################################################
#### Constructing vector infection index
# In each trial, I have Xf pops for each of the vectors. I could include in transmission model:
# total Xf population among all vectors
# proportion of vectors infected (b/c some cages have fewer than 8 vectors)
# an evenness index of Xf pops amongs vectors
# Matt's transmission parameters paper might have something to say about this
# multiple possibilities can be evaluated using AIC



####################################################################################################################
#### Import and combine culturing data, leaf data, and preference data
# Import culturing data
culturedata <- read.xlsx("data/pdr1_culturing_data.xlsx", sheetName = "Infection data")
# Remove repeat culturing data
culturedata <- culturedata %>% dplyr::filter(., notes != "repeat culture" | is.na(notes))
str(culturedata)
culturedata

# Import preference data: estimates of rate parameters from CM model for each cage
# choice 1 = source plants, choice 2 = test plants
# Ignoring variance around parameters.
# TO DO: check on convergence issues
paramDataCage <- readRDS("output/CMM_rate_parameters_per_cage.rds")
# Reshape paramDataCage
paramDataCage <- paramDataCage[,c("parameter", "estimate", "week.cage")] %>% spread(., key = parameter, value = estimate)

# Merge culturing and leaf data on week-trt-rep combination
transdata <- inner_join(culturedata, leafdata, by = c("week", "trt", "rep")) 

# Merge with preference data
transdata$week.cage <- transdata %>% with(., paste(week, trt, rep, sep=""))
transdata <- left_join(transdata, paramDataCage, by = "week.cage")
str(transdata)

saveRDS(transdata, file = "output/pdr1_transmission_preference_dataset.rds")


#############################################################################################################
#### Transmission models
transdata <- readRDS("output/pdr1_transmission_preference_dataset.rds")
transdata$source.cfu.per.g <- factor2numeric(transdata$source.cfu.per.g)

# Remove the second 12-week trials and NAs
transdata <- transdata %>% dplyr::filter(., week != 12.2) 
                                         !is.na(pd_index) & !is.na(pd_index2))

#### Model selection on PD symptom index
pdMod1 <- glm(test.plant.infection ~ week*trt + log10(source.cfu.per.g+1) + p1 + p2 + mu1 + mu2 + pd_index, data = transdata, family = "binomial")
pdMod2 <- glm(test.plant.infection ~ week*trt + log10(source.cfu.per.g+1) + p1 + p2 + mu1 + mu2 + pd_index2, data = transdata, family = "binomial")
# I think quasibinomial distribution might be better but then it doesn't calculate an AIC value. Need to look into this.
AIC(pdMod1, pdMod2)
plot(pdMod2)
summary(pdMod2)
# PD indices are essentially the same; go with Arash's index
# Re-load data set to retain data points that are NA for pd indices

transdata <- readRDS("output/pdr1_transmission_preference_dataset.rds") %>% dplyr::filter(., week != 12.2)
transdata$source.cfu.per.g <- factor2numeric(transdata$source.cfu.per.g)


#### Model selection on preference rate parameters
FullModel <- glm(test.plant.infection ~ week*trt + log10(source.cfu.per.g+1) + p1 + p2 + mu1 + mu2 + pd_index, data = transdata, family = "binomial")
prefModChoice1 <- glm(test.plant.infection ~ week*trt + log10(source.cfu.per.g+1) + p1 + mu1 + pd_index, data = transdata, family = "binomial")
prefModp1 <- glm(test.plant.infection ~ week*trt + log10(source.cfu.per.g+1) + p1 + pd_index, data = transdata, family = "binomial")
prefModmu1 <- glm(test.plant.infection ~ week*trt + log10(source.cfu.per.g+1) + mu1 + pd_index, data = transdata, family = "binomial")
prefModnull <- glm(test.plant.infection ~ week*trt + log10(source.cfu.per.g+1) + pd_index, data = transdata, family = "binomial")

AIC(FullModel, prefModChoice1, prefModp1, prefModmu1, prefModnull)
# Including only mu1 model is best
summary(prefModmu1)


#### Model selection to include source xf pop and PD index
# start from prefModmu1
noSourceMod <- glm(test.plant.infection ~ week*trt*mu1*pd_index, data = transdata, family = "binomial")
noPDMod <- glm(test.plant.infection ~ week*trt*log10(source.cfu.per.g+1)*mu1, data = transdata, family = "binomial")
trtprefMod <- glm(test.plant.infection ~ week*trt*mu1, data = transdata, family = "binomial")
trtMod <- glm(test.plant.infection ~ week*trt, data = transdata, family = "binomial")

AIC(noSourceMod, noPDMod, trtprefMod, trtMod)
# trtMod seems best
summary(trtMod)


#### symptom data
sympMod <- glm(pd_index ~ week*trt*log10(source.cfu.per.g+1), data = transdata, family = "quasipoisson")
plot(sympMod)
summary(sympMod)


#### source xf pop
boxcox(source.cfu.per.g + 1 ~ week*trt, data = transdata, lambda = seq(-2, 2, by=0.5))

sourcexfMod <- lm(log10(source.cfu.per.g+1) ~ week*trt, data = transdata)
plot(sourcexfMod)
summary(sourcexfMod)

#######################################################################################################
#### Plotting
transSummary <- transdata %>% group_by(., week, trt) %>% 
  summarise(percInfected = 100*sum(test.plant.infection, na.rm = TRUE)/sum(!is.na(test.plant.infection)),
            meanPD = mean(pd_index, na.rm = TRUE),
            sePD = sd(pd_index, na.rm = TRUE)/sqrt(sum(!is.na(pd_index))),
            meancfu = mean(log10(source.cfu.per.g+1)),
            secfu = sd(log10(source.cfu.per.g+1))/sqrt(sum(!is.na(source.cfu.per.g))))

# Transmission plot
transplot <- ggplot(data=transSummary, aes(x=week, y=percInfected)) +
  # geom_bar(position=position_dodge(), stat="identity", 
  #          aes(fill=trt)) +
  # geom_hline(aes(yintercept=50), linetype="dashed") +
  geom_line(aes(linetype=trt), size=1.25) +
  geom_point(aes(shape=trt), size=2.5) +
  #geom_errorbar(aes(ymax=meancfu+secfu, ymin=meancfu-secfu), width=0.2) +
  scale_x_continuous(name = "Weeks post inoculation", 
                     breaks = c(3,8,12)) + 
  scale_y_continuous(name = "Percent transmission",
                     limits = c(0,100)) +
  # ylab("% insects on source plant") + 
  # ylim(c(0,100)) +
  # xlab("Weeks post inoculation") +
  theme_bw(base_size=18) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank()) 

transplot
ggsave("results/figures/transmission_line_plot.jpg", plot = transplot,
       width = 7, height = 7, units = "in")


#### Symptoms plot
symptomPlot <- ggplot(data=transSummary, aes(x=week, y=meanPD)) +
  # geom_bar(position=position_dodge(), stat="identity", 
  #          aes(fill=trt)) +
  # geom_hline(aes(yintercept=50), linetype="dashed") +
  geom_line(aes(linetype=trt), size=1.25) +
  geom_point(aes(shape=trt), size=2.5) +
  geom_errorbar(aes(ymax=meanPD+sePD, ymin=meanPD-sePD), width=0.2) +
  scale_x_continuous(name = "Weeks post inoculation", 
                     breaks = c(3,8,12)) + 
  scale_y_continuous(name = "Mean PD symptoms index") +
  # ylab("% insects on source plant") + 
  # ylim(c(0,100)) +
  # xlab("Weeks post inoculation") +
  theme_bw(base_size=18) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank()) 

symptomPlot
ggsave("results/figures/pd_symptom_line_plot.jpg", plot = symptomPlot,
       width = 7, height = 7, units = "in")


#### Xf pops in source plant plot
sourcexfplot <- ggplot(data=transSummary, aes(x=week, y=meancfu)) +
  # geom_bar(position=position_dodge(), stat="identity", 
  #          aes(fill=trt)) +
  # geom_hline(aes(yintercept=50), linetype="dashed") +
  geom_line(aes(linetype=trt), size=1.25) +
  geom_point(aes(shape=trt), size=2.5) +
  geom_errorbar(aes(ymax=meancfu+secfu, ymin=meancfu-secfu), width=0.2) +
  scale_x_continuous(name = "Weeks post inoculation", 
                     breaks = c(3,8,12)) + 
  scale_y_continuous(name = "Xylella CFU/g plant tissue (log10)") +
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
ggsave("results/figures/source_xf_line_plot_log.jpg", plot = sourcexfplot,
       width = 7, height = 7, units = "in")



#### Xf pops in source plant plot scatter plot
xfscatterplot <- ggplot(data=transdata, aes(x=week, y=(source.cfu.per.g/1000000))) +
  # geom_bar(position=position_dodge(), stat="identity", 
  #          aes(fill=trt)) +
  # geom_hline(aes(yintercept=50), linetype="dashed") +
  #geom_line(aes(linetype=trt), size=1.25) +
  geom_point(aes(shape=trt), size=2.5) +
  #geom_errorbar(aes(ymax=meancfu+secfu, ymin=meancfu-secfu), width=0.2) +
  scale_x_continuous(name = "Weeks post inoculation", 
                     breaks = c(3,8,12)) + 
  scale_y_continuous(name = "Xylella CFU/g plant tissue (x 1,000,000)") +
  # ylab("% insects on source plant") + 
  # ylim(c(0,100)) +
  # xlab("Weeks post inoculation") +
  theme_bw(base_size=18) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank()) 

xfscatterplot
ggsave("results/figures/source_xf_pop_scatterplot.jpg", plot = xfscatterplot,
       width = 7, height = 7, units = "in")
