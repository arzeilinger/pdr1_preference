#### Analysis of PdR1 preference experiment data

rm(list = ls())
# libraries
my.packages <- c("openxlsx", "tidyr", "dplyr", "data.table", "ggplot2", "lattice")
lapply(my.packages, require, character.only = TRUE)

# Load preference data from .xlsx file
prefdata <- read.xlsx("data/pdr1_preference_data.xlsx", sheet = "data")
prefdata$cage <- paste(prefdata$trt, prefdata$rep, sep = "")
str(prefdata)

data1 <- prefdata[prefdata$rep == 1,]

ggplot(data=prefdata,aes(x=time_from_start_hr, y=source_plant)) + 
  geom_line(aes(group=rep, color=trt), size=1.25) +
  scale_x_continuous("Time from start (hr)") + 
  scale_y_continuous("Number of BGSS on source plant")

xyplot(source_plant ~ time_from_start_hr|trt, 
       groups = rep, 
       data = prefdata, type = "l")


#########################################################################
#### Number of bugs for PCR in each trial
lastObs <- prefdata[prefdata$time_from_start == "8 d",]
nbgss <- data.frame(cage = lastObs$cage, 
                    nbgss = 8 - lastObs$missing)
write.csv(nbgss, file = "number_of_bgss_for_pcr.csv", row.names = FALSE)


##########################################################################
#### Mean counts for each treatment, week, and time point
meanPlant <- prefdata %>% group_by(week, trt, time_from_start_hr) %>% 
  summarise(meanSourcePlant = mean(source_plant), meanTestPlant = mean(test_plant),
            sumSourcePlant = sum(source_plant), sumTestPlant = sum(test_plant),
            propTestPlant = meanTestPlant/(meanTestPlant+meanSourcePlant))

ggplot(data=meanPlant,aes(x=time_from_start_hr, y=meanSourcePlant)) + 
  geom_line(aes(color=trt), size=1.25) +
  scale_x_continuous("Time from start (hr)") + 
  scale_y_continuous("Mean number of BGSS on source plant")

# Proportion of bugs on test plant
prefTimeSeries <- ggplot(data=meanPlant[meanPlant$week != 12.2,],aes(x=time_from_start_hr, y=propTestPlant)) + 
                         geom_line(aes(color=interaction(week,trt)), size=1.25) +
                         scale_x_continuous("Time from start (hr)") + 
                         scale_y_continuous("Proportion of BGSS on test plant")
prefTimeSeries #+ scale_x_continuous(limits = c(0,30))

# Total number of bugs
prefTimeSeries <- ggplot(data=meanPlant[meanPlant$week != 12.2,],aes(x=time_from_start_hr, y=sumTestPlant)) + 
  geom_line(aes(color=interaction(week,trt)), size=1.25) +
  scale_x_continuous("Time from start (hr)") + 
  scale_y_continuous("Mean number of BGSS on test plant")
prefTimeSeries #+ scale_x_continuous(limits = c(0,30))

#############################################################################
#### Counts at the end of trials
lastObs <- prefdata[prefdata$time_from_start == "8 d",1:12]
lastObs$week.trt <- factor(paste(lastObs$week, lastObs$trt, sep = ""))
lastObs$propSource <- lastObs$source_plant/(lastObs$source_plant + lastObs$test_plant)
str(lastObs)

#### Analysis of bug counts at end of trial
lastCountaov <- lm(asin(sqrt(propSource)) ~ week*trt, data = lastObs)
plot(lastCountaov)
summary(lastCountaov)

#### t tests for differences from 0.5 null
#### for each week by trt interaction
for(i in 1:length(levels(lastObs$week.trt))){
  week.trt.i <- levels(lastObs$week.trt)[i]
  lastObs.i <- lastObs[lastObs$week.trt == week.trt.i,]
  print(week.trt.i)
  print(t.test(lastObs.i[,"propSource"], mu = 0.5))
}

#### Plotting percent
meanPropSource <- lastObs %>% group_by(week, trt) %>% 
  summarise(mean.perc = mean(propSource)*100)
meanPropSource # Proportion of bugs on the test plant, for each treatment and week

percBugsPrefplot <- ggplot(data=meanPropSource, aes(x=week, y=mean.perc)) +
  # geom_bar(position=position_dodge(), stat="identity", 
  #          aes(fill=trt)) +
  geom_hline(aes(yintercept=50), linetype="dashed") +
  geom_line(aes(linetype=trt), size=1.25) +
  geom_point(aes(shape=trt), size=2.5) +
  scale_x_continuous(name = "Weeks post inoculation", 
                     breaks = c(3,8,12)) + 
  scale_y_continuous(limits = c(0,100), name = "% insects on source plant") +
  # ylab("% insects on source plant") + 
  # ylim(c(0,100)) +
  # xlab("Weeks post inoculation") +
  theme_bw(base_size=18) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank()) 
ggsave("output/figures/prefplot.jpg", plot = percBugsPrefplot,
       width = 7, height=7, units="in")
percBugsPrefplot

