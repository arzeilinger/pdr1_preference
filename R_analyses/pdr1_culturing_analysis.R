#### Analysis of Culturing data

rm(list = ls())
# libraries
my.packages <- c("xlsx", "tidyr", "dplyr", "data.table", "ggplot2", "lattice", "MASS")
lapply(my.packages, require, character.only = TRUE)

# Load preference data from .xlsx file
plantdata <- read.xlsx("data/pdr1_culturing_data.xlsx", sheetName = "Infection data")
str(plantdata)


#### Analysis of source plant xf pops
boxcox(source.cfu.per.g + 1 ~ week*trt, data = plantdata, lambda = seq(-2, 2, by=0.5))

xf.pop.mod <- lm(source.cfu.per.g ~ week*trt, data = plantdata)
plot(xf.pop.mod)
summary(xf.pop.mod)

#### Plotting source plant xf pops
#### NEED TO EDIT THIS CODE FOR SOURCE PLANT XF ###############################################

cfumeans <- plantdata %>% group_by(week, trt) %>% 
  summarise(mean = mean(source.cfu.per.g/1000000), se = sd(source.cfu.per.g/1000000)/sqrt(length(source.cfu.per.g)))
cfumeans # Proportion of bugs on the test plant, for each treatment and week

sourcexfplot <- ggplot(data=cfumeans, aes(x=week, y=mean)) +
  # geom_bar(position=position_dodge(), stat="identity", 
  #          aes(fill=trt)) +
  # geom_hline(aes(yintercept=50), linetype="dashed") +
  geom_line(aes(linetype=trt), size=1.25) +
  geom_point(aes(shape=trt), size=2.5) +
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), width=0.2) +
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
ggsave("output/figures/source_xf_line_plot.jpg", plot = sourcexfplot,
       width = 7, height = 7, units = "in")
sourcexfplot

################################################################################
##### Transmission data
transdata <- plantdata[!is.na(plantdata$test.plant.infection),]

percTrans <- transdata %>% group_by(trt) %>% summarise(percent = sum(test.plant.infection)/length(test.plant.infection)*100)
percTrans
