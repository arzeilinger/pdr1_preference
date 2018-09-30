#### Looking at temperature data from PdR1 2017 experiment
#### Data are from data loggers placed in cages during preference trials

rm(list = ls())
# libraries
# loading dtplyr that replaces dplyr and data.table
my.packages <- c("tidyr", "dplyr", "ggplot2", "data.table", "lubridate", "openxlsx")
lapply(my.packages, require, character.only = TRUE)

# Define working directory for temperature files
tdir <- "data/2017_data/temperature_data_2017/"

## Get file names for IGH7 temperature data
#tempFiles <- list.files(path = tdir, pattern = "IGH7")
# Need to fix dates in the two .xlsx files. For now, just use the .txt files
tempFiles <- list.files(path = tdir, pattern = ".txt")
## Load data from .txt files
tList <- vector("list", length(tempFiles))

for(i in 1:length(tempFiles)){
  if(grepl(x = tempFiles[i], pattern = ".txt")){
    tList[[i]] <- read.table(paste(tdir, tempFiles[i], sep = ""), sep = ",", fill = TRUE, skip = 1, col.names = c("IGH7", "time", "temp", "Serial_number"))[,1:3]
  } else {
    tList[[i]] <- read.xlsx(paste(tdir, tempFiles[i], sep = ""), detectDates = TRUE)
  }
}

str(tList)

tdat <- tList %>% rbindlist() %>% as.data.frame()
tdat$time <- as_datetime(tdat$time)
tdat <- tdat %>% arrange(time)

with(tdat, plot(x = time, y = temp, type = "l"))

# Calculate daily max, mean, and min
tdat$day <- day(tdat$time)
tdat$month <- month(tdat$time)
daySummary <- tdat %>% group_by(month, day) %>% summarise(tmax = max(temp),
                                                          tmean = mean(temp),
                                                          tmin = min(temp))
daySummary$date <- as_date(paste("2017", daySummary$month, daySummary$day, sep = "-"))
print.data.frame(daySummary)

#### Determine the dates of each set of trials based on bgss_count data
prefdata <- read.xlsx("data/2017_data/PdR1_2017_preference-transmission_experiment_data.xlsx", sheet = "BGSS_count_data", detectDates = TRUE)
str(prefdata)
## Data.frame with the start dates for each trial set and the week-block
trials <- data.frame(startDate = unique(prefdata$start_date),
                     trialBlock = unique(paste(prefdata$week, prefdata$block, sep = "-")))

#### For each trial set, extract temperature data for each day
trialTempList <- vector("list", nrow(trials))

for(i in 1:length(trials$trialBlock)){
  startDate.i <- trials$startDate[i]
  trialDays.i <- c(startDate.i, startDate.i + 1, startDate.i + 2, startDate.i + 3, startDate.i + 4) 
  trialTemps.i <- daySummary[daySummary$date %in% trialDays.i,]
  trialTemps.i$trialBlock <- rep(trials$trialBlock[i], nrow(trialTemps.i))
  trialTempList[[i]] <- trialTemps.i
}

#### rbind into data frame then split trialBlock into week and block columns
trialTempData <- trialTempList %>% rbindlist() %>% as.data.frame() %>% 
  mutate(week = as.numeric(tstrsplit(trialBlock, "-")[[1]]), block = as.numeric(tstrsplit(trialBlock, "-")[[2]]))

#### Gather max, mean, and min temperature data and summarise for each week
tempSummary <- trialTempData %>% gather(key = measure, value = temperature, tmax:tmin) %>%
  group_by(week, measure) %>% 
  summarise(meanT = mean(temperature, na.rm = TRUE),
            seT = sd(temperature, na.rm = TRUE)/sqrt(length(temperature)))
  
tempSummary


#### Plot temperature data
tempPlot <- ggplot(data=tempSummary, aes(x=week, y=meanT)) +
  geom_point(aes(colour = measure), size=3.5) +
  geom_errorbar(aes(ymax=meanT+seT, ymin=meanT-seT), width=0.2) +
  scale_x_continuous(name = "Weeks post inoculation", 
                     breaks = c(2,5,8,14), limits = c(1,15)) + 
  scale_y_continuous(name = "Mean daily temperature (degrees C)",
                     breaks = seq(20,50,by=10), limits = c(20,50)) +
  theme_bw(base_size=18) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank()) 

tempPlot

ggsave("results/figures/2017_figures/temperature_plot_2017.jpg", plot = tempPlot,
       width = 7, height = 7, units = "in")


