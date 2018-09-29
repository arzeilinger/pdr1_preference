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
print.data.frame(daySummary)


