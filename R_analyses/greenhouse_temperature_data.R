#### Greenhouse temperature data


my.packages <- c("data.table", "openxlsx", "tidyr", "dplyr", "ggplot2")
lapply(my.packages, require, character.only = TRUE)


#### load data sets
# Most are text files, some are xlsx files

test <- read.csv("data/2017_data/temperature_data_2017/IGH7_temperature_wk2_blk1.txt", header = TRUE, check.names = FALSE)

# Get all .txt files
tempdir <- "data/2017_data/temperature_data_2017"
fileNames <- list.files(path = tempdir, pattern = "*.txt") %>% paste(tempdir, ., sep = "/")
txtfilesList <- lapply(fileNames, function(x) read.csv(x, header = TRUE, check.names = FALSE))

# Get .xlsx files
fileNames <- list.files(path = tempdir, pattern = "*.xlsx") %>% paste(tempdir, ., sep = "/")
xlsxfilesList <- lapply(fileNames, function(x) read.xlsx(x, sheet = 1, cols = 1:3))
