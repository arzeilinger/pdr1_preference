#### Functions for analysis of PdR1 experiment data and BGSS colony data

############################################################################
#### Functions for BGSS colony data

#### Function to extract sheet data from multiple Google Sheets
gsExtract <- function(x){
  sheet <- gs_read(gsbgss, ws = x, range = "A3:D20")
  # Make the function wait (time in seconds) until the next request to the Google API
  Sys.sleep(time = 8)
  # remove the header information
  #sheet <- rawsheet[-c(1,2),]
  # Sheet index as column
  sheet$sheetNum <- x
  # Because dates are a mess, make row position a surrogate for number of weeks
  #sheet$weeks <- row.names(sheet)
  return(sheet)
}


#### Function to extract data from .xlsx sheets
excelExtract <- function(x){
  rawsheet <- read.xlsx("data/BGSS Cages 2016.xlsx", sheetIndex = x)
  # remove the header information
  sheet <- rawsheet[-c(1,2),]
  # Sheet index as column
  sheet$sheetNum <- x
  # Because dates are a mess, make row position a surrogate for number of weeks
  sheet$weeks <- row.names(sheet)
  return(sheet)
}


#### Function to return only the most current colony counts from each sheet
currentFunc <- function(cagex){
  currentAdults <- cagex[nrow(cagex), "adults"]
  currentNymphs <- cagex[nrow(cagex), "nymphs"]
  return(c(currentAdults, currentNymphs))
}


#### Function to calculate mortality in each cage using linear regression 
mortality <- function(cagex){
  deaths <- cagex[which(cagex[,"adults"] == max(cagex[,"adults"], na.rm = TRUE)):nrow(cagex),c("adults","weeks")]
  decline <- tryCatch(coef(lm(adults ~ weeks, data = deaths)),
                      error = function(e) NA)
  declineSlope <- as.numeric(decline)[2]
  currentAdults <- cagex[nrow(cagex), "adults"]
  currentNymphs <- cagex[nrow(cagex), "nymphs"]
  return(c(declineSlope,currentAdults, currentNymphs))
}
