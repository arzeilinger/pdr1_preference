#####################################################################################################
#### Cleaning up Xf population data from re-infection experiment to send to Chris Wallis

pdr1DataURL <- gs_url("https://docs.google.com/spreadsheets/d/14uJLfRL6mPrdf4qABeGeip5ZkryXmMKkan3mJHeK13k/edit?usp=sharing",
                      visibility = "private")
ridataGS <- gs_read(pdr1DataURL, ws = "Re-infection_culturing_data")
ridata <- ridataGS
str(ridata)

#### Filter to only genotypes 092 and 094
ridata <- ridata %>% dplyr::filter(genotype == "092" | genotype == "094") %>%
  dplyr::filter(week <= 16) %>%
  dplyr::select(week, genotype, trt, rep, xf_plant_date_cultured, xf_cfu_per_g, notes)
ridata[ridata$week == 5, "xf_cfu_per_g"] <- NA
write.csv(ridata, file = "output/xf_populations_for_chem_analysis.csv", row.names = FALSE)
