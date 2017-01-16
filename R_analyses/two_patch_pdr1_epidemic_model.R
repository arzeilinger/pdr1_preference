#####################################################################################################
#### 2-Patch SECI Movement Model
#####################################################################################################

#### Structuring the parameter estimates, same estimates as single patch model
resistantParams <- parList$R[,1:10]
names(resistantParams) <- paste(names(resistantParams), "r", sep = "")
susceptibleParams <- parList$S[,1:10]
names(susceptibleParams) <- paste(names(susceptibleParams), "s", sep = "")
invariantParams <- parList$R[,11:13]
mrVec <- rep(0.01, nmin) # Emmigration rate from Resistant patch
msVec <- rep(0.01, nmin) # Emmigration rate from Susceptible patch
M <- rep(6000, nmin) # Total number of migrant BGSS
MuVec <- M*0
MvVec <- M*1 # Purcell 1975 estimated that 30% of BGSS in vineyards in early season were infectious
Sr0Vec <- rep(100, nmin)

patchParams <- cbind(resistantParams, susceptibleParams, invariantParams, mrVec, msVec, MuVec, MvVec, Sr0Vec)


#### Initial state variables
# Initial state variables
Ss0 <- 100;
Ur0 <- 199; Vis0 <- 1

# Testing out the model
patchDynamicsOut <- SECIMPatchDynamics(patchParams[5,])
# Plot only hosts because of larger vector numbers
plotDynamics <- patchDynamicsOut %>% dplyr::select(., -time, -starts_with("V"), -starts_with("U"))
matplot(plotDynamics[1:50,], type = "l")

# Testing the simulation function
patchSim <- SECIMPatchSimulations(patchParams[1,])

# Running a bunch of simulations
ti <- proc.time()
patchSim <- apply(patchParams, 1, SECIMPatchSimulations) %>% rbindlist() %>% as.data.frame()
tf <- proc.time()
# Time the simulations took in minutes:
(tf-ti)/60

# Re-structuring the data set, focusing on C, I, and V for the Susceptible Patch
patchSim$Vs <- with(patchSim, Vcs + Vis)
patchSim$Vr <- with(patchSim, Vcr + Vir)

saveRDS(patchSim, file = "output/simulation_results_2_patch_model.rds")


##############################################################################################################################
#### Run 2-patch model over range of Resistant patch sizes
patchParams <- patchParams[1:1000,] %>% dplyr::select(., -Sr0Vec)

Sr0Vec <- seq(1,200,length.out = 10) %>% round()

patchParList <- vector("list", length(Sr0Vec))

for(i in 1:length(Sr0Vec)){
  patchParams$Sr0Vec <- Sr0Vec[i]
  patchParList[[i]] <- patchParams
}

# Understanding what the primary spread term looks like
# sAreaProp <- Ss0/(Sr0Vec + Ss0)
# sM <- sAreaProp*MvVec
# # Number of infectious migrant vectors in susceptible patch as a function of resistant patch area
# plot(Sr0Vec, sM[1:length(Sr0Vec)])

# Running a bunch of simulations
ti <- proc.time()
patchSimList <- lapply(patchParList, function(x) apply(x, 1, SECIMPatchSimulations) %>% rbindlist() %>% as.data.frame())
tf <- proc.time()
# Time the simulations took in minutes:
(tf-ti)/60

#saveRDS(patchSimList, file = "output/simulation_results_2-patch_area.rds")

for(i in 1:length(patchSimList)){
  patchSimList[[i]]$Vs <- with(patchSimList[[i]], Vcs + Vis)
  patchSimList[[i]]$Sr0 <- Sr0Vec[i]
  patchSimList[[i]]$TI <- with(patchSimList[[i]], Cr + Ir + Cs + Is) # Total plants infected over all patches
}
# Convert to data.frame, remove Vc and Vi, and transform to "long" format
patchSimData <- patchSimList %>% rbindlist() %>% as.data.frame() %>% dplyr::select(., Cr, Ir, Cs, Is, TI, Vs, Sr0) %>% 
  gather(., key = state, value = density, Cr, Ir, Cs, Is, TI, Vs)
patchSimData$patchAreaRatio <- patchSimData$Sr0/Ss0


# Summarize simulation results
summaryPatch <- patchSimData %>% group_by(patchAreaRatio, state) %>% summarise(mean = mean(density),
                                                                               median = median(density),
                                                                               sd = sd(density),
                                                                               cil = quantile(density, 0.025),
                                                                               ciu = quantile(density, 0.975),
                                                                               max = max(density))

# Remove Vs
summaryPatch[summaryPatch$state == "Vs",]
summaryPatch <- summaryPatch %>% dplyr::filter(., state != "Vs", state != "Cr", state != "Ir", state != "TI")

#### Plotting with ggplot2
#### Mean infected density of C, I, and V
patchAreaPlot <- ggplot(data=summaryPatch, aes(x=patchAreaRatio, y=mean, group=state, shape=state)) +
  # geom_bar(position=position_dodge(), stat="identity", 
  #          aes(fill=trt)) +
  # geom_hline(aes(yintercept=50), linetype="dashed") +
  geom_line(aes(linetype=state), size=1.25) +
  geom_point(position = position_dodge(width = 0.05), aes(shape=state), size=3.5) +
  #geom_errorbar(aes(ymax=ciu, ymin=cil, width=0.2), position = position_dodge(width = 0.05)) +
  scale_x_continuous(name = "Ratio Resistant patch : Susceptible patch area") +
  #                    breaks = c(3,8,12)) + 
  scale_y_continuous(name = "Proportion infected hosts in susceptible patch") +
  # ylab("% insects on source plant") + 
  #ylim(c(0,100)) +
  # xlab("Genotype") +
  theme_bw(base_size=18) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank()) 

patchAreaPlot





ggsave("results/figures/SECI_patch_area_plot.jpg", plot = patchAreaPlot,
       width = 7, height = 7, units = "in")
