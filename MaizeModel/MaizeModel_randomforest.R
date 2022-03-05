necessaryPackages <- c("readxl", "tidyverse", "ggplot2", "grf", "sm", "foreign", "car")
new.packages <- necessaryPackages[
  !(necessaryPackages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(necessaryPackages, require, character.only = TRUE)

citation(package = "grf")
memory.size(max=TRUE)
gc()


########################## LOAD DATA #####################################

trialData <- read.dta("RandomForest/TrialData_act_lev.dta")
trialData <- as.data.frame(trialData)

trialDataF0 <- read.dta("RandomForest/TrialData_f0_lev.dta")
trialDataF0 <- as.data.frame(trialDataF0)

trialDataF1 <- read.dta("RandomForest/TrialData_f1_lev.dta")
trialDataF1 <- as.data.frame(trialDataF1)


# Y vars
Y <- as.vector(trialData$levyield)

# X vars
firstcol = which(colnames(trialData)=="fert")
lastcol = which(colnames(trialData)=="var_hy")
X <- trialData[c(firstcol:lastcol)]

# Cluster var
sitecode <- as.vector(trialData$site_year_id)

# Treatment var
W <- as.vector(trialData$fert)
Weights <- as.vector(trialData$ps_pro)

# Trial predictions
firstcol = which(colnames(trialDataF0)=="fert")
lastcol = which(colnames(trialDataF0)=="var_hy")
XT_F0 <- trialDataF0[c(firstcol:lastcol)]
XT_F1 <- trialDataF1[c(firstcol:lastcol)]


########################## Estimate RANDOM FOREST #####################################

Y.forest = regression_forest(X,Y,sample.weights=Weights,clusters=sitecode,
                             equalize.cluster.weights=FALSE,
                             tune.parameters="all",
                             seed=12345)
print(Y.forest)
varimp = variable_importance(Y.forest)

exportdata <- data.frame(varimp)
write.dta(exportdata, "RandomForest/Varimp_RF.dta")


## Predictions

Y.hat = predict(Y.forest, estimate.variance=TRUE)
RMSE <- sqrt(sum((Y.hat$predictions - Y)^2)/length(Y.hat$predictions))
print(RMSE) #2.22 [w wortmann and elev vs 2.23 before vs 2.27225 RMSE OOB with interactions]
print(RMSE/mean(Y)) #0.5531 with best spec above
print(sd(Y)) #2.625 St Dev Y

plot(Y.hat$predictions,Y, main="Pred vs Act",xlab="Y pred", ylab="Y act", pch=19)
lines(lowess(Y.hat$predictions,Y), col=2)

Y.hat.f0 = predict(Y.forest, newdata=XT_F0, estimate.variance=TRUE)
Y.hat.f1 = predict(Y.forest, newdata=XT_F1, estimate.variance=TRUE)
fertdif = Y.hat.f1$predictions-Y.hat.f0$predictions


## Describe fert differences

d <- density(fertdif)
plot(d, main="Kernel density of fert dif")
polygon(d, col="red", border="blue")

exportdata <- data.frame(trialData$id,Y,Y.hat,Y.hat.f0,Y.hat.f1,fertdif,Weights)
write.dta(exportdata, "RandomForest/Results_TrialData_rf.dta")

rm(list=c('exportdata','trialData','trialDataF0','trialDataF1','X','XT_F0','XT_F1','Y.hat','Y.hat.f0','Y.hat.f1'))



########################## Simulate RANDOM FOREST #####################################

## Trial sites sim climate dataset

trialDataSimClimate <- read.dta("RandomForest/SimData_trialsim_lev.dta")
trialDataSimClimate <- as.data.frame(trialDataSimClimate)
trialDataSimClimate <- na.omit(trialDataSimClimate)

firstcol = which(colnames(trialDataSimClimate)=="fert")
lastcol = which(colnames(trialDataSimClimate)=="var_hy")
XT_S <- trialDataSimClimate[c(firstcol:lastcol)]

Y.hat.trialsim = predict(Y.forest, newdata=XT_S, estimate.variance=TRUE)
exportdata <- data.frame(trialDataSimClimate$sitecode,trialDataSimClimate$simyr,trialDataSimClimate$fert,Y.hat.trialsim)
write.dta(exportdata, "RandomForest/Results_TrialSim_rf.dta")

rm(list=c('exportdata','trialDataSimClimate','XT_S','Y.hat.trialsim'))


## Africa sim datasets

for (i in 1:10) {
  inname <- paste0("RandomForest/SimData_robust_seg",i)
  inname <- paste0(inname,".dta")
  outname <- paste0("RandomForest/Results_SimData_rf_seg",i)
  outname <- paste0(outname,".dta")
  
  africaDataRobust <- read.dta(inname)
  africaDataRobust <- as.data.frame(africaDataRobust)
  africaDataRobust <- na.omit(africaDataRobust)
  
  N = nrow(africaDataRobust)
  split = floor(N/3)
  g1start = 1
  g1stop = split
  g2start = split + 1
  g2stop = g2start + split
  g3start = g2stop + 1
  g3stop = N
  
  firstcol = which(colnames(africaDataRobust)=="fert")
  lastcol = which(colnames(africaDataRobust)=="var_hy")

  XA_SR <- africaDataRobust[c(firstcol:lastcol)]
  
  Y.hat.sim.robust = predict(Y.forest, newdata=XA_SR[g1start:g1stop,], estimate.variance=TRUE)
  Y.hat.sim.robust2 = predict(Y.forest, newdata=XA_SR[g2start:g2stop,], estimate.variance=TRUE)
  Y.hat.sim.robust3 = predict(Y.forest, newdata=XA_SR[g3start:g3stop,], estimate.variance=TRUE)
  
  Y.hat.sim.robust = rbind(Y.hat.sim.robust,Y.hat.sim.robust2)
  Y.hat.sim.robust = rbind(Y.hat.sim.robust,Y.hat.sim.robust3)
  
  exportdata <- data.frame(africaDataRobust$cellid,africaDataRobust$iter,africaDataRobust$fert,Y.hat.sim.robust)
  write.dta(exportdata, outname)
  
  rm(list=c('exportdata','XA_SR','Y.hat.sim.robust','Y.hat.sim.robust2','Y.hat.sim.robust3','africaDataRobust'))
  }

