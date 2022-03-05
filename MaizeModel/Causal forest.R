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
trialData <- na.omit(trialData)

# Y vars
Y <- as.vector(trialData$levyield)

# X vars
firstcol = which(colnames(trialData)=="temp_p1")
lastcol = which(colnames(trialData)=="var_hy")
xvars <- c(firstcol:lastcol)
X <- trialData[xvars]

# Cluster var
sitecode <- as.vector(trialData$site_year_id)

# Treatment var
W <- as.vector(trialData$fert)
Weights <- as.vector(trialData$ps_pro)


########################## Estimate CAUSAL FOREST #####################################

# Get Y.hat in case used
Y.forest = regression_forest(X,Y,sample.weights=Weights,clusters=sitecode,
                             equalize.cluster.weights=FALSE,
                             tune.parameters="all",
                             seed=12345)
Y.hat.all = predict(Y.forest, estimate.variance=TRUE)

# Estimate W.hat using regression forest
forest.W <- regression_forest(X, W, tune.parameters = "all")
W.hat <- predict(forest.W)$predictions

# Estimate causal forest 
cf = causal_forest(X,Y,W,Y.hat=Y.hat.all$predictions,W.hat=Weights,
                   sample.weights=Weights,
                   clusters=sitecode, 
                   equalize.cluster.weights=FALSE,
                   tune.parameters="all",
                   seed=12345,
                   num.trees=500)
varimp = variable_importance(cf)
test_calibration(cf)

# Average treatment effects
#average_treatment_effect(cf,target.sample="all") #CATE
#average_treatment_effect(cf,target.sample="treated") #CATT
average_treatment_effect(cf,target.sample="all") #???
tau.hat=predict(cf, target.sample = "all",estimate.variance=TRUE)
summary(tau.hat$predict)

# Local average treatment effects
#average_late(cf, compliance.score = NULL, subset = NULL)

d <- density(tau.hat$predictions)
plot(d, main="Kernel density of CATE")
polygon(d, col="red", border="blue")

# Describe results (heterogeneity)
plot(X$temp_p1,tau.hat$predictions, main="CATE",xlab="xvar", ylab="fert response", pch=19)
lines(lowess(X$temp_p1,tau.hat$predictions), col=2)


# Export estimation data
exportdata <- data.frame(trialData$id,tau.hat)
write.dta(exportdata, "RandomForest/Results_TrialData_cf.dta")

print(cf)
exportdata <- data.frame(varimp)
write.dta(exportdata, "RandomForest/Varimp_CF.dta")

rm(list=c('exportdata','trialData','X','Y.hat.all','d','tau.hat','Y.forest','forest.W','W.hat','varimp'))
#rm(list = ls(all.names = TRUE))


## Trial sites sim climate dataset

trialDataSimClimate <- read.dta("RandomForest/SimData_trialsim_lev_cf.dta")
trialDataSimClimate <- as.data.frame(trialDataSimClimate)
trialDataSimClimate <- na.omit(trialDataSimClimate)

firstcol = which(colnames(trialDataSimClimate)=="temp_p1")
lastcol = which(colnames(trialDataSimClimate)=="var_hy")
XT_S <- trialDataSimClimate[c(firstcol:lastcol)]

tau.hat.trialsim = predict(cf, newdata=XT_S, estimate.variance=TRUE)
exportdata <- data.frame(trialDataSimClimate$sitecode,trialDataSimClimate$simyr,tau.hat.trialsim)
write.dta(exportdata, "RandomForest/Results_TrialSim_cf.dta")

rm(list=c('exportdata','trialDataSimClimate','XT_S','tau.hat.trialsim'))


## Africa sim datasets

xvarlist <- c("temp_p1","temp_p2","temp_p3","precip_p1","precip_p2","precip_p3","soilcec","soilph","acidity","soilom","soiln","claypct","siltpct","bulkdens","elevm")
hilolist <- c("hi","lo")

for (i in 1:10) {
  inname <- paste0("RandomForest/SimData_robust_cf_seg",i)
  inname <- paste0(inname,".dta")
  outname <- paste0("RandomForest/Results_SimData_cf_seg",i)
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
  
  firstcol = which(colnames(africaDataRobust)=="temp_p1")
  lastcol = which(colnames(africaDataRobust)=="var_hy")
  XA_SR <- africaDataRobust[c(firstcol:lastcol)]
  
  tau.hat.sim.robust = predict(cf, newdata=XA_SR[g1start:g1stop,], estimate.variance=TRUE)
  tau.hat.sim.robust2 = predict(cf, newdata=XA_SR[g2start:g2stop,], estimate.variance=TRUE)
  tau.hat.sim.robust3 = predict(cf, newdata=XA_SR[g3start:g3stop,], estimate.variance=TRUE)
  
  tau.hat.sim.robust = rbind(tau.hat.sim.robust,tau.hat.sim.robust2)
  tau.hat.sim.robust = rbind(tau.hat.sim.robust,tau.hat.sim.robust3)
  
  exportdata <- data.frame(africaDataRobust$cellid,africaDataRobust$iter,tau.hat.sim.robust)
  write.dta(exportdata, outname)
  
  rm(list=c('exportdata','XA_SR','tau.hat.sim.robust','tau.hat.sim.robust2','tau.hat.sim.robust3','africaDataRobust'))
  
  for (varname in xvarlist) { 
    inname <- paste0("RandomForest/SimData_SA_seg",i)
    inname <- paste0(inname,"_")
    inname <- paste0(inname,varname)
    inname <- paste0(inname,"_")
    outname <- paste0("RandomForest/Results_SimData_SA_seg",i)
    outname <- paste0(outname,"_")
    outname <- paste0(outname,varname)
    outname <- paste0(outname,"_")

    for (hilo in hilolist) {
      innamehl <- paste0(inname,hilo)
      innamehl <- paste0(innamehl,".dta")
      outnamehl <- paste0(outname,hilo)
      outnamehl <- paste0(outnamehl,".dta")
      
      africaDataRobust <- read.dta(innamehl)
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
      
      firstcol = which(colnames(africaDataRobust)=="temp_p1")
      lastcol = which(colnames(africaDataRobust)=="var_hy")
      XA_SR <- africaDataRobust[c(firstcol:lastcol)]
      
      tau.hat.sim.robust = predict(cf, newdata=XA_SR[g1start:g1stop,], estimate.variance=TRUE)
      tau.hat.sim.robust2 = predict(cf, newdata=XA_SR[g2start:g2stop,], estimate.variance=TRUE)
      tau.hat.sim.robust3 = predict(cf, newdata=XA_SR[g3start:g3stop,], estimate.variance=TRUE)
      
      tau.hat.sim.robust = rbind(tau.hat.sim.robust,tau.hat.sim.robust2)
      tau.hat.sim.robust = rbind(tau.hat.sim.robust,tau.hat.sim.robust3)
      
      exportdata <- data.frame(africaDataRobust$cellid,africaDataRobust$iter,tau.hat.sim.robust)
      write.dta(exportdata, outnamehl)
      
      rm(list=c('exportdata','XA_SR','tau.hat.sim.robust','tau.hat.sim.robust2','tau.hat.sim.robust3','africaDataRobust'))
      }
    }
  }

