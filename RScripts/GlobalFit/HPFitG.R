HPFitG<-function(x){
  #==========================Extracting sampling time from the dataset
  times <- as.numeric(HPData[, c("Time")])
  #==========================Defining initial conditions that are passed to model
  ##Eu at time zero is set to equal to m*Em/Im/yA so it can meet the maintenance costs
  Eu0 = max(x[5]*x[4]/x[1]/x[3], 1e-12)
  ##X1u is the initial chloroform labile C with Eu subtracted and corrected for incomplete extraction
  X1u0 = HPData$Cmicinit[1]/(x[8] + x[7]*Eu0)
  y0 <- c(HPData$Sinit[1], 0, 0, 0, Eu0, X1u0)
  #==========================Running simulation
  Yhat <- as.data.frame(Sparling1990ODESolv(DEBmodelIso, x[1:8], times, y0))
  colnames(Yhat) <- c("S","14CO2", "kec", "14MBC")
  Preds <- melt(Yhat)
  colnames(Preds) <- c("Variable", "Predictions")
  ##Measured data
  Preds$Observations <- melt(HPData[, c("S","CO2", "kec", "Flush")])[, 2]
  Preds$Study <- c("Santruckova et al. (2004)")
  
  return(Preds)
}