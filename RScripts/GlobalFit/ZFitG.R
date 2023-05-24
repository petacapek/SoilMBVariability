ZFitG<-function(x){
  #==========================Extracting sampling time from the dataset
  times <- as.numeric(ZData[, c("Time")])
  #==========================Defining initial conditions that are passed to model
  ##Eu at time zero is set to equal to m*Em/Im/yA so it can meet the maintenance costs
  Eu0 = 1e-12
  ##X1u is the initial chloroform labile C with Eu subtracted and corrected for incomplete extraction
  X1u0 = ZData$PLFAinit[1]/x[12]
  y0 <- c(ZData$Sinit[1], 0, 0, 0, Eu0, X1u0)
  #==========================Running simulation
  Yhat <- as.data.frame(Ziegler2005ODESolv(DEBmodelIso, x[c(1:6, 12)], times, y0))
  colnames(Yhat) <- c("S","CO2", "PLFA")
  Preds <- melt(Yhat)
  colnames(Preds) <- c("Variable", "Predictions")
  ##Measured data
  Preds$Observations <- melt(ZData[, c("S","CO2cumul", "PLFA")])[, 2]
  Preds$Study <- c("Ziegler et al. (2005)")
  
  ##Error
  return(Preds)
  
}