MFitG<-function(x){
  #==========================Extracting sampling time from the dataset
  times <- as.numeric(MData[, c("Time")])
  #==========================Defining initial conditions that are passed to model
  ##E and X1 at time zero are estimated from the Flush/DNA corrected for respective conversion factors
  X10 = MData$DNAinit[1]/x[9]
  E0 = max(1e-12, (MData$Cmicinit[1]/X10 - x[8])/x[7])
  y0 <- c(MData$Sinit[1], E0, X10, 0)
  #==========================Running simulation
  Yhat <- as.data.frame(Marstorp1999ODESolv(DEBmodel, x[1:9], times, y0))
  colnames(Yhat) <- c("S","CO2", "DNA", "Flush")
  Preds <- melt(Yhat)
  colnames(Preds) <- c("Variable", "Predictions")
  ##Measured data
  Preds$Observations <- melt(MData[, c("S","CO212cumul", "DNA", "Cmic12")])[, 2]
  Preds$Study <- c("Marstorp and Witter (1999)")
 
  ##Error
  return(Preds)
}