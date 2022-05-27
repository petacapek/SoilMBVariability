JFitG<-function(x){
  #==========================Extracting sampling time from the dataset
  times <- as.numeric(JData$Time)
  #==========================Defining initial conditions that are passed to model
  E0 = ((JData$Cmicinit[1]/JData$ATPinit[1])*x[11] - x[8])/(x[7] - (JData$Cmicinit[1]/JData$ATPinit[1])*x[10])
  X10 = JData$Cmicinit[1]/(x[8] + x[7]*E0)
  y0 <- c(as.numeric(JData$Sinit[1]), E0, X10, 0)
  #==========================Running simulation
  Yhat <- as.data.frame(Joerg2002ODESolv(DEBmodel, x[c(1:8, 10, 11)], times, y0))     
  colnames(Yhat) <- c("S", "Flush", "ATP")
  Preds <- melt(Yhat)
  colnames(Preds) <- c("Variable", "Predictions")
  ##Measured data
  Preds$Observations <- melt(JData[, c("S", "Cmic", "ATP")])[, 2]
  Preds$Study <- c("Joergensen and Raubuch (2002)")
  
  return(Preds)
}