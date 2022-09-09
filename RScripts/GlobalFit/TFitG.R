TFitG<-function(x){
  #==========================Extracting sampling time from the dataset
  times <- as.numeric(TData[TData$Treatment=="HG", c("Time")])
  #==========================Defining initial conditions that are passed to model
  E0 = max(1e-12, ((TData$Cmicinit[1]/TData$ATPinit[1])*x[11] - x[8])/(x[7] - (TData$Cmicinit[1]/TData$ATPinit[1])*x[10]))
  X10 = TData$Cmicinit[1]/(x[8] + x[7]*E0)
  y0H <- c(as.numeric(TData[TData$Treatment=="HG", c("Sinit")])[1], E0, X10, 0)
  y0L <- c(as.numeric(TData[TData$Treatment!="HG", c("Sinit")])[1], E0, X10, 0)
  #==========================Running simulation
  Yhat <- as.data.frame(rbind(Tsai1997ODESolv(DEBmodel, c(x[1:8], x[10:11]), times, y0H), 
                Tsai1997ODESolv(DEBmodel, c(x[1:8], x[10:11]), times, y0L)))  
  colnames(Yhat) <- c("CO2", "Flush", "ATP")
  Preds <- melt(Yhat)
  colnames(Preds) <- c("Variable", "Predictions")
  ##Measured data
  Preds$Observations <- melt(TData[, c("CO2cumul", "Cmic", "ATP")])[, 2]
  Preds$Study <- c("Tsai et al. (1997)")
 
  return(Preds)
}