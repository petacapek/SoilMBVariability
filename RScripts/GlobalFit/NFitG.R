NFitG<-function(x){
  #==========================Extracting sampling time from the dataset
  times <- as.numeric(NData[NData$Treatment=="A", "Time"])
  #==========================Defining initial conditions that are passed to model
  E0 = 1e-12
  X10A = with(subset(NData, Treatment == "A"), ATPinit/(x[11] + x[10]*E0))[1]
  y0A <- c(as.numeric(NData[NData$Treatment=="A", "Sinit"])[1], E0, X10A, 0)
  X10B = with(subset(NData, Treatment == "B"), ATPinit/(x[11] + x[10]*E0))[1]
  y0B <- c(as.numeric(NData[NData$Treatment=="B", "Sinit"])[1], E0, X10B, 0)
  X10C = with(subset(NData, Treatment == "C"), ATPinit/(x[11] + x[10]*E0))[1]
  y0C <- c(as.numeric(NData[NData$Treatment=="C", "Sinit"])[1], E0, X10C, 0)
  #==========================Running simulation
  Yhat <- as.data.frame(rbind(Nanni1977ODESolvATP(DEBmodel, x[c(1:6, 10, 11)], times, y0A),
                              Nanni1977ODESolvATP(DEBmodel, x[c(1:6, 10, 11)], times, y0B),
                              Nanni1977ODESolvATP(DEBmodel, x[c(1:6, 10, 11)], times, y0C)))
  colnames(Yhat) <- c("CO2", "ATP")
  Preds <- melt(Yhat)
  colnames(Preds) <- c("Variable", "Predictions")
  ##Measured data
  Preds$Observations <- melt(NData[, c("CO2", "ATP")])[, 2]
  Preds$Study <- c("Nannipieri et al. (1977)")
 
  return(Preds)
}