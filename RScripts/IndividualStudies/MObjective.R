MObjective<-function(x){
  #==========================Extracting sampling time from the dataset
  times <- as.numeric(MData[, c("Time")])
  #==========================Defining initial conditions that are passed to model
  ##E and X1 at time zero are estimated from the Flush/DNA corrected for respective conversion factors
  X10 = MData$DNAinit[1]/x[9]
  E0 = (MData$Cmicinit[1]/X10 - x[8])/x[7]
  y0 <- c(MData$Sinit[1], E0, X10, 0)
  #==========================Running simulation
  Yhat <- Marstorp1999ODESolv(DEBmodel, x, times, y0)
  #==========================Calculating error that is minimized
  ##Measured data
  Y <- data.matrix(MData[, c("S","CO212cumul", "DNA", "Flush")])
  ##Weights
  W <- matrix(rep(apply(Y, 2, sd, na.rm = T), each = length(times)), nrow = dim(Y)[1], ncol = dim(Y)[2])
  ##Error
  return(sum(((Yhat - Y)/W)^2, na.rm=T))
}