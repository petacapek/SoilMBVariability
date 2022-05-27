ZObjective<-function(x){
  #==========================Extracting sampling time from the dataset
  times <- as.numeric(ZData[, c("Time")])
  #==========================Defining initial conditions that are passed to model
  E0 = x[4]*x[5]/x[1]/x[3]
  X10 = ZData$PLFAinit[1]/(x[8] + x[7]*E0)
  y0 <- c(ZData$Sinit[1], E0, X10, 0)
  #==========================Running simulation
  Yhat <- Ziegler2005ODESolv(DEBmodel, x, times, y0)
  #==========================Calculating error that is minimized
  ##Measured data
  Y <- data.matrix(ZData[, c("S","CO2cumul", "PLFA")])
  ##Weights
  W <- matrix(rep(apply(Y, 2, sd, na.rm = T), each = length(times)), nrow = dim(Y)[1], ncol = dim(Y)[2])
  ##Error
  return(sum(((Yhat - Y)/W)^2, na.rm=T))
}