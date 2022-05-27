MObjectiveP<-function(x){
  #==========================Extracting sampling time from the dataset
  times <- as.numeric(MData[, c("Time")])
  #==========================Defining initial conditions that are passed to model
  ##B at time zero are estimated from the DNA and respective conversion factor
  B0 = MData$DNAinit[1]/x[7]
  y0 <- c(MData$Sinit[1], B0, 0)
  #==========================Running simulation
  Yhat <- Marstorp1999ODESolvP(Pirt, x, times, y0)
  #==========================Calculating error that is minimized
  ##Measured data
  Y <- data.matrix(MData[, c("S","CO212cumul", "DNA", "Flush")])
  ##Weights
  W <- matrix(rep(apply(Y, 2, sd, na.rm = T), each = length(times)), nrow = dim(Y)[1], ncol = dim(Y)[2])
  ##Error
  return(sum(((Yhat - Y)/W)^2, na.rm=T))
}