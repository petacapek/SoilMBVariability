JObjectiveP<-function(x){
  #==========================Extracting sampling time from the dataset
  times <- as.numeric(JData$Time)
  #==========================Defining initial conditions that are passed to model
  B0 = JData$Cmicinit[1]/x[6]
  y0 <- c(as.numeric(JData$Sinit[1]), B0, 0)
  #==========================Running simulation
  Yhat <- Joerg2002ODESolvP(Pirt, x, times, y0)     
  #==========================Calculating error that is minimized
  ##Measured data
  Y <- data.matrix(JData[, c("S", "Cmic", "ATP")])
  ##Weights
  W <- matrix(rep(apply(Y, 2, sd, na.rm = T), each = length(times)), nrow = dim(Y)[1], ncol = dim(Y)[2])
  ##Error
  return(sum(((Yhat - Y)/W)^2, na.rm=T))
}