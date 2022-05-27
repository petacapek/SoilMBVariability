NO<-function(x){
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
  Yhat <- rbind(Nanni1977ODESolvATP(DEBmodel, x[c(1:6, 10, 11)], times, y0A),
                Nanni1977ODESolvATP(DEBmodel, x[c(1:6, 10, 11)], times, y0B),
                Nanni1977ODESolvATP(DEBmodel, x[c(1:6, 10, 11)], times, y0C))
  #==========================Calculating error that is minimized
  ##Measured data
  Y <- data.matrix(NData[, c("CO2", "ATP")])
  ##Weights
  W <- matrix(rep(apply(Y, 2, sd, na.rm = T), each = length(times)*3), nrow = dim(Y)[1], ncol = dim(Y)[2])
  ##Error
  return(sum(((Yhat - Y)/W)^2, na.rm=T))
}