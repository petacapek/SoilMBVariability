NObjective<-function(x){
  #==========================Extracting sampling time from the dataset
  times <- as.numeric(NData[NData$Treatment=="A", "Time"])
  #==========================Defining initial conditions that are passed to model
  E0A = with(subset(NData, Treatment == "A"), ((Winit/ATPinit)*x[10] - x[8])/(x[7] - (Winit/ATPinit)*x[9]))[1] #ATP is in nmols
  X10A = with(subset(NData, Treatment == "A"), ATPinit/(x[10] + x[9]*E0A))[1]
  y0A <- c(as.numeric(NData[NData$Treatment=="A", "Sinit"])[1], E0A, X10A, 0)
  E0B = with(subset(NData, Treatment == "B"), ((Winit/ATPinit)*x[10] - x[8])/(x[7] - (Winit/ATPinit)*x[9]))[1]
  X10B = with(subset(NData, Treatment == "B"), ATPinit/(x[10] + x[9]*E0B))[1]
  y0B <- c(as.numeric(NData[NData$Treatment=="B", "Sinit"])[1], E0B, X10B, 0)
  E0C = with(subset(NData, Treatment == "C"), ((Winit/ATPinit)*x[10] - x[8])/(x[7] - (Winit/ATPinit)*x[9]))[1]
  X10C = with(subset(NData, Treatment == "C"), ATPinit/(x[10] + x[9]*E0C))[1]
  y0C <- c(as.numeric(NData[NData$Treatment=="C", "Sinit"])[1], E0C, X10C, 0)
  #==========================Running simulation
  Yhat <- rbind(Nanni1977ODESolv(DEBmodel, x, times, y0A),
                Nanni1977ODESolv(DEBmodel, x, times, y0B),
                Nanni1977ODESolv(DEBmodel, x, times, y0C))
  #==========================Calculating error that is minimized
  ##Measured data
  Y <- data.matrix(NData[, c("CO2", "W", "ATP")])
  ##Weights
  W <- matrix(rep(apply(Y, 2, sd, na.rm = T), each = length(times)*3), nrow = dim(Y)[1], ncol = dim(Y)[2])
  ##Error
  return(sum(((Yhat - Y)/W)^2, na.rm=T))
}