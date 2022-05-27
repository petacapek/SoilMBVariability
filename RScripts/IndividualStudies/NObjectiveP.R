NObjectiveP<-function(x){
  #==========================Extracting sampling time from the dataset
  times <- as.numeric(NData[NData$Treatment=="A", "Time"])
  #==========================Defining initial conditions that are passed to model
  B0A = with(subset(NData, Treatment == "A"), ATPinit/x[7])[1]
  y0A <- c(as.numeric(NData[NData$Treatment=="A", "Sinit"])[1], B0A, 0)
  B0B = with(subset(NData, Treatment == "B"), ATPinit/x[7])[1]
  y0B <- c(as.numeric(NData[NData$Treatment=="B", "Sinit"])[1], B0B, 0)
  B0C = with(subset(NData, Treatment == "C"), ATPinit/x[7])[1]
  y0C <- c(as.numeric(NData[NData$Treatment=="C", "Sinit"])[1], B0C, 0)
  #==========================Running simulation
  Yhat <- rbind(Nanni1977ODESolvP(Pirt, x, times, y0A),
                Nanni1977ODESolvP(Pirt, x, times, y0B),
                Nanni1977ODESolvP(Pirt, x, times, y0C))
  #==========================Calculating error that is minimized
  ##Measured data
  Y <- data.matrix(NData[, c("CO2", "W", "ATP")])
  ##Weights
  W <- matrix(rep(apply(Y, 2, sd, na.rm = T), each = length(times)*3), nrow = dim(Y)[1], ncol = dim(Y)[2])
  ##Error
  return(sum(((Yhat - Y)/W)^2, na.rm=T))
}