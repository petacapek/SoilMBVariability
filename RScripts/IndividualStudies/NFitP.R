NFitP<-function(x){
  #==========================Extracting sampling time from the dataset
  times <- as.numeric(NData[NData$Treatment=="A", "Time"])
  timesSim <- seq(0, max(times), length.out = 100)
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
  Sim <- as.data.frame(rbind(Nanni1977ODESolvP(Pirt, x, timesSim, y0A),
                             Nanni1977ODESolvP(Pirt, x, timesSim, y0B),
                             Nanni1977ODESolvP(Pirt, x, timesSim, y0C)))
  colnames(Sim) <- c("CO2", "W", "ATP")
  Sim$Time <- rep(timesSim, 3)
  Sim$Treatment <- rep(c("A", "B", "C"), each = length(timesSim))
  #==========================Calculating error that is minimized
  ##Measured data
  Y <- data.matrix(NData[, c("CO2", "W", "ATP")])
  ##Weights
  W <- matrix(rep(apply(Y, 2, sd, na.rm = T), each = length(times)*3), nrow = dim(Y)[1], ncol = dim(Y)[2])
  What <- matrix(rep(apply(Yhat, 2, sd, na.rm = T), each = length(times)*3), nrow = dim(Y)[1], ncol = dim(Y)[2])
  ##Means
  M <- matrix(rep(apply(Y, 2, mean, na.rm = T), each = length(times)*3), nrow = dim(Y)[1], ncol = dim(Y)[2])
  Mhat <- matrix(rep(apply(Yhat, 2, mean, na.rm = T), each = length(times)*3), nrow = dim(Y)[1], ncol = dim(Y)[2])
  ##goodness of fit
  ###R2 adjusted
  R2 = 1 - (sum((Y - Yhat)^2, na.rm = T)/sum((Y - M)^2, na.rm = T))
  R2adj = 1 - ((1 - R2)*((length(Y) - 1)/(length(Y) - 1 - length(x))))
  ###Log-Likelihood
  ll = - sum((Y - Yhat)^2, na.rm = T)/2/sd(Y, na.rm = T)^2
  ###AIC
  AIC = -2*ll + 2*length(x)
  ###normalized F
  Fnorm = sum(((Y - M)/W - (Yhat - Mhat)/What)^2, na.rm = T)
  
  errors = c(R2 = R2, R2adj = R2adj, ll = ll, AIC = AIC, Fnorm = Fnorm, n = length(Y[!is.na(Y)]), p = length(x))
  
  return(list(errors = errors, Simulation = Sim))
}