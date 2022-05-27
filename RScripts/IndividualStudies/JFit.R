JFit<-function(x){
  #==========================Extracting sampling time from the dataset
  times <- as.numeric(JData$Time)
  timesSim <- seq(0, max(times), length.out = 100)
  #==========================Defining initial conditions that are passed to model
  E0 = ((JData$Cmicinit[1]/JData$ATPinit[1])*x[10] - x[8])/(x[7] - (JData$Cmicinit[1]/JData$ATPinit[1])*x[9])
  X10 = JData$Cmicinit[1]/(x[8] + x[7]*E0)
  y0 <- c(as.numeric(JData$Sinit[1]), E0, X10, 0)
  #==========================Running simulation
  Yhat <- Joerg2002ODESolv(DEBmodel, x, times, y0)    
  Sim <- as.data.frame(Joerg2002ODESolv(DEBmodel, x, timesSim, y0))
  colnames(Sim) <- c("S", "Flush", "ATP")
  Sim$Time <- timesSim
  #==========================Calculating error that is minimized
  ##Measured data
  Y <- data.matrix(JData[, c("S", "Cmic", "ATP")])
  ##Weights
  W <- matrix(rep(apply(Y, 2, sd, na.rm = T), each = length(times)), nrow = dim(Y)[1], ncol = dim(Y)[2])
  What <- matrix(rep(apply(Yhat, 2, sd, na.rm = T), each = length(times)), nrow = dim(Y)[1], ncol = dim(Y)[2])
  ##Means
  M <- matrix(rep(apply(Y, 2, mean, na.rm = T), each = length(times)), nrow = dim(Y)[1], ncol = dim(Y)[2])
  Mhat <- matrix(rep(apply(Yhat, 2, mean, na.rm = T), each = length(times)), nrow = dim(Y)[1], ncol = dim(Y)[2])
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