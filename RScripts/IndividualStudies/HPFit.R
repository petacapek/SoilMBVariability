HPFit<-function(x){
  #==========================Extracting sampling time from the dataset
  times <- as.numeric(HPData[, c("Time")])
  timesSim <- seq(0, max(times), length.out = 100)
  #==========================Defining initial conditions that are passed to model
  ##Eu at time zero is set to equal to m*Em/Im/yA so it can meet the maintenance costs
  Eu0 = 1e-25
  ##X1u is the initial chloroform labile C with Eu subtracted and corrected for incomplete extraction
  X1u0 = HPData$Cmicinit[1]/(x[8])
  y0 <- c(HPData$Sinit[1], 0, 0, 0, Eu0, X1u0)
  #==========================Running simulation
  Yhat <- Sparling1990ODESolv(DEBmodelIso, x, times, y0)
  Sim <- as.data.frame(Sparling1990ODESolv(DEBmodelIso, x, timesSim, y0))
  colnames(Sim) <- c("S", "CO2", "kec", "Flush")
  Sim$Time <- timesSim
  #==========================Calculating error that is minimized
  ##Measured data
  Y <- data.matrix(HPData[, c("S", "CO2", "kec", "Flush")])
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
  
  ###R2 for individual variables
  R2all <- numeric()
  for(i in 1:dim(Y)[2]){
    R2all <- append(R2all, 1 - (sum((Y[, i] - Yhat[, i])^2, na.rm = T)/sum((Y[, i] - M[, i])^2, na.rm = T)))
  }
  names(R2all) <- colnames(Sim)[-5]
  
  errors = c(R2 = R2, R2adj = R2adj, ll = ll, AIC = AIC, Fnorm = Fnorm, n = length(Y), p = length(x))
  
  return(list(errors = errors, Simulation = Sim, R2all = R2all))
}