TFit<-function(x){
  #==========================Extracting sampling time from the dataset
  times <- as.numeric(TData[TData$Treatment=="HG", c("Time")])
  timesSim <- seq(0, max(times), length.out = 100)
  #==========================Defining initial conditions that are passed to model
  E0 = ((TData$Cmicinit[1]/TData$ATPinit[1])*x[10] - x[8])/(x[7] - (TData$Cmicinit[1]/TData$ATPinit[1])*x[9])
  X10 = TData$Cmicinit[1]/(x[8] + x[7]*E0)
  y0H <- c(as.numeric(TData[TData$Treatment=="HG", c("Sinit")])[1], E0, X10, 0)
  y0L <- c(as.numeric(TData[TData$Treatment!="HG", c("Sinit")])[1], E0, X10, 0)
  #==========================Running simulation
  Yhat <- rbind(Tsai1997ODESolv(DEBmodel, x, times, y0H), 
                Tsai1997ODESolv(DEBmodel, x, times, y0L))    
  Sim <- as.data.frame(rbind(Tsai1997ODESolv(DEBmodel, x, timesSim, y0H), 
               Tsai1997ODESolv(DEBmodel, x, timesSim, y0L)))
  colnames(Sim) <- c("CO2", "Flush", "ATP")
  Sim$Time <- rep(timesSim, 2)
  Sim$Treatment <- rep(c("HG", "LG"), each = length(timesSim))
  #==========================Calculating error that is minimized
  ##Measured data
  Y <- data.matrix(TData[, c("CO2cumul", "Cmic", "ATP")])
  ##Weights
  W <- matrix(rep(apply(Y, 2, sd, na.rm = T), each = length(times)*2), nrow = dim(Y)[1], ncol = dim(Y)[2])
  What <- matrix(rep(apply(Yhat, 2, sd, na.rm = T), each = length(times)*2), nrow = dim(Y)[1], ncol = dim(Y)[2])
  ##Means
  M <- matrix(rep(apply(Y, 2, mean, na.rm = T), each = length(times)*2), nrow = dim(Y)[1], ncol = dim(Y)[2])
  Mhat <- matrix(rep(apply(Yhat, 2, mean, na.rm = T), each = length(times)*2), nrow = dim(Y)[1], ncol = dim(Y)[2])
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
  names(R2all) <- colnames(Sim)[1:3]
  
  errors = c(R2 = R2, R2adj = R2adj, ll = ll, AIC = AIC, Fnorm = Fnorm, n = length(Y[!is.na(Y)]), p = length(x))
  
  return(list(errors = errors, Simulation = Sim, R2all = R2all))
}