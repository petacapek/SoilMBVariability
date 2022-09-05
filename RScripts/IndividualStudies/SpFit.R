SpFit<-function(x){
  #==========================Extracting sampling time from the dataset
  times <- seq(0, 2, by = 0.1)
  #==========================Defining initial conditions that are passed to model
  ###Eu at time zero is set to equal to m*Em/Im/yA so it can meet the maintenance costs
  Eu0 = x[4]*x[5]/x[1]/x[3]
  ###X1u is the initial chloroform labile C with Eu subtracted and corrected for incomplete extraction
  X1u0 =  Sp90$Cmicinit[1]/(x[8] + x[7]*Eu0)
  y0 <- c(Sp90$Sinit[1], 0, 0, 0, Eu0, X1u0)
  #==========================Running simulation
  Yhat <- Sparling1990ODESolv(DEBmodelIso, x, times, y0)[21, ]
  ##The other soils
  for(i in 2:nrow(Sp90)){
    Eu0 = x[4]*x[5]/x[1]/x[3]
    X1u0 = Sp90$Cmicinit[i]/(x[8] + x[7]*Eu0)
    y0 <- c(Sp90$Sinit[i], 0, 0, 0, Eu0, X1u0)
    Yhat <- rbind(Yhat, Sparling1990ODESolv(DEBmodelIso, x, times, y0)[21, ])
  }
  #==========================Calculating error
  ##Measured data
  Y <- data.matrix(Sp90[, c("S", "CO2", "kec", "Flush")])
  ##Weights
  W <- matrix(rep(apply(Y, 2, sd, na.rm = T), each = dim(Y)[1]), nrow = dim(Y)[1], ncol = dim(Y)[2])
  What <- matrix(rep(apply(Yhat, 2, sd, na.rm = T), each = dim(Y)[1]), nrow = dim(Y)[1], ncol = dim(Y)[2])
  ##Means
  M <- matrix(rep(apply(Y, 2, mean, na.rm = T), each = dim(Y)[1]), nrow = dim(Y)[1], ncol = dim(Y)[2])
  Mhat <- matrix(rep(apply(Yhat, 2, mean, na.rm = T), each = dim(Y)[1]), nrow = dim(Y)[1], ncol = dim(Y)[2])
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
  names(R2all) <- c("S", "CO2", "kec", "Flush")
  
  errors = c(R2 = R2, R2adj = R2adj, ll = ll, AIC = AIC, Fnorm = Fnorm, n = length(Y), p = length(x))
  
  return(list(errors = errors, Simulation = cbind(Y, Yhat), R2all = R2all))
}