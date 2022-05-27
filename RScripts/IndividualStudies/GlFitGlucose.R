GlFitGlucose<-function(x){
  #==========================Extracting sampling time from the dataset
  times <- as.numeric(GlData[GlData$Substrate=="Glucose", c("Time")])
  timesSim <- seq(0, max(times), length.out = 100)
  #==========================Defining initial conditions that are passed to model
  ##Eu at time zero is set to equal to m*Em/Im/yA so it can meet the maintenance costs
  Eu0 = x[4]*x[5]/x[1]/x[3]
  ##X1u is the initial chloroform labile C with Eu subtracted and corrected for incomplete extraction
  X1u0 = GlData$Cmicinit[1]/(x[8] + x[7]*Eu0)#!Cmicinit is in nmols
  y0 <- c(GlData$Sinit[1], 0, 0, 0, Eu0, X1u0)
  #==========================Running simulation
  Yhat <- Glanville2016ODESolv(DEBmodelIso, x, times, y0)
  Sim <- as.data.frame(Glanville2016ODESolv(DEBmodelIso, x, timesSim, y0))
  colnames(Sim) <- c("CO2", "kec", "Flush")
  Sim$Time <- timesSim
  #==========================Calculating error that is minimized
  ##Measured data
  Y <- data.matrix(GlData[GlData$Substrate=="Glucose", c("CO2", "kec", "Flush")])
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
  
  errors = c(R2 = R2, R2adj = R2adj, ll = ll, AIC = AIC, Fnorm = Fnorm, n = length(Y), p = length(x))
  
  return(list(errors = errors, Simulation = Sim))
}