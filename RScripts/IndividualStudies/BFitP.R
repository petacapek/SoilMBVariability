BFitP<-function(x){
  #==========================Extracting sampling time from the dataset
  timesR <- as.numeric(BData[BData$Treatment=="Rhizosphere", c("Time")])
  times <- as.numeric(BData[BData$Treatment!="Rhizosphere", c("Time")])
  timesSim <- seq(0, max(times), length.out = 100)
  #==========================Defining initial conditions that are passed to model
  ##E and X1 at time zero are estimated from the Flush/DNA corrected for respective conversion factors
  B0R = as.numeric(BData[BData$Treatment=="Rhizosphere", c("DNAinit")])[1]/x[6]
  B0 = as.numeric(BData[BData$Treatment!="Rhizosphere", c("DNAinit")])[1]/x[6]
  y0R <- c(as.numeric(BData[BData$Treatment=="Rhizosphere", c("Sinit")])[1], B0R, 0)
  y0 <- c(as.numeric(BData[BData$Treatment!="Rhizosphere", c("Sinit")])[1], B0, 0)
  #==========================Running simulation
  Yhat <- rbind(Blag2014ODESolvP(Pirt, x, timesR, y0R),
                Blag2014ODESolvP(Pirt, x, times, y0))
  Sim <- as.data.frame(rbind(Blag2014ODESolvP(Pirt, x, timesSim, y0R),
                             Blag2014ODESolvP(Pirt, x, timesSim, y0)))
  colnames(Sim) <- c("CO2", "DNA")
  Sim$Time <- rep(timesSim, 2)
  Sim$Treatment <- c(rep("Rhizosphere", times = length(timesSim)),
                     rep("Soil", times = length(timesSim)))
  #==========================Calculating error that is minimized
  ##Measured data
  Y <- data.matrix(BData[, c("CO2", "DNA")])
  ##Weights
  W <- matrix(rep(apply(Y, 2, sd, na.rm = T), each = length(c(times, timesR))), nrow = dim(Y)[1], ncol = dim(Y)[2])
  What <- matrix(rep(apply(Yhat, 2, sd, na.rm = T), each = length(c(times, timesR))), nrow = dim(Y)[1], ncol = dim(Y)[2])
  ##Means
  M <- matrix(rep(apply(Y, 2, mean, na.rm = T), each = length(c(times, timesR))), nrow = dim(Y)[1], ncol = dim(Y)[2])
  Mhat <- matrix(rep(apply(Yhat, 2, mean, na.rm = T), each = length(c(times, timesR))), nrow = dim(Y)[1], ncol = dim(Y)[2])
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