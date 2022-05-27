ChObjective<-function(data){
  objective<-function(x){
    #==========================Extracting sampling time from the dataset
    times <- as.numeric(data$Time)
    #==========================Defining initial conditions that are passed to model
    ##E and X1 at time zero are estimated from the Flush/DNA corrected for respective conversion factors
    X10 = data$DNAinit[1]/x[9]
    E0 = (data$Cmicinit[1]/X10 - x[8])/x[7]
    y0 <- c(data$Sinit[1], E0, X10, 0)
    #==========================Running simulation
    Yhat <- Chen2019ODESolv(DEBmodel, x, times, y0)
    #==========================Calculating error that is minimized
    ##Measured data
    Y <- data.matrix(data[, c("CO2", "DNA", "Cmic")])
    ##Weights
    W <- matrix(rep(apply(Y, 2, sd, na.rm = T), each = length(times)), nrow = dim(Y)[1], ncol = dim(Y)[2])
    ##Error
    return(sum(((Yhat - Y)/W)^2, na.rm=T))
  }
  goodness<-function(x){
    #==========================Extracting sampling time from the dataset
    times <- as.numeric(data$Time)
    timesSim <- seq(0, max(times), length.out = 100)
    #==========================Defining initial conditions that are passed to model
    ##E and X1 at time zero are estimated from the Flush/DNA corrected for respective conversion factors
    X10 = data$DNAinit[1]/x[9]
    E0 = (data$Cmicinit[1]/X10 - x[8])/x[7]
    y0 <- c(data$Sinit[1], E0, X10, 0)
    #==========================Running simulation
    Yhat <- Chen2019ODESolv(DEBmodel, x, times, y0)
    Sim <- as.data.frame(Chen2019ODESolv(DEBmodel, x, timesSim, y0))
    colnames(Sim) <- c("CO2", "DNA", "Flush")
    Sim$Time <- timesSim
    #==========================Calculating error that is minimized
    ##Measured data
    Y <- data.matrix(data[, c("CO2", "DNA", "Cmic")])
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
  ##First guess by MCMC 
  CP <- modMCMC(objective, p = c(Parms[,1], 0.1), 
                lower = c(Parms[, 2], 0),
                upper = c(Parms[, 3], 1), niter = 30000)
  ##Estimate
  ParmsChen2019 <- abc_optim(fn = objective, 
                                 par = as.numeric(summary(CP)[c("mean"), ]), 
                                 lb = as.numeric(summary(CP)[c("min"), ]), 
                                 ub = as.numeric(summary(CP)[c("max"), ])) 
  ChFit <- goodness(ParmsChen2019$par)
  
  return(list(Pars = ParmsChen2019$par, Fit = ChFit))
}