UObjective<-function(data){
  objective<-function(x){
    #==========================Extracting sampling time from the dataset
    times <- as.numeric(data$Time)
    #==========================Defining initial conditions that are passed to model
    ##E and X1 at time zero are estimated
    #E0 = ((data$Flush[1]/data$B[1]) - x[8]*x[9])/(x[7]*x[9] - (data$Flush[1]/data$B[1]))
    E0 = x[4]*x[5]/x[1]/x[3]
    X10 = data$B[1]/(1 + E0)
    y0 <- c(data$S0[1], E0, X10, 0)
    #==========================Running simulation
    Yhat <- HasanODESolv(DEBmodel, x, times, y0)
    #==========================Calculating error that is minimized
    ##Measured data
    Y <- data.matrix(data[, c("CO2", "B", "Walls")])
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
    #E0 = ((data$Flush[1]/data$B[1]) - x[8]*x[9])/(x[7]*x[9] - (data$Flush[1]/data$B[1]))
    E0 = x[4]*x[5]/x[1]/x[3]
    X10 = data$B[1]/(1 + E0)
    y0 <- c(data$S0[1], E0, X10, 0)
    #==========================Running simulation
    Yhat <- HasanODESolv(DEBmodel, x, times, y0)
    Sim <- as.data.frame(HasanODESolv(DEBmodel, x, timesSim, y0))
    colnames(Sim) <- c("CO2", "B", "Walls")
    Sim$Time <- timesSim
    #==========================Calculating error that is minimized
    ##Measured data
    Y <- data.matrix(data[, c("CO2", "B", "Walls")])
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
  CP <- modMCMC(objective, p = c(Parms[,1]), 
                lower = c(Parms[, 2]),
                upper = c(Parms[, 3]), niter = 30000)
  ##Estimate
  ParmsHasan <- abc_optim(fn = objective, 
                                 par = as.numeric(summary(CP)[c("mean"), ]), 
                                 lb = as.numeric(summary(CP)[c("min"), ]), 
                                 ub = as.numeric(summary(CP)[c("max"), ])) 
  UFit <- goodness(ParmsHasan$par)
  
  return(list(Pars = ParmsHasan$par, Fit = UFit))
}