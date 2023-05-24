BO<-function(x){
  #==========================Extracting sampling time from the dataset
  timesR <- as.numeric(BData[BData$Treatment=="Rhizosphere", c("Time")])
  times <- as.numeric(BData[BData$Treatment!="Rhizosphere", c("Time")])
  #==========================Defining initial conditions that are passed to model
  ##E and X1 at time zero are estimated from the Flush/DNA corrected for respective conversion factors
  X10R = as.numeric(BData[BData$Treatment=="Rhizosphere", c("DNAinit")])[1]/x[9]
  X10 = as.numeric(BData[BData$Treatment!="Rhizosphere", c("DNAinit")])[1]/x[9]
  E0 = 1e-12
  y0R <- c(as.numeric(BData[BData$Treatment=="Rhizosphere", c("Sinit")])[1], E0, X10R, 0)
  y0 <- c(as.numeric(BData[BData$Treatment!="Rhizosphere", c("Sinit")])[1], E0, X10, 0)
  #==========================Running simulation
  Yhat <- rbind(Blag2014ODESolv(DEBmodel, c(x[1:6], x[9]), timesR, y0R),
                Blag2014ODESolv(DEBmodel, c(x[1:6], x[9]), times, y0))
  #==========================Calculating error that is minimized
  ##Measured data
  Y <- data.matrix(BData[, c("CO2", "DNA")])
  ##Weights
  W <- matrix(rep(apply(Y, 2, sd, na.rm = T), each = length(c(times, timesR))), nrow = dim(Y)[1], ncol = dim(Y)[2])
  ##Error
  return(sum(((Yhat - Y)/W)^2, na.rm=T))
}