GlO<-function(x){
  #==========================Extracting sampling time from the dataset
  times <- as.numeric(GlData[GlData$Substrate=="Glucose", c("Time")])
  #==========================Defining initial conditions that are passed to model
  ##Eu is estimated
  Eu0 = max(x[5]*x[4]/x[1]/x[3], 1e-12)
  ##X1u is the initial chloroform labile C with Eu subtracted and corrected for incomplete extraction
  X1u0 = GlData$Cmicinit[1]/(x[8] + x[7]*Eu0)
  y0 <- c(GlData$Sinit[1]/1000, 0, 0, 0, Eu0, X1u0)
  #==========================Running simulation
  Yhat <- Glanville2016ODESolv(DEBmodelIso, x[1:8], times, y0)
  Yhat[, 1] <- Yhat[, 1]/1e6
  Yhat[, 3] <- Yhat[, 3]/1e6
  #==========================Calculating error that is minimized
  ##Measured data
  Y <- data.matrix(GlData[GlData$Substrate=="Glucose", c("CO2", "kec", "Flush")])
  Y[, 1] <- Y[, 1]/1e6
  Y[, 3] <- Y[, 3]/1e6
  ##Weights
  W <- matrix(rep(apply(Y, 2, sd, na.rm = T), each = length(times)), nrow = dim(Y)[1], ncol = dim(Y)[2])
  ##Error
  return(sum(((Yhat[, 2] - Y[, 2])/W[, 2])^2, na.rm=T))
}