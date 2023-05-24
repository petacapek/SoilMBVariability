HPO<-function(x){
  #==========================Extracting sampling time from the dataset
  times <- as.numeric(HPData[, c("Time")])
  #==========================Defining initial conditions that are passed to model
  ##Eu at time zero is set to equal to m*Em/Im/yA so it can meet the maintenance costs
  Eu0 = 1e-12
  ##X1u is the initial chloroform labile C with Eu subtracted and corrected for incomplete extraction
  X1u0 = HPData$Cmicinit[1]/(x[8] + x[7]*Eu0)
  y0 <- c(HPData$Sinit[1], 0, 0, 0, Eu0, X1u0)
  #==========================Running simulation
  Yhat <- Sparling1990ODESolv(DEBmodelIso, x[1:8], times, y0)
  #==========================Calculating error that is minimized
  ##Measured data
  Y <- data.matrix(HPData[, c("S","CO2", "kec", "Flush")])
  ##Weights
  W <- matrix(rep(apply(Y, 2, sd, na.rm = T), each = length(times)), nrow = dim(Y)[1], ncol = dim(Y)[2])
  ##Error
  return(sum(((Yhat - Y)/W)^2, na.rm=T))
}