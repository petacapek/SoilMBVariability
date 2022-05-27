GlObjectiveGlucose<-function(x){
  #==========================Extracting sampling time from the dataset
  times <- as.numeric(GlData[GlData$Substrate=="Glucose", c("Time")])
  #==========================Defining initial conditions that are passed to model
  ##Eu at time zero is set to equal to m*Em/Im/yA so it can meet the maintenance costs
  Eu0 = x[4]*x[5]/x[1]/x[3]
  ##X1u is the initial chloroform labile C with Eu subtracted and corrected for incomplete extraction
  X1u0 = GlData$Cmicinit[1]/(x[8] + x[7]*Eu0)#!Cmicinit is in nmols
  y0 <- c(GlData$Sinit[1], 0, 0, 0, Eu0, X1u0)
  #==========================Running simulation
  Yhat <- Glanville2016ODESolv(DEBmodelIso, x, times, y0)
  #==========================Calculating error that is minimized
  ##Measured data
  Y <- data.matrix(GlData[GlData$Substrate=="Glucose", c("CO2", "kec", "Flush")])
  ##Weights
  W <- matrix(rep(apply(Y, 2, sd, na.rm = T), each = length(times)), nrow = dim(Y)[1], ncol = dim(Y)[2])
  ##Error
  return(sum(((Yhat - Y)/W)^2, na.rm=T))
}