GlFit<-function(x){
  #==========================Extracting sampling time from the dataset
  times <- as.numeric(GlData[GlData$Substrate=="Glucose", c("Time")])
  #==========================Defining initial conditions that are passed to model
  ##Eu is estimated
  Eu0 = 1e-12
  ##X1u is the initial chloroform labile C with Eu subtracted and corrected for incomplete extraction
  X1u0 = GlData$Cmicinit[1]/(x[8] + x[7]*Eu0)
  y0 <- c(GlData$Sinit[1]/1000, 0, 0, 0, Eu0, X1u0)
  #==========================Running simulation
  Yhat <- as.data.frame(Glanville2016ODESolv(DEBmodelIso, x[c(1:8)], times, y0))
  Yhat[, 1] <- Yhat[, 1]/1e6
  Yhat[, 3] <- Yhat[, 3]/1e6
  colnames(Yhat) <- c("14CO2", "kec", "14MBC")
  Preds <- melt(Yhat)
  colnames(Preds) <- c("Variable", "Predictions")
  ##Measured data
  Preds$Observations <- c(as.numeric(GlData[GlData$Substrate=="Glucose", c("CO2")])/1e6,
                          as.numeric(GlData[GlData$Substrate=="Glucose", c("kec")]),
                          as.numeric(GlData[GlData$Substrate=="Glucose", c("Flush")])/1e6)
  #Preds$Observations <- melt(((GlData[GlData$Substrate=="Glucose", c("CO2", "kec", "Flush")])))[, 2]
  Preds$Study <- c("Glanville et al. (2016)")
  
  return(Preds)
}