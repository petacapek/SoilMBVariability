BFitG<-function(x){
  #==========================Extracting sampling time from the dataset
  timesR <- as.numeric(BData[BData$Treatment=="Rhizosphere", c("Time")])
  times <- as.numeric(BData[BData$Treatment!="Rhizosphere", c("Time")])
  #==========================Defining initial conditions that are passed to model
  ##E and X1 at time zero are estimated from the Flush/DNA corrected for respective conversion factors
  X10R = as.numeric(BData[BData$Treatment=="Rhizosphere", c("DNAinit")])[1]/x[9]
  X10 = as.numeric(BData[BData$Treatment!="Rhizosphere", c("DNAinit")])[1]/x[9]
  E0 = max(x[5]*x[4]/x[1]/x[3], 1e-12)
  y0R <- c(as.numeric(BData[BData$Treatment=="Rhizosphere", c("Sinit")])[1], E0, X10R, 0)
  y0 <- c(as.numeric(BData[BData$Treatment!="Rhizosphere", c("Sinit")])[1], E0, X10, 0)
  #==========================Running simulation
  Yhat <- as.data.frame(rbind(Blag2014ODESolv(DEBmodel, c(x[1:6], x[9]), timesR, y0R),
                Blag2014ODESolv(DEBmodel, c(x[1:6], x[9]), times, y0)))
  colnames(Yhat) <- c("CO2", "DNA")
  Preds <- melt(Yhat)
  colnames(Preds) <- c("Variable", "Predictions")
  Preds$Observations <- melt(BData[, c("CO2", "DNA")])[, 2]
  Preds$Study <- c("Blagodatskaya et al. (2014)")
  
  return(Preds)
}