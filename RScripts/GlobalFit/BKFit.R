BKFit<-function(x){
  #==========================Extracting sampling time from the dataset
  times <- as.numeric(BKData[BKData$Treatment=="LC", c("Time")])
  #==========================Defining initial conditions that are passed to model
  ##Eu at time zero approaches zero
  Eu0 = 1e-12
  ##X1u is the initial chloroform labile C with Eu subtracted and corrected for incomplete extraction
  X1u0 = BKData$Cmicinit[1]/(x[8])
  y0HC <- c(as.numeric(BKData[BKData$Treatment=="HC", c("Sinit")])[1], 0, 0, 0, Eu0, X1u0)
  y0HN <- c(as.numeric(BKData[BKData$Treatment=="HCHN", c("Sinit")])[1], 0, 0, 0, Eu0, X1u0)
  #==========================Running simulation
  Yhat <- as.data.frame(rbind(Bremer1990ODESolv(DEBmodelIso, x[1:8], times, y0HC),
                              Bremer1990ODESolv(DEBmodelIso, x[1:8], times, y0HN)))
  colnames(Yhat) <- c("14S","14CO2", "kec", "14MBC")
  Preds <- melt(Yhat)
  colnames(Preds) <- c("Variable", "Predictions")
  ##Measured data
  Preds$Observations <- melt(BKData[BKData$Treatment!="LC", c("S","CO214cumul", "kec", "Cmic14")])[, 2]
  Preds$Study <- c("Bremer and van Kessel (1990)")
  
  return(Preds)
}