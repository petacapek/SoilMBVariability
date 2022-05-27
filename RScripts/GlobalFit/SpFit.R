SpFit<-function(x){
  #==========================Time 2 days is reported in the paper only so model is supplemented with continuous
  #variable ending at time 2 (length = 21)
  times <- seq(0, 2, by = 0.1)
  #Model simulations are performed in the for loop because the different initial conditions established for each
  #soil
  ##First soil (i.e. first row of the data frame)
  #==========================Defining initial conditions that are passed to model
  ###Eu at time zero is set to equal to m*Em/Im/yA so it can meet the maintenance costs
  Eu0 = 1e-12
  ###X1u is the initial chloroform labile C with Eu subtracted and corrected for incomplete extraction
  X1u0 = Sp90$Cmicinit[1]/(x[8] + x[7]*Eu0)
  y0 <- c(Sp90$Sinit[1], 0, 0, 0, Eu0, X1u0)
  #==========================Running simulation
  Yhat <- Sparling1990ODESolv(DEBmodelIso, x[1:8], times, y0)[21, ]
  ##The other soils
  for(i in 2:nrow(Sp90)){
    Eu0 = 1e-12
    X1u0 =  Sp90$Cmicinit[i]/(x[8] + x[7]*Eu0)
    y0 <- c(Sp90$Sinit[i], 0, 0, 0, Eu0, X1u0)
    Yhat <- rbind(Yhat, Sparling1990ODESolv(DEBmodelIso, x[1:8], times, y0)[21, ])
  }
  Yhat <- as.data.frame(Yhat)
  colnames(Yhat) <- c("14S", "14CO2", "kec", "14MBC")
  Preds <- melt(Yhat)
  colnames(Preds) <- c("Variable", "Predictions")
  ##Measured data
  Preds$Observations <- melt(Sp90[, c("S", "CO2", "kec", "Flush")])[, 2]
  Preds$Study <- c("Sparling et al. (1988, 1990)")
  
  return(Preds)
}