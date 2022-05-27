SpObjective<-function(x){
  #==========================Time 2 days is reported in the paper only so model is supplemented with continuous
  #variable ending at time 2 (length = 21)
  times <- seq(0, 2, by = 0.1)
  #Model simulations are performed in the for loop because the different initial conditions established for each
  #soil
  ##First soil (i.e. first row of the data frame)
  #==========================Defining initial conditions that are passed to model
  ###Eu at time zero is set to equal to m*Em/Im/yA so it can meet the maintenance costs
  Eu0 = x[4]*x[5]/x[1]/x[3]
  ###X1u is the initial chloroform labile C with Eu subtracted and corrected for incomplete extraction
  X1u0 = Sp90$Cmicinit[1]/(x[8] + x[7]*Eu0)
  y0 <- c(Sp90$Sinit[1], 0, 0, 0, Eu0, X1u0)
  #==========================Running simulation
  Yhat <- Sparling1990ODESolv(DEBmodelIso, x, times, y0)[21, ]
  ##The other soils
  for(i in 2:nrow(Sp90)){
    Eu0 = x[4]*x[5]/x[1]/x[3]
    X1u0 =  Sp90$Cmicinit[i]/(x[8] + x[7]*Eu0)
    y0 <- c(Sp90$Sinit[i], 0, 0, 0, Eu0, X1u0)
    Yhat <- rbind(Yhat, Sparling1990ODESolv(DEBmodelIso, x, times, y0)[21, ])
  }
  #==========================Calculating error that is minimized
  ##Measured data
  Y <- data.matrix(Sp90[, c("S", "CO2", "kec", "Flush")])
  ##Weights
  W <- matrix(rep(apply(Y, 2, sd, na.rm = T), each = dim(Y)[1]), nrow = dim(Y)[1], ncol = dim(Y)[2])
  ##Error
  return(sum(((Yhat - Y)/W)^2, na.rm=T))
}