TObjectiveM<-function(x){
  #==========================Extracting sampling time from the dataset
  times <- as.numeric(TData[TData$Treatment=="HG", c("Time")])
  #==========================Defining initial conditions that are passed to model
  B0 = TData$Cmicinit[1]/x[5]
  y0H <- c(as.numeric(TData[TData$Treatment=="HG", c("Sinit")])[1], B0, 0)
  y0L <- c(as.numeric(TData[TData$Treatment!="HG", c("Sinit")])[1], B0, 0)
  #==========================Running simulation
  Yhat <- rbind(Tsai1997ODESolvM(Monod, x, times, y0H), 
                Tsai1997ODESolvM(Monod, x, times, y0L))       
  #==========================Calculating error that is minimized
  ##Measured data
  Y <- data.matrix(TData[, c("CO2cumul", "Cmic", "ATP")])
  ##Weights
  W <- matrix(rep(apply(Y, 2, sd, na.rm = T), each = length(times)*2), nrow = dim(Y)[1], ncol = dim(Y)[2])
  ##Error
  return(sum(((Yhat - Y)/W)^2, na.rm=T))
}