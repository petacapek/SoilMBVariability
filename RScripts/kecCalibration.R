#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~kec coefficient calibration~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~applying DEB theory~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Libraries
library(ggplot2)
library(dplyr)
library(reshape)
library(minpack.lm)
library(bbmle)
library(deSolve)
library(FME)
library(DEoptim)
library(ABCoptim)
library(rcompanion)
library(optimx)
library(reticulate)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#ggplot theme
theme_min<-theme(axis.text.x=element_text(vjust=0.2, size=18, colour="black"),
                 axis.text.y=element_text(hjust=0.2, size=18, colour="black"),
                 axis.title=element_text(size=18, colour="black"),
                 axis.line=element_line(size=0.5, colour="black"),
                 strip.text=element_text(size=18, face="bold"),
                 axis.ticks=element_line(size=1, colour="black"),
                 #axis.ticks.length=unit(-0.05, "cm"),
                 panel.background=element_rect(colour="black", fill="white"),
                 panel.grid=element_line(linetype=0),
                 legend.text=element_text(size=14, colour="black"),
                 legend.title=element_text(size=14, colour="black"),
                 legend.position=c("right"),
                 legend.key.size=unit(1, "cm"),
                 strip.background=element_rect(fill="grey98", colour="black"),
                 legend.key=element_rect(fill="white", size=1.2),
                 legend.spacing=unit(0.5, "cm"),
                 plot.title=element_text(size=18, face="bold", hjust=-0.05))
#===========================
# All estimated parameters are stored here
parsAll <- data.frame(Study = character(), Treatment = character(), 
                      Im = numeric(), Km = numeric(), yA = numeric(),
                      Em = numeric(), m = numeric(), g = numeric(), ne = numeric(), nX1 = numeric(),
                      iX1 = numeric(), te = numeric(), tX1 = numeric(), re = numeric(), rX1 = numeric(), 
                      pe = numeric(), pX1 = numeric(), le = numeric(), lX1 = numeric())
#===========================
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Study of Glanville et al. (2016)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Data
GlData = read.csv("../SoilMBVariabilityData/Glanville2016.csv")
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
##All C pools are in pmol(C)/g(DW) except of initial chloroform labile C, which is in nmol(C)/g(DW)

#Visualizing the data
##kec
ggplot(GlData, aes(Time, kec)) + geom_point(cex=6, pch=21, aes(fill = Substrate)) +
  theme_min + ylab(expression(paste(italic(k[ec])~factor))) + xlab("Time (days)") + 
  theme(legend.position = c(0.8, 0.6), legend.title = element_blank()) + 
  scale_fill_manual(values = c("white", "grey")) + scale_y_continuous(limits = c(0.2, 0.43)) +
  geom_hline(yintercept = 0.38, lwd = 1, lty = 2)
##14C Flush
ggplot(GlData, aes(Time, Flush)) + geom_point(cex=6, pch=21, aes(fill = Substrate)) +
  theme_min + ylab(expression(paste("Flush (pmo", l^{14},C~g(DW)^{-1}, ")"))) + xlab("Time (days)") + 
  theme(legend.position = c(0.8, 0.8), legend.title = element_blank()) + 
  scale_fill_manual(values = c("white", "grey"))
##14C CO2
ggplot(GlData, aes(Time, CO2)) + geom_point(cex=6, pch=21, aes(fill = Substrate)) +
  theme_min + ylab(expression(paste(CO[2], " (pmo", l^{14},C~g(DW)^{-1}, ")"))) + xlab("Time (days)") + 
  theme(legend.position = c(0.8, 0.4), legend.title = element_blank()) + 
  scale_fill_manual(values = c("white", "grey"))

#Model calibration
##Reading respective python scripts
source_python("../PythonScripts/IndividualStudies/DEBmodelIso.py")
source_python("../PythonScripts/IndividualStudies/Glanville2016ODESolv.py")
#=====================TEST FOR GLUCOSE not run
# #time
# times<-as.numeric(GlData[GlData$Substrate=="Glucose", c("Time")])
# #parameters
# p<-c(Im = 1, Km = 50, yA=0.8, Em=1, m=1e-5, g=0.3, ne=0.8, nX1=0.4)
# #initial states
# y0<-c(10, 0, 0, 0, 1, 3)
# testout<-Glanville2016ODESolv(DEBmodelIso, p, times, y0)
#=====================GLUCOSE
#Objective function
source("IndividualStudies/GlObjectiveGlucose.R")
##Optimization
#==================================================#
##Parameters (Initial guess, lower and upper bound) 
Im = c(1, 1e-2, 20)
Km = c(25, 0.1, 3000)
yA = c(0.9, 0, 1)
Em = c(1, 1e-3, 1e3)
m = c(1e-3, 1e-8, 1)
g = c(0.3, 0.01, 10)
ne = c(0.8, 0, 1)
nX1 = c(0.3, 0, 1)

Parms = rbind(Im, Km, yA, Em, m, g, ne, nX1)
#==================================================#
##First guess by MCMC 
GlP1 <- modMCMC(GlObjectiveGlucose, p = Parms[,1], lower = Parms[, 2], upper = Parms[, 3], niter = 30000)
summary(GlP1)
##Estimate
ParmsGlanvilleGlucose <- abc_optim(fn = GlObjectiveGlucose, 
                                   par = as.numeric(summary(GlP1)[c("mean"), ]), 
                                   lb = as.numeric(summary(GlP1)[c("min"), ]), 
                                   ub = as.numeric(summary(GlP1)[c("max"), ]))
ParmsGlanvilleGlucose$par
##Uncertainty
GlPsd <- modMCMC(GlObjectiveGlucose, p = ParmsGlanvilleGlucose$par, lower = Parms[1:8, 2], 
                 upper = Parms[1:8, 3], niter = 5000)
summary(GlPsd)
#==================================================#
#Goodness of fit and simulations
source("IndividualStudies/GlFitGlucose.R")
SimGlanvilleGlucose <- GlFitGlucose(ParmsGlanvilleGlucose$par)
SimGlanvilleGlucose$errors
SimGlanvilleGlucose$R2all
#=====================Alanin
#Objective function
source("IndividualStudies/GlObjectiveAla.R")
##Optimization
##First guess by MCMC 
GlP2 <- modMCMC(GlObjectiveAla, p = Parms[,1], lower = Parms[, 2], upper = Parms[, 3], niter = 30000)
summary(GlP2)
##Estimate
ParmsGlanvilleAla <- abc_optim(fn = GlObjectiveAla, 
                               par = as.numeric(summary(GlP2)[c("mean"), ]), 
                               lb = as.numeric(summary(GlP2)[c("min"), ]), 
                               ub = as.numeric(summary(GlP2)[c("max"), ]))
ParmsGlanvilleAla$par
##Uncertainty
GlPsd2 <- modMCMC(GlObjectiveAla, p = ParmsGlanvilleAla$par, lower = Parms[1:8, 2], 
                 upper = Parms[1:8, 3], niter = 5000)
summary(GlPsd2)
#==================================================#
#Goodness of fit and simulations
source("IndividualStudies/GlFitAla.R")
SimGlanvilleAla <- GlFitAla(ParmsGlanvilleAla$par)
SimGlanvilleAla$errors
SimGlanvilleAla$R2all
#=====================Visualizing the data with simulations
Gluc <- SimGlanvilleGlucose$Simulation
Gluc$Substrate <- "Glucose"
Ala <- SimGlanvilleAla$Simulation
Ala$Substrate <- "Alanin"
SimGlanville <- rbind(Gluc, Ala)
write.csv(SimGlanville, "../Manuscript/figure_data/SimGlanville.csv", row.names = F)

##kec
ggplot(GlData, aes(Time, kec)) + geom_point(cex=6, pch=21, aes(fill = Substrate)) +
  theme_min + ylab(expression(paste(italic(k[ec])~factor))) + xlab("Time (days)") + 
  theme(legend.position = c(0.8, 0.6), legend.title = element_blank()) + 
  scale_fill_manual(values = c("white", "grey")) + scale_y_continuous(limits = c(0.2, 0.43)) +
  geom_hline(yintercept = 0.38, lwd = 1, lty = 2) +
  geom_line(data = SimGlanville, aes(Time, kec, color = Substrate), lwd = 1.2) +
  scale_color_manual(values = c("grey30", "black"))
##14C Flush
ggplot(GlData, aes(Time, Flush)) + geom_point(cex=6, pch=21, aes(fill = Substrate)) +
  theme_min + ylab(expression(paste("Flush (pmo", l^{14},C~g(DW)^{-1}, ")"))) + xlab("Time (days)") + 
  theme(legend.position = c(0.8, 0.8), legend.title = element_blank()) + 
  scale_fill_manual(values = c("white", "grey")) +
  geom_line(data = SimGlanville, aes(Time, Flush, color = Substrate), lwd = 1.2) +
  scale_color_manual(values = c("grey30", "black"))
##14C CO2
ggplot(GlData, aes(Time, CO2)) + geom_point(cex=6, pch=21, aes(fill = Substrate)) +
  theme_min + ylab(expression(paste(CO[2], " (pmo", l^{14},C~g(DW)^{-1}, ")"))) + xlab("Time (days)") + 
  theme(legend.position = c(0.8, 0.4), legend.title = element_blank()) + 
  scale_fill_manual(values = c("white", "grey")) +
  geom_line(data = SimGlanville, aes(Time, CO2, color = Substrate), lwd = 1.2) +
  scale_color_manual(values = c("grey30", "black"))
#=====================Save parameters
parsAll <- rbind(parsAll, data.frame(Study = c("Glanville et al. (2016)"),
                                     Treatment = c("Glucose"),
                                     Im = ParmsGlanvilleGlucose$par[1], Km = ParmsGlanvilleGlucose$par[2], 
                                     yA = ParmsGlanvilleGlucose$par[3],
                                     Em = ParmsGlanvilleGlucose$par[4], m = ParmsGlanvilleGlucose$par[5], 
                                     g = ParmsGlanvilleGlucose$par[6], 
                                     ne = ParmsGlanvilleGlucose$par[7], nX1 = ParmsGlanvilleGlucose$par[8], 
                                     iX1 = NA, te = NA, tX1 = NA, re = NA, rX1 = NA, 
                                     pe = NA, pX1 = NA, le = NA, lX1 = NA))
parsAll <- rbind(parsAll, data.frame(Study = c("Glanville et al. (2016)"),
                                     Treatment = c("Alanine"),
                                     Im = ParmsGlanvilleAla$par[1], Km = ParmsGlanvilleAla$par[2], 
                                     yA = ParmsGlanvilleAla$par[3],
                                     Em = ParmsGlanvilleAla$par[4], m = ParmsGlanvilleAla$par[5], 
                                     g = ParmsGlanvilleAla$par[6], 
                                     ne = ParmsGlanvilleAla$par[7], nX1 = ParmsGlanvilleAla$par[8], 
                                     iX1 = NA, te = NA, tX1 = NA, re = NA, rX1 = NA, 
                                     pe = NA, pX1 = NA, le = NA, lX1 = NA))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Study of Tessier et al. (1998)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#NOT RUN - WRONG DATA OR METHOD DESCRIPTION
#Studies also not involved in the analyses due to various issues in reprted data:
##Bremer and Kuikman, 1994; Dictor et al., 1998; Zagal, 1993; Gregorich et al., 1990 and 1991; Nguyen and 
##Guckert, 2001; Luna-Guido et al., 2001 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Study of Sparling et al. (1988, 1990)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Data
Sp90 <- read.csv("../SoilMBVariabilityData/Sparling1990.csv")
#Visualizing the data
##kec
ggplot(Sp90, aes(Cmicinit, kec)) + geom_point(cex=6, pch=21) +
  theme_min + ylab(expression(paste(italic(k[ec])~factor))) + 
  xlab(expression(paste("Initial Flush (", mu, "mol (C) ", g(DW)^{-1}, ")")))
##CO2
ggplot(Sp90, aes(Cmicinit, CO2)) + geom_point(cex=6, pch=21) +
  theme_min + ylab(expression(paste(CO[2], " (", mu, "mo", l^{14}, C~g(DW)^{-1}, ")"))) + 
  xlab(expression(paste("Initial Flush (", mu, "mol (C) ", g(DW)^{-1}, ")")))
##Flush
ggplot(Sp90, aes(Cmicinit, Flush)) + geom_point(cex=6, pch=21) +
  theme_min + ylab(expression(paste("Flush (", mu, "mo", l^{14}, C~g(DW)^{-1}, ")"))) + 
  xlab(expression(paste("Initial Flush (", mu, "mol (C) ", g(DW)^{-1}, ")")))
##Consumed substrate
ggplot(Sp90, aes(Cmicinit, Sinit-S)) + geom_point(cex=6, pch=21) +
  theme_min + ylab(expression(paste("Consumed glucose (", mu, "mo", l^{14}, C~g(DW)^{-1}, ")"))) + 
  xlab(expression(paste("Initial Flush (", mu, "mol (C) ", g(DW)^{-1}, ")")))

#Model calibration
##Reading respective python scripts
source_python("../PythonScripts/IndividualStudies/Sparling1990ODESolv.py")
##Reading objective function
source("IndividualStudies/SpObjective.R")
##Optimization
##First guess by MCMC 
SpP <- modMCMC(SpObjective, p = Parms[,1], lower = Parms[, 2], upper = Parms[, 3], niter = 30000)
summary(SpP)
##Estimate
ParmsSparling90 <- abc_optim(fn = SpObjective, 
                             par = as.numeric(summary(SpP)[c("mean"), ]), 
                             lb = as.numeric(summary(SpP)[c("min"), ]), 
                             ub = as.numeric(summary(SpP)[c("max"), ]))

ParmsSparling90$par
#Uncertainty
SpPu <- modMCMC(SpObjective, p = ParmsSparling90$par, lower = Parms[1:8, 2], upper = Parms[1:8, 3], niter = 5000)
summary(SpPu)
#==================================================#
#Goodness of fit and simulations
source("IndividualStudies/SpFit.R")
SimSparling90 <- SpFit(ParmsSparling90$par)
SimSparling90$errors
SimSparling90$R2all
#=====================Visualizing the data with simulations
SparlingSim <- as.data.frame(SimSparling90$Simulation)
SparlingSimY <- melt(SparlingSim[, 1:4])
colnames(SparlingSimY) <- c("Variable", "Measured")
SparlingSimY$Predicted <- melt(SparlingSim[, 5:8])[, 2]
write.csv(SparlingSimY, "../Manuscript/figure_data/SimSparling.csv", row.names = F)

ggplot(SparlingSimY, aes(Measured, Predicted)) + geom_point(cex = 6, pch = 21, fill = "grey") +
  theme_min + facet_wrap(~Variable, scales = "free") +
  geom_abline(intercept = 0, slope = 1)

#=====================Save parameters
parsAll <- rbind(parsAll, data.frame(Study = c("Sparling et al. (1988, 1990)"),
                                     Treatment = c("Glucose"),
                                     Im = ParmsSparling90$par[1], Km = ParmsSparling90$par[2], 
                                     yA = ParmsSparling90$par[3],
                                     Em = ParmsSparling90$par[4], m = ParmsSparling90$par[5], 
                                     g = ParmsSparling90$par[6], 
                                     ne = ParmsSparling90$par[7], nX1 = ParmsSparling90$par[8], 
                                     iX1 = NA, te = NA, tX1 = NA, re = NA, rX1 = NA, 
                                     pe = NA, pX1 = NA, le = NA, lX1 = NA))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Study of Santruckova et al. (2004)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Data
HPData = read.csv("../SoilMBVariabilityData/Santruckova2004.csv")
#Visualizing the data
##kec factor
ggplot(HPData, aes(Time, kec)) + geom_point(cex=6, pch=21, fill = "grey") +
  theme_min + ylab(expression(paste(italic(k[ec])~factor))) + xlab("Time (days)") 
##Cmic 14
ggplot(HPData, aes(Time, Flush)) + geom_point(cex=6, pch=21, fill = "grey") +
  theme_min + ylab(expression(paste("Flush (", mu, "mol(C) ", g^{-1}, ")"))) + xlab("Time (days)") 
##Cumulative respiration
ggplot(HPData, aes(Time, CO2)) + geom_point(cex=6, pch=21, fill = "grey") +
  theme_min + ylab(expression(paste("Cumulative respiration (", mu, "mol(C) ", g^{-1}, ")"))) + xlab("Time (days)") 
##Substrate concentration
ggplot(HPData, aes(Time, S)) + geom_point(cex=6, pch=21, fill = "grey") +
  theme_min + ylab(expression(paste("Glucose (", mu, "mol(C) ", g^{-1}, ")"))) + xlab("Time (days)") 

#Model calibration
##Reading respective python scripts
#source_python("../PythonScripts/IndividualStudies/Sparling1990ODESolv.py") !!!!the same calculation is performed 
##Reading objective function
source("IndividualStudies/HPObjective.R")
##Optimization
##First guess by MCMC 
HPP <- modMCMC(HPObjective, p = Parms[,1], lower = Parms[, 2], upper = Parms[, 3], niter = 30000)
summary(HPP)
##Estimate
ParmsSantruckova2004 <- abc_optim(fn = HPObjective, 
                                  par = as.numeric(summary(HPP)[c("mean"), ]), 
                                  lb = as.numeric(summary(HPP)[c("min"), ]), 
                                  ub = as.numeric(summary(HPP)[c("max"), ]))
ParmsSantruckova2004$par
#Uncertainty
HPPU <- modMCMC(HPObjective, p = ParmsSantruckova2004$par, 
                lower = Parms[1:8, 2], upper = Parms[1:8, 3], niter = 5000)
summary(HPPU)
#==================================================#
#Goodness of fit and simulations
source("IndividualStudies/HPFit.R")
SimSantruckova2004 <- HPFit(ParmsSantruckova2004$par)
SimSantruckova2004$errors
SimSantruckova2004$R2all
write.csv(SimSantruckova2004$Simulation, "../Manuscript/figure_data/SimHasan.csv", row.names = F)
#=====================Visualizing the data with simulations
##Cumulative respiration
ggplot(HPData, aes(Time, CO2)) + geom_point(cex=6, pch=21, fill = "grey") +
  theme_min + ylab(expression(paste("Cumul. resp. (", mu, "mol ", atop(14, ), "C " , g^{-1}, ")"))) + 
  xlab("Time (days)") +
  geom_line(data = (SimSantruckova2004$Simulation), aes(Time, CO2), lwd = 1.2) +
  ggtitle("A)")

##Substrate concentration
ggplot(HPData, aes(Time, S)) + geom_point(cex=6, pch=21, fill = "grey") +
  theme_min + ylab(expression(paste("Glucose (", mu, "mol ", atop(14, ), "C " , g^{-1}, ")"))) + 
  xlab("Time (days)") +
  geom_line(data = (SimSantruckova2004$Simulation), aes(Time, S), lwd = 1.2) +
  ggtitle("B)")

##Cmic 14
ggplot(HPData, aes(Time, Flush)) + geom_point(cex=6, pch=21, fill = "grey") +
  theme_min + ylab(expression(paste(CHCl[3] ," flush (", mu, "mol ", atop(14, ), "C " , g^{-1}, ")"))) + 
  xlab("Time (days)") +
  geom_line(data = (SimSantruckova2004$Simulation), aes(Time, Flush), lwd = 1.2) +
  ggtitle("C)")

##kec factor
ggplot(HPData, aes(Time, kec)) + geom_point(cex=6, pch=21, fill = "grey") +
  theme_min + ylab(expression(paste(italic(k[ec])~factor))) + xlab("Time (days)") + #ylim(0.20, 0.35) +
  geom_line(data = (SimSantruckova2004$Simulation), aes(Time, kec), lwd = 1.2) +
  scale_y_continuous(limits = c(0.2, 0.5)) +
  ggtitle("D)")
#=====================Save parameters
parsAll <- rbind(parsAll, data.frame(Study = c("Santruckova et al. (2004)"),
                                     Treatment = c("Glucose"),
                                     Im = ParmsSantruckova2004$par[1], Km = ParmsSantruckova2004$par[2], 
                                     yA = ParmsSantruckova2004$par[3],
                                     Em = ParmsSantruckova2004$par[4], m = ParmsSantruckova2004$par[5], 
                                     g = ParmsSantruckova2004$par[6], 
                                     ne = ParmsSantruckova2004$par[7], nX1 = ParmsSantruckova2004$par[8], 
                                     iX1 = NA, te = NA, tX1 = NA, re = NA, rX1 = NA, 
                                     pe = NA, pX1 = NA, le = NA, lX1 = NA))
#=====================
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Study of Bremer and van Kessel (1990)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Data
BKData = read.csv("../SoilMBVariabilityData/BremerKessel1990.csv")
#Visualizing the data
##Substrate concentration
ggplot(BKData, aes(Time, S)) + geom_point(cex=6, pch=21, fill = "grey") +
  theme_min + ylab(expression(paste("Glucose (", mu, "mol(C) ", g^{-1}, ")"))) + xlab("Time (days)") +
  facet_wrap(~Treatment, scales="free")
##Cumulative respiration
ggplot(BKData, aes(Time, CO214cumul)) + geom_point(cex=6, pch=21, fill = "grey") +
  theme_min + ylab(expression(paste("Cumulative respiration (", mu, "mol(C) ", g^{-1}, ")"))) + xlab("Time (days)") +
  facet_wrap(~Treatment, scales="free")
##Biomass 14 !This is really biomass (without applying kec factor)
ggplot(BKData, aes(Time, Cmic14)) + geom_point(cex=6, pch=21, fill = "grey") +
  theme_min + ylab(expression(paste(MB^{14}, "C(", mu, "mol(C) ", g^{-1}, ")"))) + xlab("Time (days)") +
  facet_wrap(~Treatment, scales="free")
#kec factor
ggplot(BKData, aes(Time, kec)) + geom_point(cex=6, pch=21, fill = "grey") +
  theme_min + ylab(expression(paste(italic(k[ec])~factor))) + xlab("Time (days)") +
  facet_wrap(~Treatment, scales="free")

#Model calibration
##Reading respective python scripts
source_python("../PythonScripts/IndividualStudies/Bremer1990ODESolv.py") 
source_python("../PythonScripts/IndividualStudies/DEBmodelIso.py") 
##Reading objective function
source("IndividualStudies/BKObjective.R")
##Optimization
##First guess by MCMC 
BKP <- modMCMC(BKObjective, p = Parms[,1], lower = Parms[, 2], upper = Parms[, 3], niter = 30000)
summary(BKP)
##Estimate
ParmsBremer1990 <- abc_optim(fn = BKObjective, 
                             par = as.numeric(summary(BKP)[c("mean"), ]), 
                             lb = as.numeric(summary(BKP)[c("min"), ]), 
                             ub = as.numeric(summary(BKP)[c("max"), ])) 
ParmsBremer1990$par
##Uncertainty
BKPU <- modMCMC(BKObjective, p = ParmsBremer1990$par, 
                lower = Parms[1:8, 2], upper = Parms[1:8, 3], niter = 5000)
summary(BKPU)
#==================================================#
#Goodness of fit and simulations
source("IndividualStudies/BKFit.R")
SimBremer1990 <- BKFit(ParmsBremer1990$par)
SimBremer1990$errors
SimBremer1990$R2all
write.csv(SimBremer1990$Simulation, "../Manuscript/figure_data/SimBremer.csv", row.names = F)
#=====================Visualizing the data with simulations

##Substrate concentration
ggplot(BKData, aes(Time, S)) + geom_point(cex=6, pch=21, fill = "grey") +
  theme_min + ylab(expression(paste("Glucose (", mu, "mol(C) ", g^{-1}, ")"))) + xlab("Time (days)") +
  geom_line(data = (SimBremer1990$Simulation), aes(Time, S), lwd = 1.5) +
  facet_wrap(~Treatment, scales="free")
##Cumulative respiration
ggplot(BKData, aes(Time, CO214cumul)) + geom_point(cex=6, pch=21, fill = "grey") +
  theme_min + ylab(expression(paste("Cumulative respiration (", mu, "mol(C) ", g^{-1}, ")"))) + xlab("Time (days)") +
  geom_line(data = (SimBremer1990$Simulation), aes(Time, CO2), lwd = 1.5) +
  facet_wrap(~Treatment, scales="free")
##Biomass 14
ggplot(BKData, aes(Time, Cmic14)) + geom_point(cex=6, pch=21, fill = "grey") +
  theme_min + ylab(expression(paste(MB^{14}, "C(", mu, "mol(C) ", g^{-1}, ")"))) + xlab("Time (days)") +
  geom_line(data = (SimBremer1990$Simulation), aes(Time, Cmic14), lwd = 1.5) +
  facet_wrap(~Treatment, scales="free")
#kec factor
ggplot(BKData, aes(Time, kec)) + geom_point(cex=6, pch=21, fill = "grey") +
  theme_min + ylab(expression(paste(italic(k[ec])~factor))) + xlab("Time (days)") +
  geom_line(data = (SimBremer1990$Simulation), aes(Time, kec), lwd = 1.5) +
  facet_grid(.~Treatment)
#=====================Save parameters
parsAll <- rbind(parsAll, data.frame(Study = c("Bremer and van Kessel (1990)"),
                                     Treatment = c("Glucose"),
                                     Im = ParmsBremer1990$par[1], Km = ParmsBremer1990$par[2], 
                                     yA = ParmsBremer1990$par[3],
                                     Em = ParmsBremer1990$par[4], m = ParmsBremer1990$par[5], 
                                     g = ParmsBremer1990$par[6], 
                                     ne = ParmsBremer1990$par[7], nX1 = ParmsBremer1990$par[8], 
                                     iX1 = NA, te = NA, tX1 = NA, re = NA, rX1 = NA, 
                                     pe = NA, pX1 = NA, le = NA, lX1 = NA))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Marstorp and Witter 1999~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Data
MData = read.csv("../SoilMBVariabilityData/Marstorp1999.csv")
#Visualizing the data
##Substrate concentration
ggplot(MData, aes(Time, S)) + geom_point(cex=6, pch=21, fill = "grey") +
  theme_min + ylab(expression(paste("Glucose (", mu, "mol(C) ", g^{-1}, ")"))) + xlab("Time (days)") 
##Cmic 12
ggplot(MData, aes(Time, Flush)) + geom_point(cex=6, pch=21, fill = "grey") +
  theme_min + ylab(expression(paste("Flush (", mu, "mol(C) ", g^{-1}, ")"))) + xlab("Time (days)") +
  ylim(0, 15)
##Cmic 14
ggplot(MData, aes(Time, Cmic14)) + geom_point(cex=6, pch=21, fill = "grey") +
  theme_min + ylab(expression(paste("Flush (", mu, "mol(C) ", g^{-1}, ")"))) + xlab("Time (days)") 
##Cumulative respiration
ggplot(MData, aes(Time, CO212cumul)) + geom_point(cex=6, pch=21, fill = "grey") +
  theme_min + ylab(expression(paste("Cumulative respiration (", mu, "mol(C) ", g^{-1}, ")"))) + xlab("Time (days)") 
##DNA
ggplot(MData, aes(Time, DNA)) + geom_point(cex=6, pch=21, fill = "grey") +
  theme_min + ylab(expression(paste("DNA (", mu, "mol(C) ", g^{-1}, ")"))) + xlab("Time (days)") 

#Model calibration
##Reading respective python scripts
source_python("../PythonScripts/IndividualStudies/DEBmodel.py") 
source_python("../PythonScripts/IndividualStudies/Marstorp1999ODESolv.py") 
##Reading objective function
source("IndividualStudies/MObjective.R")
##Optimization
##First guess by MCMC 
MP <- modMCMC(MObjective, p = c(Parms[,1], 0.1), 
               lower = c(Parms[, 2], 0),
               upper = c(Parms[, 3], 1), niter = 30000)
summary(MP)
##Estimate
ParmsMarstorp1999 <- abc_optim(fn = MObjective, 
                               par = as.numeric(summary(MP)[c("mean"), ]), 
                               lb = as.numeric(summary(MP)[c("min"), ]), 
                               ub = as.numeric(summary(MP)[c("max"), ])) 
ParmsMarstorp1999$par
##Uncertainty
MPU <- modMCMC(MObjective, p = ParmsMarstorp1999$par, 
              lower = c(Parms[1:8, 2], 0),
              upper = c(Parms[1:8, 3], 1), niter = 5000)
summary(MPU)
#==================================================#
#Goodness of fit and simulations
source("IndividualStudies/MFit.R")
SimMarstorp1999 <- MFit(ParmsMarstorp1999$par)
SimMarstorp1999$errors
SimMarstorp1999$R2all
write.csv(SimMarstorp1999$Simulation, "../Manuscript/figure_data/SimMarstorp.csv", row.names = F)
#=====================Save parameters
parsAll <- rbind(parsAll, data.frame(Study = c("Marstorp and Witter (1999)"),
                                     Treatment = c("Glucose"),
                                     Im = ParmsMarstorp1999$par[1], Km = ParmsMarstorp1999$par[2], 
                                     yA = ParmsMarstorp1999$par[3],
                                     Em = ParmsMarstorp1999$par[4], m = ParmsMarstorp1999$par[5], 
                                     g = ParmsMarstorp1999$par[6], 
                                     ne = ParmsMarstorp1999$par[7], nX1 = ParmsMarstorp1999$par[8], 
                                     iX1 = ParmsMarstorp1999$par[9], te = NA, tX1 = NA, re = NA, rX1 = NA, 
                                     pe = NA, pX1 = NA, le = NA, lX1 = NA))
#=====================

#=================================
# Models comparison - F statistic
#=================================
##Monod model
source_python("../PythonScripts/IndividualStudies/Monod.py")
source_python("../PythonScripts/IndividualStudies/Marstorp1999ODESolvM.py")
source("./IndividualStudies/MObjectiveM.R")
##Optimization
##First guess by MCMC 
MPM <- modMCMC(MObjectiveM, p = c(1, 25, 0.5, 1e-3, 0.3, 0.01), 
              lower = c(1e-3, 0.1, 0.05, 1e-8, 0, 0),
              upper = c(10, 3000, 1, 1, 1, 1), niter = 30000)
##Estimate
ParmsMarstorp1999M <- abc_optim(fn = MObjectiveM, par = as.numeric(summary(MPM)["mean", ]), 
                                lb = as.numeric(summary(MPM)["min", ]), 
                                ub = as.numeric(summary(MPM)["max", ])) 
ParmsMarstorp1999M$par
#Goodness of fit and simulations
source("IndividualStudies/MFitM.R")
SimMarstorp1999M <- MFitM(ParmsMarstorp1999M$par)
SimMarstorp1999M$errors
write.csv(SimMarstorp1999M$Simulation, "../Manuscript/figure_data/SimMarstorpM.csv", row.names = F)
##Pirt model
source_python("../PythonScripts/IndividualStudies/Pirt.py")
source_python("../PythonScripts/IndividualStudies/Marstorp1999ODESolvP.py")
source("./IndividualStudies/MObjectiveP.R")
##Optimization
##First guess by MCMC 
MPP <- modMCMC(MObjectiveP, p = c(1, 25, 0.5, 1e-3, 1e-3, 0.3, 0.01), 
               lower = c(1e-3, 0.1, 0.05, 1e-8, 1e-8, 0, 0),
               upper = c(10, 3000, 1, 1, 1, 1, 1), niter = 30000)
##Estimate
ParmsMarstorp1999P <- abc_optim(fn = MObjectiveP, par = as.numeric(summary(MPP)["mean", ]), 
                                lb = as.numeric(summary(MPP)["min", ]), 
                                ub = as.numeric(summary(MPP)["max", ])) 
ParmsMarstorp1999P$par
#Goodness of fit and simulations
source("IndividualStudies/MFitP.R")
SimMarstorp1999P <- MFitP(ParmsMarstorp1999P$par)
SimMarstorp1999P$errors
write.csv(SimMarstorp1999P$Simulation, "../Manuscript/figure_data/SimMarstorpP.csv", row.names = F)
#No. of observations
nt = SimMarstorp1999$errors[6]
#DEB model
DEBSSRes = SimMarstorp1999$errors[5] 
DEBpar = SimMarstorp1999$errors[7]
#Monod model
MSSRes = SimMarstorp1999M$errors[5] 
Mpar = SimMarstorp1999M$errors[7]
#Pirt model
PSSRes = SimMarstorp1999P$errors[5] 
Ppar = SimMarstorp1999P$errors[7]

#DEB vs Monod
####F value =  (M1ss - M2ss)*(nt - M2p)/M2ss/(M2p - M1p)
(MSSRes - DEBSSRes)*(nt - DEBpar)/DEBSSRes/(DEBpar - Mpar)

####associated p value
pf(q=(MSSRes - DEBSSRes)*(nt - DEBpar)/DEBSSRes/(DEBpar - Mpar), 
   df1=(DEBpar - Mpar), 
   df2=(nt - DEBpar), 
   lower.tail=F)

#DEB vs Pirt
####F value =  (M1ss - M2ss)*(nt - M2p)/M2ss/(M2p - M1p)
(PSSRes - DEBSSRes)*(nt - DEBpar)/DEBSSRes/(DEBpar - Ppar)

####associated p value
pf(q=(PSSRes - DEBSSRes)*(nt - DEBpar)/DEBSSRes/(DEBpar - Ppar), 
   df1=(DEBpar - Ppar), 
   df2=(nt - DEBpar), 
   lower.tail=F)

#=====================Visualizing the data with simulations

##Substrate concentration
ggplot(MData, aes(Time, S)) + geom_point(cex=6, pch=21, fill = "grey") +
  theme_min + ylab(expression(paste("Glucose (", mu, "mol(C) ", g^{-1}, ")"))) + xlab("Time (days)") +
  geom_line(data = (SimMarstorp1999$Simulation), aes(Time, S), lwd = 1.5) +
  geom_line(data = (SimMarstorp1999M$Simulation), aes(Time, S), lwd = 1.5, color = "grey") +
  geom_line(data = (SimMarstorp1999P$Simulation), aes(Time, S), lwd = 1.5, color = "grey30", lty = 2)
##Cumulative respiration
ggplot(MData, aes(Time, CO212cumul)) + geom_point(cex=6, pch=21, fill = "grey") +
  theme_min + ylab(expression(paste("Cumulative respiration (", mu, "mol(C) ", g^{-1}, ")"))) + xlab("Time (days)") +
  geom_line(data = (SimMarstorp1999$Simulation), aes(Time, CO2), lwd = 1.5) +
  geom_line(data = (SimMarstorp1999M$Simulation), aes(Time, CO2), lwd = 1.5, color = "grey") +
  geom_line(data = (SimMarstorp1999P$Simulation), aes(Time, CO2), lwd = 1.5, color = "grey30", lty = 2)
##Cmic
ggplot(MData, aes(Time, Flush)) + geom_point(cex=6, pch=21, fill = "grey") +
  theme_min + ylab(expression(paste("Flush (", mu, "mol(C) ", g^{-1}, ")"))) + xlab("Time (days)") +
  ylim(0, 15) + geom_line(data = (SimMarstorp1999$Simulation), aes(Time, Flush), lwd = 1.5) +
  geom_line(data = (SimMarstorp1999M$Simulation), aes(Time, Flush), lwd = 1.5, color = "grey") +
  geom_line(data = (SimMarstorp1999P$Simulation), aes(Time, Flush), lwd = 1.5, color = "grey30", lty = 2)
##DNA
ggplot(MData, aes(Time, DNA)) + geom_point(cex=6, pch=21, fill = "grey") +
  theme_min + ylab(expression(paste("DNA (", mu, "mol(C) ", g^{-1}, ")"))) + xlab("Time (days)") +
  geom_line(data = (SimMarstorp1999$Simulation), aes(Time, DNA), lwd = 1.5) +
  geom_line(data = (SimMarstorp1999M$Simulation), aes(Time, DNA), lwd = 1.5, color = "grey") +
  geom_line(data = (SimMarstorp1999P$Simulation), aes(Time, DNA), lwd = 1.5, color = "grey30", lty = 2)
##DNA/Flush
ggplot(MData, aes(Time, DNA*1000/Flush)) + geom_point(cex=6, pch=21, fill = "grey") +
  theme_min + ylab(expression(paste(frac(DNA, CHCl[3]~Flush), "   (", frac(mmol, mol), ")"))) + xlab("Time (days)") +
  geom_line(data = (SimMarstorp1999$Simulation), aes(Time, DNA*1000/Flush), lwd = 1.2) +
  ggtitle("E)")
  #geom_line(data = (SimMarstorp1999M$Simulation), aes(Time, DNA*1000/Flush), lwd = 1.5, color = "grey") +
  #geom_line(data = (SimMarstorp1999P$Simulation), aes(Time, DNA*1000/Flush), lwd = 1.5, color = "grey30", lty = 2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Ziegler et al. (2005)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Data
ZData = read.csv("../SoilMBVariabilityData/Ziegler2005.csv")
#Visualizing the data
##MBC
ggplot(ZData, aes(Time, S)) + geom_point(cex=6, pch=21, fill = "grey") +
  theme_min + ylab(expression(paste("Glucose (", mu, "mol(C) ", g^{-1}, ")"))) + xlab("Time (days)") 
##Cumulative respiration
ggplot(ZData, aes(Time, CO2cumul)) + geom_point(cex=6, pch=21, fill = "grey") +
  theme_min + ylab(expression(paste("Cumulative respiration (", mu, "mol(C) ", g^{-1}, ")"))) + xlab("Time (days)") 
##PLFA
ggplot(ZData, aes(Time, PLFA)) + geom_point(cex=6, pch=21, fill = "grey") +
  theme_min + ylab(expression(paste("PLFA (", mu, "mol(C) ", g^{-1}, ")"))) + xlab("Time (days)")

#Model calibration
##Reading respective python scripts
source_python("../PythonScripts/IndividualStudies/Ziegler2005ODESolv.py") 
##Reading objective function
source("IndividualStudies/ZObjective.R")
##Optimization
##First guess by MCMC 
ZP <- modMCMC(ZObjective, p = c(Parms[1:6, 1], 0.1), 
              lower = c(Parms[1:6, 2], 0),
              upper = c(Parms[1:6, 3], 1), niter = 30000)
summary(ZP)
##Estimate
ParmsZiegler2014 <- abc_optim(fn = ZObjective, 
                               par = as.numeric(summary(ZP)[c("mean"), ]), 
                               lb = as.numeric(summary(ZP)[c("min"), ]), 
                               ub = as.numeric(summary(ZP)[c("max"), ])) 
ParmsZiegler2014$par
##Uncertainty
ZPU <- modMCMC(ZObjective, p = ParmsZiegler2014$par, 
              lower = c(Parms[1:6, 2], 0),
              upper = c(Parms[1:6, 3], 1), niter = 5000)
summary(ZPU)
#==================================================#
#Goodness of fit and simulations
source("IndividualStudies/ZFit.R")
SimZiegler2014 <- ZFit(ParmsZiegler2014$par)
SimZiegler2014$errors
SimZiegler2014$R2all
write.csv(SimZiegler2014$Simulation, "../Manuscript/figure_data/SimZiegler.csv", row.names = F)
#=====================Save parameters
parsAll <- rbind(parsAll, data.frame(Study = c("Ziegler et al. (2005)"),
                                     Treatment = c("Glucose"),
                                     Im = ParmsZiegler2014$par[1], Km = ParmsZiegler2014$par[2], 
                                     yA = ParmsZiegler2014$par[3],
                                     Em = ParmsZiegler2014$par[4], m = ParmsZiegler2014$par[5], 
                                     g = ParmsZiegler2014$par[6], 
                                     ne = NA, nX1 = NA, 
                                     iX1 = NA, te = NA, tX1 = NA, re = NA, rX1 = NA, 
                                     pe = NA, pX1 = NA, 
                                     le = ParmsZiegler2014$par[7], lX1 = ParmsZiegler2014$par[8]))
#=====================

#=================================
# Models comparison - F statistic
#=================================
##Monod model
source_python("../PythonScripts/IndividualStudies/MonodIso.py")
source_python("../PythonScripts/IndividualStudies/Ziegler2005ODESolvM.py")
source("./IndividualStudies/ZObjectiveM.R")
##Optimization
##First guess by MCMC 
ZPM <- modMCMC(ZObjectiveM, p = c(1, 25, 0.5, 1e-3, 0.01), 
               lower = c(1e-3, 0.1, 0.05, 1e-8, 0),
               upper = c(100, 3000, 1, 1, 1), niter = 30000)
##Estimate
ParmsZiegler2005M <- abc_optim(fn = ZObjectiveM, par = as.numeric(summary(ZPM)["mean", ]), 
                                lb = as.numeric(summary(ZPM)["min", ]), 
                                ub = as.numeric(summary(ZPM)["max", ])) 
ParmsZiegler2005M$par
#Goodness of fit and simulations
source("IndividualStudies/ZFitM.R")
SimZiegler2005M <- ZFitM(ParmsZiegler2005M$par)
SimZiegler2005M$errors
write.csv(SimZiegler2005M$Simulation, "../Manuscript/figure_data/SimZieglerM.csv", row.names = F)
##Pirt model
source_python("../PythonScripts/IndividualStudies/PirtIso.py")
source_python("../PythonScripts/IndividualStudies/Ziegler2005ODESolvP.py")
source("./IndividualStudies/ZObjectiveP.R")
##Optimization
##First guess by MCMC 
ZPP <- modMCMC(ZObjectiveP, p = c(1, 25, 0.5, 1e-3, 1e-3, 0.01), 
               lower = c(1e-3, 0.1, 0.05, 1e-8, 1e-8, 0),
               upper = c(100, 3000, 1, 1, 1, 1), niter = 30000)
##Estimate
ParmsZiegler2005P <- abc_optim(fn = ZObjectiveP, par = as.numeric(summary(ZPP)["mean", ]), 
                                lb = as.numeric(summary(ZPP)["min", ]), 
                                ub = as.numeric(summary(ZPP)["max", ])) 
ParmsZiegler2005P$par
#Goodness of fit and simulations
source("IndividualStudies/ZFitP.R")
SimZiegler2005P <- ZFitP(ParmsZiegler2005P$par)
SimZiegler2005P$errors
write.csv(SimZiegler2005P$Simulation, "../Manuscript/figure_data/SimZieglerP.csv", row.names = F)
#No. of observations
nt = SimZiegler2005M$errors[6]
#DEB model
DEBSSRes = SimZiegler2014$errors[5]
DEBpar = SimZiegler2014$errors[7]
#Monod model
MSSRes = SimZiegler2005M$errors[5]
Mpar = SimZiegler2005M$errors[7]
#Pirt model
PSSRes = SimZiegler2005P$errors[5]#Weighted residual sum of squares from python (all variables standardized to 0 mean and unit variance)
Ppar = SimZiegler2005P$errors[7]

#DEB vs Monod
####F value =  (M1ss - M2ss)*(nt - M2p)/M2ss/(M2p - M1p)
(MSSRes - DEBSSRes)*(nt - DEBpar)/DEBSSRes/(DEBpar - Mpar)

####associated p value
pf(q=(MSSRes - DEBSSRes)*(nt - DEBpar)/DEBSSRes/(DEBpar - Mpar), 
   df1=(DEBpar - Mpar), 
   df2=(nt - DEBpar), 
   lower.tail=F)

#DEB vs Pirt
####F value =  (M1ss - M2ss)*(nt - M2p)/M2ss/(M2p - M1p)
(PSSRes - DEBSSRes)*(nt - DEBpar)/DEBSSRes/(DEBpar - Ppar)

####associated p value
pf(q=(PSSRes - DEBSSRes)*(nt - DEBpar)/DEBSSRes/(DEBpar - Ppar), 
   df1=(DEBpar - Ppar), 
   df2=(nt - DEBpar), 
   lower.tail=F)

##Visualizing the fit
##Glucose
ggplot(ZData, aes(Time, S)) + geom_point(cex=6, pch=21, fill = "grey") +
  theme_min + ylab(expression(paste("Glucose (", mu, "mol(C) ", g^{-1}, ")"))) + xlab("Time (days)") +
  geom_line(data = (SimZiegler2014$Simulation), aes(Time, S), lwd = 1.5) +
  geom_line(data = (SimZiegler2005M$Simulation), aes(Time, S), lwd = 1.5, color = "grey") +
  geom_line(data = (SimZiegler2005P$Simulation), aes(Time, S), lwd = 1.5, color = "grey30", lty = 2) 
##Cumulative respiration
ggplot(ZData, aes(Time, CO2cumul)) + geom_point(cex=6, pch=21, fill = "grey") +
  theme_min + ylab(expression(paste("Cumulative respiration (", mu, "mol(C) ", g^{-1}, ")"))) + xlab("Time (days)") +
  geom_line(data = (SimZiegler2014$Simulation), aes(Time, CO2), lwd = 1.5) +
  geom_line(data = (SimZiegler2005M$Simulation), aes(Time, CO2), lwd = 1.5, color = "grey") +
  geom_line(data = (SimZiegler2005P$Simulation), aes(Time, CO2), lwd = 1.5, color = "grey30", lty = 2) 
##PLFA
ggplot(ZData, aes(Time, PLFA)) + geom_point(cex=6, pch=21, fill = "grey") +
  theme_min + ylab(expression(paste("PLFA (", mu, "mol(C) ", g^{-1}, ")"))) + xlab("Time (days)") +
  geom_line(data = (SimZiegler2014$Simulation), aes(Time, PLFA), lwd = 1.5) +
  geom_line(data = (SimZiegler2005M$Simulation), aes(Time, PLFA), lwd = 1.5, color = "grey") +
  geom_line(data = (SimZiegler2005P$Simulation), aes(Time, PLFA), lwd = 1.5, color = "grey30", lty = 2) 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Tsai et al. (1997)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Data
TData = read.csv("../SoilMBVariabilityData/Tsai1997.csv")
#Visualizing the data
##MBC
ggplot(TData, aes(Time, Cmic)) + geom_point(cex=6, pch=21, fill = "grey") +
  theme_min + ylab(expression(paste("Flush (", mu, "mol(C) ", g^{-1}, ")"))) + xlab("Time (days)") +
  facet_grid(.~Treatment)
##Cumulative respiration
ggplot(TData, aes(Time, CO2cumul)) + geom_point(cex=6, pch=21, fill = "grey") +
  theme_min + ylab(expression(paste("Cumulative respiration (", mu, "mol(C) ", g^{-1}, ")"))) + xlab("Time (days)") +
  facet_grid(.~Treatment)
##ATP
ggplot(TData, aes(Time, ATP*1000)) + geom_point(cex=6, pch=21, fill = "grey") +
  theme_min + ylab(expression(paste("ATP (nmol(ATP) ", g^{-1}, ")"))) + xlab("Time (days)") +
  facet_grid(.~Treatment)
##ATP to Cmic
ggplot(TData, aes(Time, ATP/Cmic)) + geom_point(cex=6, pch=21, fill = "grey") +
  theme_min + ylab(expression(paste("ATP/Flush (mmol(ATP)/mol(C))"))) + xlab("Time (days)") +
  facet_grid(.~Treatment)

#Model calibration
##Reading respective python scripts
source_python("../PythonScripts/IndividualStudies/Tsai1997ODESolv.py") 
source_python("../PythonScripts/IndividualStudies/DEBmodel.py") 
##Reading objective function
source("IndividualStudies/TObjective.R")
##Optimization
##First guess by MCMC 
TP <- modMCMC(TObjective, p = c(Parms[,1], 1e-4, 1e-4), 
              lower = c(Parms[, 2], 0, 0),
              upper = c(Parms[, 3], 0.1, 0.1), niter = 30000)
summary(TP)
##Estimate
ParmsTsai1997 <- abc_optim(fn = TObjective, 
                              par = as.numeric(summary(TP)[c("mean"), ]), 
                              lb = as.numeric(summary(TP)[c("min"), ]), 
                              ub = as.numeric(summary(TP)[c("max"), ])) 
ParmsTsai1997$par
##Uncertainty
TPU <- modMCMC(TObjective, p = ParmsTsai1997$par, 
              lower = c(Parms[1:8, 2], 0, 0),
              upper = c(Parms[1:8, 3], 0.1, 0.1), niter = 5000)
summary(TPU)
#==================================================#
#Goodness of fit and simulations
source("IndividualStudies/TFit.R")
SimTsai1997 <- TFit(ParmsTsai1997$par)
SimTsai1997$errors
SimTsai1997$R2all
write.csv(SimTsai1997$Simulation, "../Manuscript/figure_data/SimTsai.csv", row.names = F)
#=====================Save parameters
parsAll <- rbind(parsAll, data.frame(Study = c("Tsai et al. (1997)"),
                                     Treatment = c("Glucose"),
                                     Im = ParmsTsai1997$par[1], Km = ParmsTsai1997$par[2], 
                                     yA = ParmsTsai1997$par[3],
                                     Em = ParmsTsai1997$par[4], m = ParmsTsai1997$par[5], 
                                     g = ParmsTsai1997$par[6], 
                                     ne = ParmsTsai1997$par[7], nX1 = ParmsTsai1997$par[8], 
                                     iX1 = NA, te = ParmsTsai1997$par[9], 
                                     tX1 = ParmsTsai1997$par[10], re = NA, rX1 = NA, 
                                     pe = NA, pX1 = NA, 
                                     le = NA, lX1 = NA))
#=====================

#=================================
# Models comparison - F statistic
#=================================
##Monod model
source_python("../PythonScripts/IndividualStudies/Tsai1997ODESolvM.py")
source_python("../PythonScripts/IndividualStudies/Monod.py")
source("./IndividualStudies/TObjectiveM.R")
##Optimization
##First guess by MCMC 
TPM <- modMCMC(TObjectiveM, p = c(1, 25, 0.5, 1e-3, 0.4, 0.01), 
               lower = c(1e-3, 0.1, 0.05, 1e-8, 0, 0),
               upper = c(100, 3000, 1, 1, 1, 1), niter = 30000)
##Estimate
ParmsTsai1997M <- abc_optim(fn = TObjectiveM, par = as.numeric(summary(TPM)["mean", ]), 
                               lb = as.numeric(summary(TPM)["min", ]), 
                               ub = as.numeric(summary(TPM)["max", ])) 
ParmsTsai1997M$par
#Goodness of fit and simulations
source("IndividualStudies/TFitM.R")
SimTsai1997M <- TFitM(ParmsTsai1997M$par)
SimTsai1997M$errors
write.csv(SimTsai1997M$Simulation, "../Manuscript/figure_data/SimTsaiM.csv", row.names = F)
##Pirt model
source_python("../PythonScripts/IndividualStudies/Tsai1997ODESolvP.py")
source_python("../PythonScripts/IndividualStudies/Pirt.py")
source("./IndividualStudies/TObjectiveP.R")
##Optimization
##First guess by MCMC 
TPP <- modMCMC(TObjectiveP, p = c(1, 25, 0.5, 1e-3, 1e-3, 0.3, 0.01), 
               lower = c(1e-3, 0.1, 0.05, 1e-8, 1e-8, 0, 0),
               upper = c(100, 3000, 1, 1, 1, 1, 1), niter = 30000)
##Estimate
ParmsTsai1997P <- abc_optim(fn = TObjectiveP, par = as.numeric(summary(TPP)["mean", ]), 
                               lb = as.numeric(summary(TPP)["min", ]), 
                               ub = as.numeric(summary(TPP)["max", ])) 
ParmsTsai1997P$par
#Goodness of fit and simulations
source("IndividualStudies/TFitP.R")
SimTsai1997P <- TFitP(ParmsTsai1997P$par)
SimTsai1997P$errors
write.csv(SimTsai1997P$Simulation, "../Manuscript/figure_data/SimTsaiP.csv", row.names = F)
#No. of observations
nt = SimTsai1997P$errors[6]
#DEB model
DEBSSRes = SimTsai1997$errors[5]
DEBpar = SimTsai1997$errors[7]
#Monod model
MSSRes = SimTsai1997M$errors[5]
Mpar = SimTsai1997M$errors[7]
#Pirt model
PSSRes = SimTsai1997P$errors[5]#Weighted residual sum of squares from python (all variables standardized to 0 mean and unit variance)
Ppar = SimTsai1997P$errors[7]

#DEB vs Monod
####F value =  (M1ss - M2ss)*(nt - M2p)/M2ss/(M2p - M1p)
(MSSRes - DEBSSRes)*(nt - DEBpar)/DEBSSRes/(DEBpar - Mpar)

####associated p value
pf(q=(MSSRes - DEBSSRes)*(nt - DEBpar)/DEBSSRes/(DEBpar - Mpar), 
   df1=(DEBpar - Mpar), 
   df2=(nt - DEBpar), 
   lower.tail=F)

#DEB vs Pirt
####F value =  (M1ss - M2ss)*(nt - M2p)/M2ss/(M2p - M1p)
(PSSRes - DEBSSRes)*(nt - DEBpar)/DEBSSRes/(DEBpar - Ppar)

####associated p value
pf(q=(PSSRes - DEBSSRes)*(nt - DEBpar)/DEBSSRes/(DEBpar - Ppar), 
   df1=(DEBpar - Ppar), 
   df2=(nt - DEBpar), 
   lower.tail=F)

##Visualizing the fit
##MBC
ggplot(TData, aes(Time, Cmic)) + geom_point(cex=6, pch=21, fill = "grey") +
  theme_min + ylab(expression(paste("Flush (", mu, "mol(C) ", g^{-1}, ")"))) + xlab("Time (days)") +
  geom_line(data = (SimTsai1997$Simulation), aes(Time, Flush), lwd = 1.5) +
  geom_line(data = (SimTsai1997M$Simulation), aes(Time, Flush), lwd = 1.5, color = "grey") +
  geom_line(data = (SimTsai1997P$Simulation), aes(Time, Flush), lwd = 1.5, color = "grey30", lty = 2) +
  facet_grid(.~Treatment)
##Cumulative respiration
ggplot(TData, aes(Time, CO2cumul)) + geom_point(cex=6, pch=21, fill = "grey") +
  theme_min + ylab(expression(paste("Cumulative respiration (", mu, "mol(C) ", g^{-1}, ")"))) + xlab("Time (days)") +
  geom_line(data = (SimTsai1997$Simulation), aes(Time, CO2), lwd = 1.5) +
  geom_line(data = (SimTsai1997M$Simulation), aes(Time, CO2), lwd = 1.5, color = "grey") +
  geom_line(data = (SimTsai1997P$Simulation), aes(Time, CO2), lwd = 1.5, color = "grey30", lty = 2) +
  facet_grid(.~Treatment)
##ATP
ggplot(TData, aes(Time, ATP*1000)) + geom_point(cex=6, pch=21, fill = "grey") +
  theme_min + ylab(expression(paste("ATP (nmol(ATP) ", g^{-1}, ")"))) + xlab("Time (days)") +
  geom_line(data = (SimTsai1997$Simulation), aes(Time, ATP*1000), lwd = 1.5) +
  geom_line(data = (SimTsai1997M$Simulation), aes(Time, ATP*1000), lwd = 1.5, color = "grey") +
  geom_line(data = (SimTsai1997P$Simulation), aes(Time, ATP*1000), lwd = 1.5, color = "grey30", lty = 2) +
  facet_grid(.~Treatment)

##ATP to Cmic
ggplot(TData, aes(Time, ATP*1000/Cmic)) + geom_point(cex=6, pch=21, fill = "grey") +
  theme_min + ylab(expression(paste(frac(ATP, CHCl[3]~flush),"    (", frac(mmol, mol),")"))) + xlab("Time (days)") +
  geom_line(data = (SimTsai1997$Simulation), aes(Time, ATP*1000/Flush), lwd = 1.5) +
  #geom_line(data = (SimTsai1997M$Simulation), aes(Time, ATP*1000/Flush), lwd = 1.5, color = "grey") +
  #geom_line(data = (SimTsai1997P$Simulation), aes(Time, ATP*1000/Flush), lwd = 1.5, color = "grey30", lty = 2) +
  facet_grid(.~Treatment)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Joergensen and Raubuch (2002)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Data
JData = read.csv("../SoilMBVariabilityData/Joergensen2002.csv")
#Visualizing the data
##Glucose
ggplot(JData, aes(Time, S)) + geom_point(cex=6, pch=21, fill = "grey") +
  theme_min + ylab(expression(paste("Glucose (", mu, "mol(C) ", g^{-1}, ")"))) + xlab("Time (days)") 
##MBC
ggplot(JData, aes(Time, Cmic)) + geom_point(cex=6, pch=21, fill = "grey") +
  theme_min + ylab(expression(paste("Flush (", mu, "mol(C) ", g^{-1}, ")"))) + xlab("Time (days)") 
##ATP
ggplot(JData, aes(Time, ATP*1000)) + geom_point(cex=6, pch=21, fill = "grey") +
  theme_min + ylab(expression(paste("ATP (nmol(ATP) ", g^{-1}, ")"))) + xlab("Time (days)") 

##ATP to Cmic
ggplot(JData, aes(Time, ATP*1000/Cmic)) + geom_point(cex=6, pch=21, fill = "grey") +
  theme_min + ylab(expression(paste("ATP/Flush (mmol(ATP)/mol(C))"))) + xlab("Time (days)") 

#Model calibration
##Reading respective python scripts
source_python("../PythonScripts/IndividualStudies/Joerg2002ODESolv.py") 
source_python("../PythonScripts/IndividualStudies/DEBmodel.py") 
##Reading objective function
source("IndividualStudies/JObjective.R")
##Optimization
##First guess by MCMC 
JP <- modMCMC(JObjective, p = c(Parms[,1], 1e-4, 1e-4), 
              lower = c(Parms[, 2], 0, 0),
              upper = c(Parms[, 3], 0.1, 0.1), niter = 30000)
summary(JP)
##Estimate
ParmsJoerg2002 <- abc_optim(fn = JObjective, 
                           par = as.numeric(summary(JP)[c("mean"), ]), 
                           lb = as.numeric(summary(JP)[c("min"), ]), 
                           ub = as.numeric(summary(JP)[c("max"), ])) 
ParmsJoerg2002$par
##Uncertainty
JPU <- modMCMC(JObjective, p = ParmsJoerg2002$par, 
              lower = c(Parms[1:8, 2], 0, 0),
              upper = c(Parms[1:8, 3], 0.1, 0.1), niter = 5000)
summary(JPU)
#==================================================#
#Goodness of fit and simulations
source("IndividualStudies/JFit.R")
SimJoerg2002 <- JFit(ParmsJoerg2002$par)
SimJoerg2002$errors
SimJoerg2002$R2all
write.csv(SimJoerg2002$Simulation, "../Manuscript/figure_data/SimJoerg.csv", row.names = F)
#=====================Save parameters
parsAll <- rbind(parsAll, data.frame(Study = c("Joergensen and Raubuch (2002)"),
                                     Treatment = c("Glucose"),
                                     Im = ParmsJoerg2002$par[1], Km = ParmsJoerg2002$par[2], 
                                     yA = ParmsJoerg2002$par[3],
                                     Em = ParmsJoerg2002$par[4], m = ParmsJoerg2002$par[5], 
                                     g = ParmsJoerg2002$par[6], 
                                     ne = ParmsJoerg2002$par[7], nX1 = ParmsJoerg2002$par[8], 
                                     iX1 = NA, te = ParmsJoerg2002$par[9], 
                                     tX1 = ParmsJoerg2002$par[10], re = NA, rX1 = NA, 
                                     pe = NA, pX1 = NA, 
                                     le = NA, lX1 = NA))
#=====================

#=================================
# Models comparison - F statistic
#=================================
##Monod model
source_python("../PythonScripts/IndividualStudies/Joerg2002ODESolvM.py")
source_python("../PythonScripts/IndividualStudies/Monod.py")
source("./IndividualStudies/JObjectiveM.R")
##Optimization
##First guess by MCMC 
JPM <- modMCMC(JObjectiveM, p = c(1, 25, 0.5, 1e-3, 0.4, 0.01), 
               lower = c(1e-3, 0.1, 0.05, 1e-8, 0, 0),
               upper = c(100, 3000, 1, 1, 1, 1), niter = 30000)
##Estimate
ParmsJoerg2002M <- abc_optim(fn = JObjectiveM, par = as.numeric(summary(JPM)["mean", ]), 
                            lb = as.numeric(summary(JPM)["min", ]), 
                            ub = as.numeric(summary(JPM)["max", ])) 
ParmsJoerg2002M$par
#Goodness of fit and simulations
source("IndividualStudies/JFitM.R")
SimJoerg2002M <- JFitM(ParmsJoerg2002M$par)
SimJoerg2002M$errors
write.csv(SimJoerg2002M$Simulation, "../Manuscript/figure_data/SimJoergM.csv", row.names = F)
##Pirt model
source_python("../PythonScripts/IndividualStudies/Joerg2002ODESolvP.py")
source_python("../PythonScripts/IndividualStudies/Pirt.py")
source("./IndividualStudies/JObjectiveP.R")
##Optimization
##First guess by MCMC 
JPP <- modMCMC(JObjectiveP, p = c(1, 25, 0.5, 1e-3, 1e-3, 0.3, 0.01), 
               lower = c(1e-3, 0.1, 0.05, 1e-8, 1e-8, 0, 0),
               upper = c(100, 3000, 1, 1, 1, 1, 1), niter = 30000)
##Estimate
ParmsJoerg2002P <- abc_optim(fn = JObjectiveP, par = as.numeric(summary(JPP)["mean", ]), 
                            lb = as.numeric(summary(JPP)["min", ]), 
                            ub = as.numeric(summary(JPP)["max", ])) 
ParmsJoerg2002P$par
#Goodness of fit and simulations
source("IndividualStudies/JFitP.R")
SimJoerg2002P <- JFitP(ParmsJoerg2002P$par)
SimJoerg2002P$errors
write.csv(SimJoerg2002P$Simulation, "../Manuscript/figure_data/SimJoergP.csv", row.names = F)
#No. of observations
nt = SimJoerg2002P$errors[6]
#DEB model
DEBSSRes = SimJoerg2002$errors[5]
DEBpar = SimJoerg2002$errors[7]
#Monod model
MSSRes = SimJoerg2002M$errors[5]
Mpar = SimJoerg2002M$errors[7]
#Pirt model
PSSRes = SimJoerg2002P$errors[5]#Weighted residual sum of squares from python (all variables standardized to 0 mean and unit variance)
Ppar = SimJoerg2002P$errors[7]

#DEB vs Monod
####F value =  (M1ss - M2ss)*(nt - M2p)/M2ss/(M2p - M1p)
(MSSRes - DEBSSRes)*(nt - DEBpar)/DEBSSRes/(DEBpar - Mpar)

####associated p value
pf(q=(MSSRes - DEBSSRes)*(nt - DEBpar)/DEBSSRes/(DEBpar - Mpar), 
   df1=(DEBpar - Mpar), 
   df2=(nt - DEBpar), 
   lower.tail=F)

#DEB vs Pirt
####F value =  (M1ss - M2ss)*(nt - M2p)/M2ss/(M2p - M1p)
(PSSRes - DEBSSRes)*(nt - DEBpar)/DEBSSRes/(DEBpar - Ppar)

####associated p value
pf(q=(PSSRes - DEBSSRes)*(nt - DEBpar)/DEBSSRes/(DEBpar - Ppar), 
   df1=(DEBpar - Ppar), 
   df2=(nt - DEBpar), 
   lower.tail=F)

##Visualizing the fit
##Glucose
ggplot(JData, aes(Time, S)) + geom_point(cex=6, pch=21, fill = "grey") +
  theme_min + ylab(expression(paste("Glucose (", mu, "mol(C) ", g^{-1}, ")"))) + xlab("Time (days)") +
  geom_line(data = (SimJoerg2002$Simulation), aes(Time, S), lwd = 1.5) +
  geom_line(data = (SimJoerg2002M$Simulation), aes(Time, S), lwd = 1.5, color = "grey") +
  geom_line(data = (SimJoerg2002P$Simulation), aes(Time, S), lwd = 1.5, color = "grey30", lty = 2)
##MBC
ggplot(JData, aes(Time, Cmic)) + geom_point(cex=6, pch=21, fill = "grey") +
  theme_min + ylab(expression(paste("Flush (", mu, "mol(C) ", g^{-1}, ")"))) + xlab("Time (days)") +
  geom_line(data = (SimJoerg2002$Simulation), aes(Time, Flush), lwd = 1.5) +
  geom_line(data = (SimJoerg2002M$Simulation), aes(Time, Flush), lwd = 1.5, color = "grey") +
  geom_line(data = (SimJoerg2002P$Simulation), aes(Time, Flush), lwd = 1.5, color = "grey30", lty = 2)
##ATP
ggplot(JData, aes(Time, ATP*1000)) + geom_point(cex=6, pch=21, fill = "grey") +
  theme_min + ylab(expression(paste("ATP (nmol(ATP) ", g^{-1}, ")"))) + xlab("Time (days)") +
  geom_line(data = (SimJoerg2002$Simulation), aes(Time, ATP*1000), lwd = 1.5) +
  geom_line(data = (SimJoerg2002M$Simulation), aes(Time, ATP*1000), lwd = 1.5, color = "grey") +
  geom_line(data = (SimJoerg2002P$Simulation), aes(Time, ATP*1000), lwd = 1.5, color = "grey30", lty = 2)

##ATP to Cmic
ggplot(JData, aes(Time, ATP*1000/Cmic)) + geom_point(cex=6, pch=21, fill = "grey") +
  theme_min + ylab(expression(paste(frac(ATP, CHCl[3]~flush), "   (", frac(mmol, mol), ")"))) + xlab("Time (days)") +
  geom_line(data = (SimJoerg2002$Simulation), aes(Time, ATP*1000/Flush), lwd = 1.2) #+
  geom_line(data = (SimJoerg2002M$Simulation), aes(Time, ATP*1000/Flush), lwd = 1.5, color = "grey") +
  geom_line(data = (SimJoerg2002P$Simulation), aes(Time, ATP*1000/Flush), lwd = 1.5, color = "grey30", lty = 2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Nannipieri et al. (1977)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Data
NData = read.csv("../SoilMBVariabilityData/Nannipieri1977.csv")
#Visualizing the data
##CO2
ggplot(NData, aes(Time, CO2)) + geom_point(cex=6, pch=21, fill = "grey") + facet_grid(.~Treatment) +
  theme_min + ylab(expression(paste("Cumulative respoiration (", mu, "mol(C) ", g^{-1}, ")"))) + xlab("Time (days)") 
##Biomass
ggplot(NData, aes(Time, W/1000)) + geom_point(cex=6, pch=21, fill = "grey") + facet_grid(.~Treatment) +
  theme_min + ylab(expression(paste("Biomass (", "mg ", g^{-1}, ")"))) + xlab("Time (days)") 
##ATP
ggplot(NData, aes(Time, ATP*1000)) + geom_point(cex=6, pch=21, fill = "grey") + facet_grid(.~Treatment) +
  theme_min + ylab(expression(paste("ATP (nmol(ATP) ", g^{-1}, ")"))) + xlab("Time (days)") 

#Model calibration
##Reading respective python scripts
source_python("../PythonScripts/IndividualStudies/Nanni1977ODESolv.py") 
source_python("../PythonScripts/IndividualStudies/DEBmodel.py") 
##Reading objective function
source("IndividualStudies/NObjective.R")
##Optimization
##First guess by MCMC 
NP <- modMCMC(NObjective, p = c(Parms[1:6,1], 100, 100, 1e-4, 1e-4), 
              lower = c(Parms[1:6, 2], 1, 1, 0, 0),
              upper = c(Parms[1:6, 3], 1000, 1000, 1, 1), niter = 100000)
summary(NP)
##Estimate
ParmsNanni1977 <- abc_optim(fn = NObjective, 
                            par = as.numeric(summary(NP)[c("mean"), ]), 
                            lb = as.numeric(summary(NP)[c("min"), ]), 
                            ub = as.numeric(summary(NP)[c("max"), ])) 
ParmsNanni1977$par
##Uncertainty
NPU <- modMCMC(NObjective, p = ParmsNanni1977$par, 
              lower = c(Parms[1:6, 2], 1, 1, 0, 0),
              upper = c(Parms[1:6, 3], 1000, 1000, 1, 1), niter = 5000)
summary(NPU)
#==================================================#
#Goodness of fit and simulations
source("IndividualStudies/NFit.R")
SimNanni1977 <- NFit(ParmsNanni1977$par)
SimNanni1977$errors
SimNanni1977$R2all
write.csv(SimNanni1977$Simulation, "../Manuscript/figure_data/SimNannipieri.csv", row.names = F)
#=====================Save parameters
parsAll <- rbind(parsAll, data.frame(Study = c("Nannipieri et al. (1977)"),
                                     Treatment = c("Glucose"),
                                     Im = ParmsNanni1977$par[1], Km = ParmsNanni1977$par[2], 
                                     yA = ParmsNanni1977$par[3],
                                     Em = ParmsNanni1977$par[4], m = ParmsNanni1977$par[5], 
                                     g = ParmsNanni1977$par[6], 
                                     ne = NA, nX1 = NA, 
                                     iX1 = NA, te = ParmsNanni1977$par[9], 
                                     tX1 = ParmsNanni1977$par[10], re = NA, rX1 = NA, 
                                     pe = NA, pX1 = NA, 
                                     le = NA, lX1 = NA))
#=====================

#=================================
# Models comparison - F statistic
#=================================
##Monod model
source_python("../PythonScripts/IndividualStudies/Nanni1977ODESolvM.py")
source_python("../PythonScripts/IndividualStudies/Monod.py")
source("./IndividualStudies/NObjectiveM.R")
##Optimization
##First guess by MCMC 
NPM <- modMCMC(NObjectiveM, p = c(1, 25, 0.5, 1e-3, 10, 0.01), 
               lower = c(1e-3, 0.1, 0.05, 1e-8, 1, 0),
               upper = c(100, 3000, 1, 1, 500, 1), niter = 30000)
##Estimate
ParmsNanni1977M <- abc_optim(fn = NObjectiveM, par = as.numeric(summary(NPM)["mean", ]), 
                             lb = as.numeric(summary(NPM)["min", ]), 
                             ub = as.numeric(summary(NPM)["max", ])) 
ParmsNanni1977M$par
#Goodness of fit and simulations
source("IndividualStudies/NFitM.R")
SimNanni1977M <- NFitM(ParmsNanni1977M$par)
SimNanni1977M$errors
write.csv(SimNanni1977M$Simulation, "../Manuscript/figure_data/SimNannipieriM.csv", row.names = F)
##Pirt model
source_python("../PythonScripts/IndividualStudies/Nanni1977ODESolvP.py")
source_python("../PythonScripts/IndividualStudies/Pirt.py")
source("./IndividualStudies/NObjectiveP.R")
##Optimization
##First guess by MCMC 
NPP <- modMCMC(NObjectiveP, p = c(1, 25, 0.5, 1e-3, 1e-3, 10, 0.01), 
               lower = c(1e-3, 0.1, 0.05, 1e-8, 1e-8, 1, 0),
               upper = c(100, 3000, 1, 1, 1, 1000, 1), niter = 30000)
##Estimate
ParmsNanni1977P <- abc_optim(fn = NObjectiveP, par = as.numeric(summary(NPP)["mean", ]), 
                             lb = as.numeric(summary(NPP)["min", ]), 
                             ub = as.numeric(summary(NPP)["max", ])) 
ParmsNanni1977P$par
#Goodness of fit and simulations
source("IndividualStudies/NFitP.R")
SimNanni1977P <- NFitP(ParmsNanni1977P$par)
SimNanni1977P$errors
write.csv(SimNanni1977P$Simulation, "../Manuscript/figure_data/SimNannipieriP.csv", row.names = F)
#No. of observations
nt = SimNanni1977P$errors[6]
#DEB model
DEBSSRes = SimNanni1977$errors[5]
DEBpar = SimNanni1977$errors[7]
#Monod model
MSSRes = SimNanni1977M$errors[5]
Mpar = SimNanni1977M$errors[7]
#Pirt model
PSSRes = SimNanni1977P$errors[5]#Weighted residual sum of squares from python (all variables standardized to 0 mean and unit variance)
Ppar = SimNanni1977P$errors[7]

#DEB vs Monod
####F value =  (M1ss - M2ss)*(nt - M2p)/M2ss/(M2p - M1p)
(MSSRes - DEBSSRes)*(nt - DEBpar)/DEBSSRes/(DEBpar - Mpar)

####associated p value
pf(q=(MSSRes - DEBSSRes)*(nt - DEBpar)/DEBSSRes/(DEBpar - Mpar), 
   df1=(DEBpar - Mpar), 
   df2=(nt - DEBpar), 
   lower.tail=F)

#DEB vs Pirt
####F value =  (M1ss - M2ss)*(nt - M2p)/M2ss/(M2p - M1p)
(PSSRes - DEBSSRes)*(nt - DEBpar)/DEBSSRes/(DEBpar - Ppar)

####associated p value
pf(q=(PSSRes - DEBSSRes)*(nt - DEBpar)/DEBSSRes/(DEBpar - Ppar), 
   df1=(DEBpar - Ppar), 
   df2=(nt - DEBpar), 
   lower.tail=F)

##Visualizing the fit
##CO2
ggplot(NData, aes(Time, CO2)) + geom_point(cex=6, pch=21, fill = "grey") + facet_grid(.~Treatment) +
  theme_min + ylab(expression(paste("Cumulative respiration (", mu, "mol(C) ", g^{-1}, ")"))) + xlab("Time (days)") +
  geom_line(data = (SimNanni1977$Simulation), aes(Time, CO2), lwd = 1.5) +
  geom_line(data = (SimNanni1977M$Simulation), aes(Time, CO2), lwd = 1.5, color = "grey") +
  geom_line(data = (SimNanni1977P$Simulation), aes(Time, CO2), lwd = 1.5, color = "grey30", lty = 2)
##Biomass
ggplot(NData, aes(Time, W/1000)) + geom_point(cex=6, pch=21, fill = "grey") + facet_grid(.~Treatment) +
  theme_min + ylab(expression(paste("Biomass (mg ", g^{-1}, ")"))) + xlab("Time (days)") +
  geom_line(data = (SimNanni1977$Simulation), aes(Time, W/1000), lwd = 1.5) +
  geom_line(data = (SimNanni1977M$Simulation), aes(Time, W/1000), lwd = 1.5, color = "grey") +
  geom_line(data = (SimNanni1977P$Simulation), aes(Time, W/1000), lwd = 1.5, color = "grey30", lty = 2)
##ATP
ggplot(NData, aes(Time, ATP)) + geom_point(cex=6, pch=21, fill = "grey") + facet_grid(.~Treatment) +
  theme_min + ylab(expression(paste("ATP (", mu, "mol(ATP) ", g^{-1}, ")"))) + xlab("Time (days)") +
  geom_line(data = (SimNanni1977$Simulation), aes(Time, ATP), lwd = 1.5) +
  geom_line(data = (SimNanni1977M$Simulation), aes(Time, ATP), lwd = 1.5, color = "grey") +
  geom_line(data = (SimNanni1977P$Simulation), aes(Time, ATP), lwd = 1.5, color = "grey30", lty = 2)
##ATP/W
ggplot(NData[NData$Treatment=="A", ], aes(Time, ATP/W*1000)) + geom_point(cex=6, pch=21, fill = "grey") + 
  theme_min + ylab(expression(paste(frac(ATP, W), "   (", frac(mmol, g(Cells)),")"))) + xlab("Time (days)") +
  geom_line(data = subset(SimNanni1977$Simulation, Treatment == "A"), aes(Time, ATP/W*1000), lwd = 1.2) #+
  geom_line(data = (SimNanni1977M$Simulation), aes(Time, ATP/W), lwd = 1.5, color = "grey") +
  geom_line(data = (SimNanni1977P$Simulation), aes(Time, ATP/W), lwd = 1.5, color = "grey30", lty = 2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Blagodatskaya et al. (2014)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Data
BData = read.csv("../SoilMBVariabilityData/Blagodatskaya2014.csv")
#Visualizing the data
##CO2
ggplot(BData, aes(Time, CO2)) + geom_point(cex=6, pch=21, fill = "grey") + facet_grid(.~Treatment) +
  theme_min + ylab(expression(paste("Cumulative respoiration (", mu, "mol(C) ", g^{-1}, ")"))) + xlab("Time (days)") 
##DNA
ggplot(BData, aes(Time, DNA)) + geom_point(cex=6, pch=21, fill = "grey") + facet_grid(.~Treatment) +
  theme_min + ylab(expression(paste("DNA (nmol(C) ", g^{-1}, ")"))) + xlab("Time (days)") 

#Model calibration
##Reading respective python scripts
source_python("../PythonScripts/IndividualStudies/Blag2014ODESolv.py") 
source_python("../PythonScripts/IndividualStudies/DEBmodel.py") 
##Reading objective function
source("IndividualStudies/BObjective.R")
##Optimization
##First guess by MCMC 
BP <- modMCMC(BObjective, p = c(Parms[1:6,1], 1e-3), 
              lower = c(Parms[1:6, 2], 0),
              upper = c(Parms[1:6, 3], 1), niter = 30000)
summary(BP)
##Estimate
ParmsBlag2014 <- abc_optim(fn = BObjective, 
                            par = as.numeric(summary(BP)[c("mean"), ]), 
                            lb = as.numeric(summary(BP)[c("min"), ]), 
                            ub = as.numeric(summary(BP)[c("max"), ])) 
ParmsBlag2014$par
##Uncertainty
BPU <- modMCMC(BObjective, p = ParmsBlag2014$par, 
              lower = c(Parms[1:6, 2], 0),
              upper = c(Parms[1:6, 3], 1), niter = 5000)
summary(BPU)
#==================================================#
#Goodness of fit and simulations
source("IndividualStudies/BFit.R")
SimBlag2014 <- BFit(ParmsBlag2014$par)
SimBlag2014$errors
SimBlag2014$R2all
write.csv(SimBlag2014$Simulation, "../Manuscript/figure_data/SimBlag.csv", row.names = F)
#=====================Save parameters
parsAll <- rbind(parsAll, data.frame(Study = c("Blagodatskaya et al. (2014)"),
                                     Treatment = c("Glucose"),
                                     Im = ParmsBlag2014$par[1], Km = ParmsBlag2014$par[2], 
                                     yA = ParmsBlag2014$par[3],
                                     Em = ParmsBlag2014$par[4], m = ParmsBlag2014$par[5], 
                                     g = ParmsBlag2014$par[6], 
                                     ne = NA, nX1 = NA, 
                                     iX1 = ParmsBlag2014$par[7], te = NA, 
                                     tX1 = NA, re = NA, rX1 = NA, 
                                     pe = NA, pX1 = NA, 
                                     le = NA, lX1 = NA))
#=====================

#=================================
# Models comparison - F statistic
#=================================
##Monod model
source_python("../PythonScripts/IndividualStudies/Blag2014ODESolvM.py")
source_python("../PythonScripts/IndividualStudies/Monod.py")
source("./IndividualStudies/BObjectiveM.R")
##Optimization
##First guess by MCMC 
BPM <- modMCMC(BObjectiveM, p = c(1, 25, 0.5, 1e-3, 0.001), 
               lower = c(1e-3, 0.1, 0.05, 1e-8, 0),
               upper = c(100, 3000, 1, 1, 1), niter = 30000)
##Estimate
ParmsBlag2014M <- abc_optim(fn = BObjectiveM, par = as.numeric(summary(BPM)["mean", ]), 
                             lb = as.numeric(summary(BPM)["min", ]), 
                             ub = as.numeric(summary(BPM)["max", ])) 
ParmsBlag2014M$par
#Goodness of fit and simulations
source("IndividualStudies/BFitM.R")
SimBlag2014M <- BFitM(ParmsBlag2014M$par)
SimBlag2014M$errors
write.csv(SimBlag2014M$Simulation, "../Manuscript/figure_data/SimBlagM.csv", row.names = F)
##Pirt model
source_python("../PythonScripts/IndividualStudies/Blag2014ODESolvP.py")
source_python("../PythonScripts/IndividualStudies/Pirt.py")
source("./IndividualStudies/BObjectiveP.R")
##Optimization
##First guess by MCMC 
BPP <- modMCMC(BObjectiveP, p = c(1, 25, 0.5, 1e-3, 1e-3, 0.001), 
               lower = c(1e-3, 0.1, 0.05, 1e-8, 1e-8, 0),
               upper = c(100, 3000, 1, 1, 1, 1), niter = 30000)
##Estimate
ParmsBlag2014P <- abc_optim(fn = BObjectiveP, par = as.numeric(summary(BPP)["mean", ]), 
                             lb = as.numeric(summary(BPP)["min", ]), 
                             ub = as.numeric(summary(BPP)["max", ])) 
ParmsBlag2014P$par
#Goodness of fit and simulations
source("IndividualStudies/BFitP.R")
SimBlag2014P <- BFitP(ParmsBlag2014P$par)
SimBlag2014P$errors
write.csv(SimBlag2014P$Simulation, "../Manuscript/figure_data/SimBlagP.csv", row.names = F)
#No. of observations
nt = SimBlag2014P$errors[6]
#DEB model
DEBSSRes = SimBlag2014$errors[5]
DEBpar = SimBlag2014$errors[7]
#Monod model
MSSRes = SimBlag2014M$errors[5]
Mpar = SimBlag2014M$errors[7]
#Pirt model
PSSRes = SimBlag2014P$errors[5]#Weighted residual sum of squares from python (all variables standardized to 0 mean and unit variance)
Ppar = SimBlag2014P$errors[7]

#DEB vs Monod
####F value =  (M1ss - M2ss)*(nt - M2p)/M2ss/(M2p - M1p)
(MSSRes - DEBSSRes)*(nt - DEBpar)/DEBSSRes/(DEBpar - Mpar)

####associated p value
pf(q=(MSSRes - DEBSSRes)*(nt - DEBpar)/DEBSSRes/(DEBpar - Mpar), 
   df1=(DEBpar - Mpar), 
   df2=(nt - DEBpar), 
   lower.tail=F)

#DEB vs Pirt
####F value =  (M1ss - M2ss)*(nt - M2p)/M2ss/(M2p - M1p)
(PSSRes - DEBSSRes)*(nt - DEBpar)/DEBSSRes/(DEBpar - Ppar)

####associated p value
pf(q=(PSSRes - DEBSSRes)*(nt - DEBpar)/DEBSSRes/(DEBpar - Ppar), 
   df1=(DEBpar - Ppar), 
   df2=(nt - DEBpar), 
   lower.tail=F)

##Visualizing the fit
##CO2
ggplot(BData, aes(Time, CO2)) + geom_point(cex=6, pch=21, fill = "grey") + facet_grid(.~Treatment) +
  theme_min + ylab(expression(paste("Cumulative respoiration (", mu, "mol(C) ", g^{-1}, ")"))) + xlab("Time (days)") +
  geom_line(data = (SimBlag2014$Simulation), aes(Time, CO2), lwd = 1.5) +
  geom_line(data = (SimBlag2014M$Simulation), aes(Time, CO2), lwd = 1.5, color = "grey") +
  geom_line(data = (SimBlag2014P$Simulation), aes(Time, CO2), lwd = 1.5, color = "grey30", lty = 2)
##DNA
ggplot(BData, aes(Time, DNA)) + geom_point(cex=6, pch=21, fill = "grey") + facet_grid(.~Treatment) +
  theme_min + ylab(expression(paste("DNA (nmol(C) ", g^{-1}, ")"))) + xlab("Time (days)") +
  geom_line(data = (SimBlag2014$Simulation), aes(Time, DNA), lwd = 1.5) +
  geom_line(data = (SimBlag2014M$Simulation), aes(Time, DNA), lwd = 1.5, color = "grey") +
  geom_line(data = (SimBlag2014P$Simulation), aes(Time, DNA), lwd = 1.5, color = "grey30", lty = 2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Global optimization~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#==============Model parameters and predictions are derived across subset of studies analyzed above
#==================================================#
##Parameters, which will be optimized (Initial guess, lower and upper bound) 
Im0 = c(mean(parsAll$Im), min(parsAll$Im), max(parsAll$Im)) #1
Km0 = c(mean(parsAll$Km), min(parsAll$Km), max(parsAll$Km)) #2
yA0 = c(mean(parsAll$yA[-10]), min(parsAll$yA), max(parsAll$yA)) #3
Em0 = c(mean(parsAll$Em), min(parsAll$Em), max(parsAll$Em)) #4
m0 = c(mean(parsAll$m), min(parsAll$m), max(parsAll$m))     #5
g0 = c(mean(parsAll$g), min(parsAll$g), max(parsAll$g))     #6
ne0 = c(mean(parsAll$ne, na.rm = T), min(parsAll$ne, na.rm = T), max(parsAll$ne, na.rm = T)) #7
nX10 = c(mean(parsAll$nX1, na.rm = T), min(parsAll$nX1, na.rm = T), max(parsAll$nX1, na.rm = T)) #8
iX10 = c(mean(parsAll$iX1, na.rm = T), 1e-4, 0.1)           #9
te0 = c(mean(parsAll$te, na.rm = T), min(parsAll$te, na.rm = T), max(parsAll$te, na.rm = T)) #10
tX10 = c(mean(parsAll$tX1, na.rm = T), min(parsAll$tX1, na.rm = T), max(parsAll$tX1, na.rm = T)) #11
#le0 = c(mean(parsAll$le, na.rm = T), 0, 0.1)                
#lX10 = c(mean(parsAll$lX1, na.rm = T), 0, 0.1)             
#GlE0 = c(0.5, 0, 1)                                         
#HPE0 = c(0.5, 0, 1)                                         #12
#BKE0 = c(0.5, 0, 1)                                         #13
#NE0 = c(0.5, 0, 1)                                         #14
#ZE0 = c(0.5, 0, 1)                                         
#Me0 = c(100, 1, 1000)                                       
#MX10 = c(100, 1, 1000)                                      
#BlE0 = c(0.5, 0, 1)                                         

#============================================= Assuming E0 = 0 ===============================================#
                
Parms = rbind(Im0, Km0, yA0, Em0, m0, g0, ne0, nX10, 
               iX10, te0, tX10)
#==================================================#
##Uploading respective functions
source_python("../PythonScripts/IndividualStudies/DEBmodelIso.py")
source_python("../PythonScripts/IndividualStudies/DEBmodel.py")
source_python("../PythonScripts/IndividualStudies/Glanville2016ODESolv.py")
source("GlobalFit/GlO.R")
source_python("../PythonScripts/IndividualStudies/Sparling1990ODESolv.py")
source("GlobalFit/HPO.R")
source("GlobalFit/SpO.R")
source_python("../PythonScripts/IndividualStudies/Bremer1990ODESolv.py") 
source("GlobalFit/BKO.R")
source_python("../PythonScripts/IndividualStudies/Marstorp1999ODESolv.py") 
source("GlobalFit/MO.R")
source_python("../PythonScripts/IndividualStudies/Tsai1997ODESolv.py") 
source("GlobalFit/TO.R")
source_python("../PythonScripts/IndividualStudies/Joerg2002ODESolv.py")
source("GlobalFit/JO.R")
source_python("../PythonScripts/IndividualStudies/Nanni1977ODESolvATP.py") 
source("GlobalFit/NO.R")
source_python("../PythonScripts/IndividualStudies/Blag2014ODESolv.py") 
source("GlobalFit/BO.R")
#==================================================#
##Objective function
GlobalFit <- function(x){
  return(sum(
    BO(x), GlO(x), HPO(x), MO(x), TO(x), JO(x), NO(x) #  ZO(x), SpO2(x), BKO(x),
  ))
}
#==================================================#
##Optimization
###First guess by MCMC 
GF <- modMCMC(GlobalFit, p = Parms[, 1], 
              lower = Parms[, 2],
              upper = Parms[, 3], updatecov = 100, burninlength = 1000, niter = 100000)
##Estimate
###Artificial bee-colony algorithm
GlobalParmsABC <- abc_optim(fn = GlobalFit, par = as.numeric(summary(GF)["mean", ]), 
                            lb = as.numeric(summary(GF)["min", ]), 
                            ub = as.numeric(summary(GF)["max", ])) 
GlobalParmsABC$par
#==================================================#
##Visualization
###Uploading respective functions
source("GlobalFit/BKFit.R")
source("GlobalFit/HPFit.R")
source("GlobalFit/GlFit.R")
source("GlobalFit/NFitG.R")
source("GlobalFit/JFitG.R")
source("GlobalFit/TFitG.R")
source("GlobalFit/MFitG.R")
source("GlobalFit/SpFit.R")
source("GlobalFit/BFit.R")
##Fit function
GlobalFitViz <- function(x){
  return(rbind(
    BFit(x), GlFit(x),BKFit(x), HPFit(x), NFitG(x), TFitG(x), JFitG(x), MFitG(x) #SpFit2(x), ZFit(x),  
  ))
}
VizOut <- GlobalFitViz(GlobalParmsABC$par)

#Calculate R squared
GlobalGoddness <- VizOut %>% group_by(Variable) %>% 
  summarize(R2 = 1-(sum((Observations-Predictions)^2, na.rm = T)/
              (sum((Observations-(mean(Observations, na.rm = T)))^2, na.rm = T))))
GlobalGoddness$R2adj <- with(GlobalGoddness, 1 - ((1 - R2)*((nrow(VizOut) - 1)/(nrow(VizOut) - 1 - 11))))

#Labels and annotation
VizLab <- c('CO2' = "CO[2]", 'DNA' = "DNA", '14CO2' = "atop(14, )~CO[2]", 'kec' = "italic(k[ec])",
            '14MBC' = "atop(14, )~MBC", '14S' = "atop(14, )~S", 'S' = "S", 'ATP' = "ATP",
            'Flush' = "CHCl[3]~flush")
VizAnot <- data.frame(Variable = as.character(GlobalGoddness$Variable), label = c("R^{2}==0.82", "R^{2}==0.93", "R^{2}==0.14",
                                                         "R^{2}==0.17", "R^{2}==0.61", "R^{2}==0.65",
                                                         "R^{2}==0.89", "R^{2}==0.50", "R^{2}==0.64"))

ggplot(VizOut, aes(Observations, Predictions)) + geom_point(cex = 6, aes(color = Study)) + theme_min + 
  geom_text(data = VizAnot, mapping = aes(x = -Inf, y = -Inf, label = label),
            hjust = -1.6, vjust = -1, parse = T, cex = 6) +
  facet_wrap(~Variable, scales = "free", labeller = as_labeller(VizLab, label_parsed)) + 
  geom_abline(intercept = 0, slope = 1) +
  theme()

#Parameters bounds
GFf <- modMCMC(GlobalFit, p = GlobalParmsABC$par, 
               lower = Parms2[, 2],
               upper = Parms2[, 3], updatecov = 100, burninlength = 1000, niter = 100000)
summary(GFf)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Herbert 1961~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
HEData <- read.csv("../SoilMBVariabilityData/Herbert1961.csv")
#Visualizing the data
##Cell mass
ggplot(HEData, aes(SGR, Mass)) + geom_point(cex = 6, pch = 21, fill = "grey") +
  theme_min + 
  xlab(expression(paste("Specific growth rate (", h^{-1}, ")"))) +
  ylab(expression(paste("Cell mass (picograms)")))
##DNA
ggplot(HEData, aes(SGR, DNA)) + geom_point(cex = 6, pch = 21, fill = "grey") +
  theme_min + 
  xlab(expression(paste("Specific growth rate (", h^{-1}, ")"))) +
  ylab(expression(paste("DNA abundance (unitless)")))
##DNA
ggplot(HEData, aes(SGR, RNA)) + geom_point(cex = 6, pch = 21, fill = "grey") +
  theme_min + 
  xlab(expression(paste("Specific growth rate (", h^{-1}, ")"))) +
  ylab(expression(paste("RNA abundance (unitless)")))
##Proteins
ggplot(HEData, aes(SGR, Protein)) + geom_point(cex = 6, pch = 21, fill = "grey") +
  theme_min + 
  xlab(expression(paste("Specific growth rate (", h^{-1}, ")"))) +
  ylab(expression(paste("Proteins abundance (unitless)")))

#==================================================================================================
##Parameters estimation
###Mass = X1*(1 + (SGR + SGR*g + m)/(z - SGR)); z = Im*yA/Em
###DNA = iX1/(1 + (SGR + SGR*g + m)/(z - SGR))
###RNA = (rX1 + re*F)/(1 + F); F = (SGR + SGR*g + m)/(z - SGR)
###Protein = (pX1 + pe*F)/(1 + F); F = (SGR + SGR*g + m)/(z - SGR)
###Coefficients: X1, z, g, m, iX1, rX1, re, pX1, pe
#==================================================================================================
###Starting with the mass
HEmass <- function(...){
  objective<-function(x){
   return(sum((HEData$Mass - x[1]*(1 + (HEData$SGR + HEData$SGR*x[2] + x[3])/(x[4] - HEData$SGR)))^2, na.rm = T))
  }
  ###Parameters estimation using MCMC
  mcmc<-modMCMC(f=objective, p=c(1, 0.1, 1e-5, 0.8), 
                lower=c(1e-3, 1e-3, 1e-15, 1e-3), 
                upper=c(1e6, 10, 1, 10), niter=100000)
  ###Improving the estimation using differential evolution algorithm
  # parDE<-DEoptim(fn=objective, lower=as.numeric(summary(mcmc)["min",]), upper=as.numeric(summary(mcmc)["max",]), 
  #                    control = c(itermax = 10000, steptol = 50, reltol = 1e-8, trace=FALSE, strategy=3, NP=250))
  # return(parDE$optim$bestmem)
  ##or artificial bee colony algorithm
  parABC<-abc_optim(fn=objective, par=as.numeric(summary(mcmc)["mean",]), 
                    lb=as.numeric(summary(mcmc)["min",]), 
                    ub=as.numeric(summary(mcmc)["max",]), maxCycle = 1e6)
  return(parABC$par)
}

HEmasspar <- HEmass()

##Visualize predictions
HEPreds <- data.frame(SGR = seq(0, max(HEData$SGR), length.out = 50))
HEPreds$Mass <- with(HEPreds, HEmasspar[1]*(1 + (SGR + SGR*HEmasspar[2] + HEmasspar[3])/(HEmasspar[4] - SGR)))

##Cell mass
ggplot(HEData, aes(SGR, Mass)) + geom_point(cex = 6, pch = 21, fill = "grey") +
  theme_min + 
  xlab(expression(paste("Specific growth rate (", h^{-1}, ")"))) + ylim(0, 0.7) +
  ylab(expression(paste("Cell mass (picograms)"))) +
  geom_line(data = HEPreds, aes(SGR, Mass), lwd = 1.5)

#===========
# DNA
#===========
nlsDNA<-nls(DNA~iX1/(1 + (SGR + SGR*HEmasspar[2] + HEmasspar[3])/(HEmasspar[4] - SGR)),
            HEData, start = list(iX1 = 0.005))
summary(nlsDNA)
##Visualize predictions
HEPreds$DNA <- with(HEPreds, coef(nlsDNA)/(1 + (SGR + SGR*HEmasspar[2] + HEmasspar[3])/(HEmasspar[4] - SGR)))

##DNA
ggplot(HEData, aes(SGR, DNA)) + geom_point(cex = 6, pch = 21, fill = "grey") +
  theme_min + 
  xlab(expression(paste("Specific growth rate (", h^{-1}, ")"))) +
  ylab(expression(paste("DNA abundance (unitless)"))) + 
  geom_line(data = HEPreds, aes(SGR, DNA), lwd = 1.5) 

#===========
# RNA
#===========
nlsRNA<-nls(RNA~(rX1 + re*(SGR + SGR*HEmasspar[2] + HEmasspar[3])/(HEmasspar[4] - SGR))/
              (1 +  (SGR + SGR*HEmasspar[2] + HEmasspar[3])/(HEmasspar[4] - SGR)),
            HEData, start = list(rX1 = 0.1, re = 0.4))
summary(nlsRNA)
##Visualize predictions
HEPreds$RNA <- with(HEPreds, (coef(nlsRNA)[1] + coef(nlsRNA)[2]*(SGR + SGR*HEmasspar[2] + HEmasspar[3])/(HEmasspar[4] - SGR))/
                      (1 +  (SGR + SGR*HEmasspar[2] + HEmasspar[3])/(HEmasspar[4] - SGR)))

##RNA
ggplot(HEData, aes(SGR, RNA)) + geom_point(cex = 6, pch = 21, fill = "grey") +
  theme_min + 
  xlab(expression(paste("Specific growth rate (", h^{-1}, ")"))) +
  ylab(expression(paste("RNA abundance (unitless)"))) + 
  geom_line(data = HEPreds, aes(SGR, RNA), lwd = 1.5)

#===========
# Proteins
#===========
nlsP<-nls(Protein~(pX1 + pe*(SGR + SGR*HEmasspar[2] + HEmasspar[3])/(HEmasspar[4] - SGR))/
            (1 +  (SGR + SGR*HEmasspar[2] + HEmasspar[3])/(HEmasspar[4] - SGR)),
            HEData, start = list(pX1 = 0.8, pe = 0.1))
summary(nlsP)
##Visualize predictions
HEPreds$Protein <- with(HEPreds, (coef(nlsP)[1] + coef(nlsP)[2]*(SGR + SGR*HEmasspar[2] + HEmasspar[3])/(HEmasspar[4] - SGR))/
                          (1 +  (SGR + SGR*HEmasspar[2] + HEmasspar[3])/(HEmasspar[4] - SGR)))

##Protein
ggplot(HEData, aes(SGR, Protein)) + geom_point(cex = 6, pch = 21, fill = "grey") +
  theme_min + 
  xlab(expression(paste("Specific growth rate (", h^{-1}, ")"))) +
  ylab(expression(paste("Proteins abundance (unitless)"))) + 
  geom_line(data = HEPreds, aes(SGR, Protein), lwd = 1.5)

#=================
# All in one graph
#=================
ggplot(HEData, aes(x = SGR)) + geom_point(cex = 6, pch = 21, aes(y = Mass, fill = "Cell mass")) +
  geom_point(cex = 6, pch = 21, aes(y = RNA*5, fill = "RNA")) + theme_min +
  geom_point(cex = 6, pch = 21, aes(y = Protein, fill = "Proteins")) +
  geom_point(cex = 6, pch = 21, aes(y = DNA*5, fill = "DNA")) +
  scale_y_continuous(sec.axis = sec_axis(~./5, name = "RNA and DNA abundance")) +
  labs(y = "Cell Mass (mg) and abundance of Proteins",
       x = expression(paste("Specific growth rate (", h^{-1}, ")"))) +
  scale_fill_manual(values = c("grey60", "black", "white", "grey90")) + 
  xlim(0, 1) + theme(legend.title = element_blank(), legend.position = c(0.2, 0.85), 
                     legend.background = element_blank()) +
  geom_line(data = HEPreds, aes(SGR, Mass), lwd = 1, color = "grey60") +
  geom_line(data = HEPreds, aes(SGR, DNA*5), lwd = 1, color = "black") +
  geom_line(data = HEPreds, aes(SGR, RNA*5), lwd = 1, color = "grey") +
  geom_line(data = HEPreds, aes(SGR, Protein), lwd = 1, color = "black") +
  ggtitle("A)")

#==================Growth rate - E/Em relationship
#Parameters
p <- GlobalParmsABC$par
names(p) <- c("Im", "Km", "yA", "Em", "m", "g", "ne", "nX1", "iX1", "te", "tX1") 
#=======
EEm <- seq(0, 1, 0.05)
SGR <- (p[1]*p[3]*EEm-p[5])/(1 + p[6] + EEm*p[4])

ggplot(data.frame(EEm, SGR), aes(EEm, SGR)) + geom_line(lwd = 0.9) + theme_min +
  ylab(expression(paste("SGR (", day^{-1}, ")"))) + xlab(expression(paste(frac(E, E[m])))) + 
  ggtitle("A)")

#==================kec/kdna - E/Em relationship
kec <- (p[7]*EEm*p[4] + p[8])/(1 + EEm*p[4])
kDNA <- (p[9])/(1 + EEm*p[4])

ggplot(data.frame(EEm, kec, kDNA), aes(x = EEm)) + geom_line(lwd = 0.9, aes(y = kec, color = "kec")) + 
  geom_line(lwd = 0.9, aes(y = kDNA*5, color = "kDNA")) + theme_min +
  scale_y_continuous(sec.axis = sec_axis(~./5, 
                                         name = expression(paste(italic(k[DNA]))))) + 
  ylab(expression(paste(italic(k[ec])))) + xlab(expression(paste(frac(E, E[m])))) + 
  ggtitle("B)") + theme(legend.title = element_blank(), legend.text.align = 0, legend.position = c(0.7, 0.6)) + 
  scale_color_manual(values = c("grey60", "black"), labels = c(expression(italic(k[DNA])), expression(italic(k[ec]))))

#==================kec - E/Em theoretical relationship
##Six different combinations of nX1 and ne
kecT <- data.frame(EEm = rep(seq(0, 1, 0.05), 6),
                   Legend = rep(c("A", "B", "C", "D", "E", "F"), each = length(seq(0, 1, 0.05))),
                   kec = c((0.8*seq(0, 1, 0.05)*p[4] + 0.3)/(1 + seq(0, 1, 0.05)*p[4]),
                           (0.8*seq(0, 1, 0.05)*p[4] + 0.1)/(1 + seq(0, 1, 0.05)*p[4]),
                           (0.6*seq(0, 1, 0.05)*p[4] + 0.3)/(1 + seq(0, 1, 0.05)*p[4]),
                           (0.6*seq(0, 1, 0.05)*p[4] + 0.1)/(1 + seq(0, 1, 0.05)*p[4]),
                           (0.4*seq(0, 1, 0.05)*p[4] + 0.3)/(1 + seq(0, 1, 0.05)*p[4]),
                           (0.4*seq(0, 1, 0.05)*p[4] + 0.1)/(1 + seq(0, 1, 0.05)*p[4])))

ggplot(kecT, aes(EEm, kec)) + geom_line(aes(color = Legend, linetype = Legend)) +
  theme_min + scale_color_manual(values = c("black", "black", "grey", "grey", "grey60", "grey60"),
                                 labels = c(expression(n[G]==0.8~~n[B]==0.3), 
                                            expression(n[G]==0.8~~n[B]==0.1),
                                            expression(n[G]==0.6~~n[B]==0.3), 
                                            expression(n[G]==0.6~~n[B]==0.1),
                                            expression(n[G]==0.4~~n[B]==0.3), 
                                            expression(n[G]==0.4~~n[B]==0.1)))+
  scale_linetype_manual(values = c("solid", "dashed", "solid", "dashed", "solid", "dashed"),
                     labels = c(expression(n[G]==0.8~~n[B]==0.3), 
                                expression(n[G]==0.8~~n[B]==0.1),
                                expression(n[G]==0.6~~n[B]==0.3), 
                                expression(n[G]==0.6~~n[B]==0.1),
                                expression(n[G]==0.4~~n[B]==0.3), 
                                expression(n[G]==0.4~~n[B]==0.1))) +
  xlab(expression(paste(frac(G, G[m])))) + ylab(expression(paste(italic(k[ec])))) +
  scale_y_continuous(limits = c(0.1, 0.65), breaks = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6)) +
  geom_hline(yintercept = 0.45, color = "red") + theme(legend.title = element_blank()) +
  ggtitle("C)")
  
#==================Composite figure to sho in main text
grid.arrange(
  ##Substrate concentration
  ggplot(HPData, aes(Time, S)) + geom_point(cex=6, pch=21, fill = "grey") +
    theme_min + ylab(expression(paste({}^{14},"S (", mu, "mol(C) " , g^{-1}, "(DW))"))) + 
    xlab("Time (days)") +
    geom_line(data = HPSim, aes(Time, S), lwd = 1.2) +
    ggtitle("A)"),
  ##Cmic 14
  ggplot(HPData, aes(Time, Flush)) + geom_point(cex=6, pch=21, fill = "grey") +
    theme_min + ylab(expression(paste({}^{14}," CLC (", mu, "mol(C) " , g^{-1}, "(DW))"))) + 
    xlab("Time (days)") +
    geom_line(data = HPSim, aes(Time, Flush), lwd = 1.2) +
    ggtitle("B)"),
  ##Cumulative respiration
  ggplot(HPData, aes(Time, CO2)) + geom_point(cex=6, pch=21, fill = "grey") +
    theme_min + ylab(expression(paste({}^{14},"C", O[2]," (", mu, "mol(C) " , g^{-1}, "(DW))"))) + 
    xlab("Time (days)") +
    geom_line(data = HPSim, aes(Time, CO2), lwd = 1.2) +
    ggtitle("C)"),
  ##kec factor
  ggplot(HPData, aes(Time, kec)) + geom_point(cex=6, pch=21, fill = "grey") +
    theme_min + ylab(expression(paste(italic(k[ec])))) + xlab("Time (days)") + #ylim(0.20, 0.35) +
    geom_line(data = HPSim, aes(Time, kec), lwd = 1.2) +
    scale_y_continuous(limits = c(0.2, 0.5)) +
    ggtitle("D)"),
  nrow = 2
)

