#########################################################################################################
##############################Soil microbial biomass composition variability#############################
###################Using Dynamic Energy Budget theory (DEB) as an explanatory concept####################
#########################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Libraries~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
library(dplyr)
library(ggplot2)
library(bbmle)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~GGplot theme~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
theme_min<-theme(axis.text.x=element_text(vjust=0.2, size=18, colour="black"),
                 axis.text.y=element_text(hjust=0.2, size=18, colour="black"),
                 axis.title=element_text(size=18, colour="black"),
                 axis.line=element_line(size=0.5, colour="black"),
                 strip.text=element_text(size=18, face="bold"),
                 axis.ticks=element_line(size=1, colour="black"),
                 axis.ticks.length=unit(-0.05, "cm"),
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
#########################################Glanville et al. (2016)##########################################
##Reading the data
GlanData = read.csv("/mnt/580CBE2464C5F83D/pracovni/data_statistika/SoilMBVariability/SoilMBVariabilityData/Glanville2016.csv")

##Plotting the data
ggplot(GlanData, aes(Time/24, kec_original)) + geom_point(cex=6, pch=21, aes(fill = Substrate)) +
  theme_min + xlab("Time (days)") + ylab(expression(paste(italic(k[ec])))) +
  theme(legend.position = c(0.8, 0.8)) + scale_fill_manual(values = c("grey", "black"))

##Estimating the parameters
##Scaled energy density (e) changes over time following the equation
###e(t) = e0*exp(-v*t); where t is time and v is energy conductance
##kec is then described by the equation 
###kec = (nX1*MX1 + ne*ce*e)/(MX1 + ce*e)
####ne is assumed to be 1 (all energy reserves are released by chloroform and extracted from soil),
####ce is assumed to be 4*MX1 (Hanegraaf and Muller, 2001)
##The fitted equation is thus
###kec = (nX1 + 4*e0*exp(-v*t))/(1 + 4*e0*exp(-v*t))
####3 parameters are fitted

####Equation fit
kecFit1 = nls(kec_original~(nX1 + 4*e0*exp(-v*Time))/(1 + 4*e0*exp(-v*Time)), 
              data = GlanData, start = list(nX1 = 0.1, v=0.3, e0 = 0.5))
summary(kecFit1)

##Visualizing the fit
ggplot(GlanData, aes(Time, kec_original)) + geom_point(cex=6, pch=21, aes(fill = Substrate)) +
  theme_min + xlab("Time (days)") + ylab(expression(paste(italic(k[ec])))) +
  theme(legend.position = c(0.8, 0.8)) + scale_fill_manual(values = c("grey", "black")) +
  stat_function(fun = function(x){(0.263059 + 4*0.046645*exp(-0.008484*x))/(1 + 4*0.046645*exp(-0.008484*x))},
                color="black")

####Let's assume that e0 is different for different substrates
kecFit2 = nls(kec_original~(nX1 + 4*e0[Substrate]*exp(-v*Time))/(1 + 4*e0[Substrate]*exp(-v*Time)), 
              data = GlanData, start = list(nX1 = 0.1, v=0.3, e0 = c(0.5, 0.5)))
summary(kecFit2)
anova(kecFit1, kecFit2)

##Visualizing the fit
ggplot(GlanData, aes(Time, kec_original)) + geom_point(cex=6, pch=21, aes(fill = Substrate)) +
  theme_min + xlab("Time (days)") + ylab(expression(paste(italic(k[ec])))) +
  theme(legend.position = c(0.8, 0.8)) + scale_fill_manual(values = c("grey", "black")) +
  stat_function(fun = function(x){(0.262187 + 4*0.036832*exp(-0.007665*x))/(1 + 4*0.036832*exp(-0.007665*x))},
                color="grey30", lwd = 1.2)+
  stat_function(fun = function(x){(0.262187 + 4*0.054898*exp(-0.007665*x))/(1 + 4*0.054898*exp(-0.007665*x))},
                color="black", lwd = 1.2)

#Comparing the model with 2-pool model
kecFit3 = nls(kec_original~X1[Substrate]*exp(-k1[Substrate]*Time) + X2[Substrate]*exp(-k2[Substrate]*Time), 
              data = GlanData, start = list(X1 = c(3.926e-01, 3.926e-01), 
                                            X2 = c(3.926e-01, 3.926e-01),
                                            k1 = c(3.470e-02, 3.470e-02), 
                                            k2 = c(1.962e-04, 1.962e-04)),
              control = c(maxiter = 1000))
summary(kecFit3)
anova(kecFit2, kecFit3)
AICtab(kecFit2, kecFit3, sort = T, weights = T, base = T, logLik = T)

##Visualizing the fit
ggplot(GlanData, aes(Time, kec_original)) + geom_point(cex=6, pch=21, aes(fill = Substrate)) +
  theme_min + xlab("Time (days)") + ylab(expression(paste(italic(k[ec])))) +
  theme(legend.position = c(0.8, 0.8)) + scale_fill_manual(values = c("grey", "black")) +
  stat_function(fun = function(x){(0.262187 + 4*0.036832*exp(-0.007665*x))/(1 + 4*0.036832*exp(-0.007665*x))},
                color="grey30", lwd = 1.2)+
  stat_function(fun = function(x){(0.262187 + 4*0.054898*exp(-0.007665*x))/(1 + 4*0.054898*exp(-0.007665*x))},
                color="black", lwd = 1.2)+
  stat_function(fun = function(x){6.835e-02*exp(-9.365e-02*x) + 3.190e-01*exp(-2.564e-04*x)},
                color="grey30", lwd = 1.2, lty = 2)+
  stat_function(fun = function(x){1.367e-01*exp(-7.628e-03*x) + 2.604e-01*exp(2.052e-05*x)},
                color="black", lwd = 1.2, lty = 2)

#########################################Marstorp et al. (1999)##########################################
##Reading the data
MarData = read.csv("/mnt/580CBE2464C5F83D/pracovni/data_statistika/SoilMBVariability/SoilMBVariabilityData/Marstorp1999.csv")

##Estimating the parameters
##In this case, more complicated equation have to be solved
###de/dt = v*f - v*e