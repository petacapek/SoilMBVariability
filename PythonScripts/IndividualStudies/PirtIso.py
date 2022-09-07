def PirtIso (y, t, pars):
    #define initial pools
    Sl=y[0];    Bl=y[1];	  CO2l=y[2]; 
    Su=y[3];    Bu=y[4]; 
    #define parameters
    Im=pars[0]; 
    Km=pars[1];     
    yA=pars[2];
    d=pars[3];
    m=pars[4] 
    #Scaling function for substrate uptake
    f=Sl/(Km+Sl)
    
    #Isotope signals
    Satm = Sl/(Sl + Su)
    Batm = Bl/(Bl + Bu)
        
    #Sums of isotope pools
    S = Sl + Su
    B = Bl + Bu
    
    #Fluxes
    uptake=Im*B*f
    growth = uptake*yA
    death = B*d
    maintenance = B*m
    
    #Define derivatives
    dSldt = -uptake*Satm + death*Batm
    dBldt = growth*Satm - death*Batm - maintenance*Batm
    dCO2ldt = uptake*(1 - yA)*Satm + maintenance*Batm
    dSudt = -uptake*(1-Satm) + death*(1-Batm)
    dBudt = growth*(1-Satm) - death*(1-Batm) - maintenance*(1-Batm)
         
    return dSldt, dBldt, dCO2ldt, dSudt, dBudt;
