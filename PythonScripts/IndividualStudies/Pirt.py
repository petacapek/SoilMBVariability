def Pirt (y, t, pars):
    #define initial pools
    S=y[0];    B=y[1];	  CO2=y[2];  
    #define parameters
    Im=pars[0]; 
    Km=pars[1];     
    yA=pars[2];
    d=pars[3];
    m=pars[4] 
    #Scaling function for substrate uptake
    f=S/(Km+S)
    
    #Fluxes
    uptake=Im*B*f
    growth = uptake*yA
    death = B*d
    maintenance = B*m
    
    #Define derivatives
    dSdt = -uptake + death
    dBdt = growth - death - maintenance
    dCO2dt = uptake*(1 - yA) + maintenance
         
    return dSdt, dBdt, dCO2dt;
