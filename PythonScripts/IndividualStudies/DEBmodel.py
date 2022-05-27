def DEBmodel (y, t, pars):
    #define initial pools
    S=y[0];    E=y[1];    X1=y[2];  CO2=y[3];  
    #define parameters
    Im=pars[0]; 
    Km=pars[1];     
    yA=pars[2];
    Em=pars[3]; 
    m=pars[4]; 
    g=pars[5];
    #Scaling function for substrate uptake
    f=S/(Km+S)
    
    #Fluxes
    uptake=Im*X1*f
    assimilation = Im*yA*f #X1 specific
    mobilization = Im*yA*E/Em
    growth = (mobilization - m)/(1 + g + E) 
    
    #Define derivatives
    dSdt = -uptake
    dEdt = assimilation - mobilization
    dX1dt = X1*(max(0, growth) + min(0, growth*(1+g)))
    dCO2dt = uptake*(1 - yA) + X1*(max(growth*g, 0) + m)
         
    return dSdt, dEdt, dX1dt, dCO2dt;
