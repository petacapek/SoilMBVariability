def DEBmodel (y, t, pars):
    #define initial pools (eu is assumed to be zero)
    S=y[0];    e=y[1];    X1=y[2];     CO2=y[3];
    
    #define parameters
    yA=pars[0]; 
    Km=pars[1];     
    v=pars[2];
    m=pars[3]; 
    g=pars[4]; 
    #k=pars[4];
    ce=pars[5];
    MX1=ce/4;
    
    #Scaling function for substrate uptake
    f=S/(Km+S)
    #Fluxes
    uptake=(v*ce/yA)*X1*f
    growth = (v*e-m*g)/(e + g)
    
    #Define derivatives
    dSdt = -uptake
    dedt = v*(f - e)
    dX1dt = growth*X1 #- k*X1
    dCO2dt = uptake*(1 - yA) + ce*(X1*e*(v-growth)) - growth*X1*MX1 
           
    return dSdt, dedt, dX1dt, dCO2dt;
