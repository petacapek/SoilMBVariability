def DEBmodel (y, t, pars):
    #define initial pools (eu is assumed to be zero)
    S=y[0];    e=y[1];    X1=y[2];     CO2=y[3];
    
    #define parameters
    yA=pars[0]; 
    Km=pars[1];   
    v1=pars[2];
    v2=pars[3]; 
    m=pars[4]; 
    g=pars[5]; 
    #k=pars[4];
    ce=pars[6];
    MX1=ce/4;
    
    #Scaling function for substrate uptake
    f1=S/(Km+S)
    f2=S/(Km/500+S)
    #Fluxes
    uptake1=(v1*ce/yA)*X1*f1
    uptake2=(v2*ce)*X1*f2
    growth = ((v1+v2)*e-m*g)/(e + g)
    
    #Define derivatives
    dSdt = -uptake1-uptake2
    dedt = v1*f1 + v2*f2 - (v1+v2)*e
    dX1dt = growth*X1 #- k*X1
    dCO2dt = uptake1*(1 - yA) + ce*(X1*e*(v1+v2-growth)) - growth*X1*MX1 
           
    return dSdt, dedt, dX1dt, dCO2dt;
