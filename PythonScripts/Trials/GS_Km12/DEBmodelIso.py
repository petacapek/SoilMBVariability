def DEBmodelIso (y, t, pars):
    #define initial pools (eu is assumed to be zero)
    Sl=y[0];    el=y[1];    X1l=y[2];     CO2l=y[3];
    X1u=y[4]
    #define parameters
    yA=pars[0]; 
    Km1=pars[1];   
    Km2=pars[2];
    v=pars[3]; 
    m=pars[4]; 
    g=pars[5]; 
    #k=pars[4];
    ce=pars[6];
    MX1=ce/4;
    
    #Scaling function for substrate uptake
    f1=Sl/(Km1+Sl) #labelled substrate only
    f2=Sl/(Km2+Sl)
    
    #Isotope signals
    eatm = 1
    X1atm = X1l/(X1l + X1u)
    
    #Fluxes
    uptake1=(v*ce/yA)*(X1l+X1u)*f1 #labelled substrate only
    uptake2=(v*ce)*(X1l+X1u)*f2 #labelled substrate only
    growth = (v*el-m*g)/(el + g)
    
    #Define derivatives
    ##Labelled pools
    dSldt = -uptake1 - uptake2
    deldt = v*(f1 + f2) - v*el
    dX1ldt = max(0, (X1l + X1u)*growth) + min(0, (X1l + X1u)*growth*X1atm) 
    dCO2ldt = uptake1*(1 - yA) + ce*((X1l + X1u)*el*(v - growth)) - max(0, (X1l + X1u)*growth)*MX1 - min(0, (X1l + X1u)*growth*X1atm)*MX1 
    ##Unabelled pools
    #dSudt
    #deudt = - v*eu
    dX1udt = min(0, (X1l + X1u)*growth*(1-X1atm)) 
        
    return dSldt, deldt, dX1ldt, dCO2ldt, dX1udt;
