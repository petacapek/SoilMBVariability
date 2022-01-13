def DEBmodelIso (y, t, pars):
    #define initial pools (eu is assumed to be zero)
    Sl=y[0];    el=y[1];    X1l=y[2];     CO2l=y[3];
    X1u=y[4]
    #define parameters
    yA=pars[0]; 
    Km=pars[1];     
    v=pars[2];
    m=pars[3]; 
    g=pars[4]; 
    ce=pars[5];
    MX1=ce/4;
    #Scaling function for substrate uptake
    f=Sl/(Km+Sl) #labelled substrate only
    
    #Isotope signals
    eatm = 1
    X1atm = X1l/(X1l + X1u)
    
    #Fluxes
    uptake=(v*ce/yA)*(X1l+X1u)*f #labelled substrate only
    growth = (v*el-m*g)/(el + g)
    
    #Define derivatives
    ##Labelled pools
    dSldt = -uptake
    deldt = v*(f - el)
    dX1ldt = max(0, (X1l + X1u)*growth) + min(0, (X1l + X1u)*growth*X1atm) 
    dCO2ldt = uptake*(1 - yA) + ce*((X1l + X1u)*el*(v-growth)) - max(0, (X1l + X1u)*growth)*MX1 - min(0, (X1l + X1u)*growth*X1atm)*MX1 
    ##Unabelled pools
    #dSudt
    #deudt = - v*eu
    dX1udt = min(0, (X1l + X1u)*growth*(1-X1atm)) 
        
    return dSldt, deldt, dX1ldt, dCO2ldt, dX1udt;
