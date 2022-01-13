def DEBmodelIso (y, t, pars):
    #define initial pools (eu is assumed to be zero)
    Sl=y[0];    el=y[1];    X1l=y[2];     CO2l=y[3];
    eu=y[4];	X1u=y[5]
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
    eatm = el/(el + eu)
    X1atm = X1l/(X1l + X1u)
    
    #Fluxes
    uptake=(v*ce/yA)*(X1l+X1u)*f #labelled substrate only
    growth = (v*(el + eu)-m*g)/(el + eu + g)
    
    #Define derivatives
    ##Labelled pools
    dSldt = -uptake
    deldt = v*f - v*(el + eu)*eatm
    dX1ldt = max(0, (X1l + X1u)*growth*eatm) + min(0, (X1l + X1u)*growth*X1atm) 
    dCO2ldt = uptake*(1 - yA) + ce*((X1l + X1u)*(el + eu)*(v-growth)*eatm) - max(0, (X1l + X1u)*growth*eatm) - min(0, (X1l + X1u)*growth*X1atm) 
    ##Unabelled pools
    #dSudt
    deudt = - v*(el + eu)*(1-eatm)
    dX1udt = max(0, (X1l + X1u)*growth*(1 - eatm)) + min(0, (X1l + X1u)*growth*(1 - X1atm))
        
    return dSldt, deldt, dX1ldt, dCO2ldt, deudt, dX1udt;
