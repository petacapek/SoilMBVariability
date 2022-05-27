def DEBmodelIso (y, t, pars):
    #define initial pools
    Sl=y[0];    El=y[1];    X1l=y[2];     
    CO2l=y[3];  Eu=y[4];    X1u=y[5];
    #define parameters
    Im=pars[0]; 
    Km=pars[1];     
    yA=pars[2];
    Em=pars[3]; 
    m=pars[4]; 
    g=pars[5];
    #Scaling function for substrate uptake
    f=Sl/(Km+Sl) #labelled substrate only
    
    #Isotope signals
    Eatm = El/(El + Eu)
    X1atm = X1l/(X1l + X1u)
    
    #Sums of isotope pools
    X1 = X1l + X1u
    E = El + Eu
    
    #Fluxes
    uptake=Im*X1*f #labelled substrate only
    assimilation = Im*yA*f #X1 specific
    mobilization = Im*yA*E/Em
    growth = (mobilization - m)/(1 + g + E) 
    
    #Define derivatives
    ##Labelled pools
    dSldt = -uptake
    dEldt = assimilation - mobilization*Eatm
    dX1ldt = X1*(max(0, growth*Eatm) + min(0, growth*(1+g)*X1atm))
    dCO2ldt = uptake*(1 - yA) + X1*(max(growth*g*Eatm, 0) - min(0, growth*(1+g)*X1atm) + (m + min(0, growth*(1+g)*X1atm))*Eatm)
    ##Unabelled pools
    #dSudt
    dEudt = - mobilization*(1 - Eatm)
    dX1udt = X1*(max(0, growth*(1 - Eatm)) + min(0, growth*(1+g)*(1 - X1atm))) 
        
    return dSldt, dEldt, dX1ldt, dCO2ldt, dEudt, dX1udt;
