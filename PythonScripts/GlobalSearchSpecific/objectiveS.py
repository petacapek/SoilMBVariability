import numpy as np
import pandas as pd
from scipy.integrate import odeint
from calcDEB_S import calcDEB_S
from DEBmodelIso import DEBmodelIso

def objectiveS (x):
    #define parameters
    ##yA, Km, v, m, g, ce, nX1, eu_i = x[16]
    pars = x[0:7]
    
    #read data
    d = pd.read_csv('/mnt/580CBE2464C5F83D/pracovni/data_statistika/SoilMBVariability/SoilMBVariabilityData/Santruckova2004.csv', sep=',')

    #initial conditions
    Sl_i = d.Sinit[0]
    #eu_i = x[16]  
    X1u_i = d.Cmicinit[0]/(pars[5]*0.25*pars[6])
        
    y0 = np.array([Sl_i, 0, 0, 0, X1u_i])
    #times
    t = d.Time
    #model simulations
    yhat_full = calcDEB_S(DEBmodelIso, pars, t, y0)
    #observations
    obs=np.concatenate((np.array([d.S]).reshape(len(d.Time),1),
                        np.array([d.CO214cumul]).reshape(len(d.Time),1),
                        np.array([d.kec]).reshape(len(d.Time),1),
                        np.array([d.Cmic14]).reshape(len(d.Time),1)), axis=1)
    #weights
    weights=np.concatenate((np.nanstd(d.S).repeat(len(d.Time)).reshape(len(d.Time),1),
                            np.nanstd(d.CO214cumul).repeat(len(d.Time)).reshape(len(d.Time),1),
                            np.nanstd(d.kec).repeat(len(d.Time)).reshape(len(d.Time),1),
                            np.nanstd((d.Cmic14)).repeat(len(d.Time)).reshape(len(d.Time),1)),
                       axis=1)
    out=np.nansum(((yhat_full-obs)/weights)**2)
    return out
