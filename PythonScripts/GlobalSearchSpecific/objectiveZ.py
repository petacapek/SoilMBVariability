import numpy as np
import pandas as pd
from scipy.integrate import odeint
from calcDEB_Z import calcDEB_Z
from DEBmodel import DEBmodel

def objectiveZ (x):
    #define parameters
    ##yA, Km, v, m, g, ce, le, lX1
    x0 = np.array(x)
    pars = x0[[0, 1, 2, 3, 4, 5, 10, 11]]
    
    #reading data
    d = pd.read_csv('/mnt/580CBE2464C5F83D/pracovni/data_statistika/SoilMBVariability/SoilMBVariabilityData/Ziegler2005.csv', sep=',')

    #initial conditions
    S_i = d.Sinit[0]
            
    e_i = x[16]
    X1_i = d.PLFAinit[0]/(pars[5]*(0.25*pars[7] + pars[6]*e_i))
    
    y0 = np.array([S_i, e_i, X1_i, 0])
    #times
    t = d.Time

    #model simulations
    yhat_full = calcDEB_Z(DEBmodel, pars, t, y0)
    
     
    #observations
    obs=np.concatenate((np.array([d.S]).reshape(len(d.Time),1),
                        np.array([d.CO2cumul]).reshape(len(d.Time),1),
                        np.array([d.PLFA]).reshape(len(d.Time),1)),
                     axis=1)

    #weights
    weights=np.concatenate((np.nanstd(d.S).repeat(len(d.Time)).reshape(len(d.Time),1),
                            np.nanstd(d.CO2cumul).repeat(len(d.Time)).reshape(len(d.Time),1),
                            np.nanstd((d.PLFA)).repeat(len(d.Time)).reshape(len(d.Time),1)),
                       axis=1)
    
    out=np.nansum(((yhat_full-obs)/weights)**2)

    return out
