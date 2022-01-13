import numpy as np
import pandas as pd
from scipy.integrate import odeint
from calcDEB_N import calcDEB_N
from DEBmodel import DEBmodel

def objectiveN_A (x):
    #define parameters
    ##yA, Km, v, m, g, ce, te, tX1, MwX1, Mwe
    x0 = np.array(x)
    pars = x0[[0, 1, 2, 3, 4, 5, 8, 9, 10, 11]]
    
    #reading data for treatment A
    dall = pd.read_csv('/mnt/580CBE2464C5F83D/pracovni/data_statistika/SoilMBVariability/SoilMBVariabilityData/Nannipieri1977.csv', sep=',')
    d = dall[dall.Treatment == "A"]

    #initial conditions
    S_i = d.Sinit[0]
            
    e_i = 0.25*(pars[7] - d.ATPinit[0]/d.Winit[0]*pars[8])/(d.ATPinit[0]/d.Winit[0]*(pars[9] - pars[6]))
    X1_i = d.Winit[0]/(pars[5]*(0.25*pars[8] + pars[9]*e_i))
    
    y0 = np.array([S_i, e_i, X1_i, 0])
    #times
    t = d.Time

    #model simulations
    yhat_full = calcDEB_N(DEBmodel, pars, t, y0)
    
     
    #observations
    obs=np.concatenate((np.array([d.CO2]).reshape(len(d.Time),1),
                        np.array([d.W]).reshape(len(d.Time),1),
                        np.array([d.ATP]).reshape(len(d.Time),1)),
                     axis=1)

    #weights
    weights=np.concatenate((np.nanstd(d.CO2).repeat(len(d.Time)).reshape(len(d.Time),1),
                            np.nanstd((d.W)).repeat(len(d.Time)).reshape(len(d.Time),1),
                            np.nanstd((d.ATP)).repeat(len(d.Time)).reshape(len(d.Time),1)),
                       axis=1)
    
    out=np.nansum(((yhat_full-obs)/weights)**2)

    return out
