import numpy as np
import pandas as pd
from scipy.integrate import odeint
from calcDEB_J import calcDEB_J
from DEBmodel import DEBmodel

def objectiveJ (x):
    #define parameters
    ##yA, Km, v, m, g, ce, nX1, te, tX1
    x0 = np.array(x)
    pars = x0[[0,1,2,3,4,5,6, 8, 9]]
    
    #read data
    d = pd.read_csv('/mnt/580CBE2464C5F83D/pracovni/data_statistika/SoilMBVariability/SoilMBVariabilityData/Joergensen2002.csv', sep=',')

    #initial conditions
    S_i = d.Sinit[0]
            
    e_i = 0.25*((d.ATPinit[0]/d.Cmicinit[0])*pars[6] - pars[8])/(pars[7] - (d.ATPinit[0]/d.Cmicinit[0]))
    X1_i = d.Cmicinit[0]/((0.25*pars[6] + e_i)*pars[5])
    
    y0 = np.array([S_i, e_i, X1_i, 0])
    #times
    t = d.Time

    #model simulations
    yhat_full = calcDEB_J(DEBmodel, pars, t, y0)
    
     
    #observations
    obs=np.concatenate((np.array([d.S]).reshape(len(d.Time),1),
                        np.array([d.Cmic]).reshape(len(d.Time),1),
                        np.array([d.ATP]).reshape(len(d.Time),1)),
                     axis=1)

    #weights
    weights=np.concatenate((np.nanstd(d.S).repeat(len(d.Time)).reshape(len(d.Time),1),
                            np.nanstd((d.Cmic)).repeat(len(d.Time)).reshape(len(d.Time),1),
                            np.nanstd((d.ATP)).repeat(len(d.Time)).reshape(len(d.Time),1)),
                       axis=1)
    
    out=np.nansum(((yhat_full-obs)/weights)**2)

    return out
