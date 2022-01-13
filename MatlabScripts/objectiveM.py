import numpy as np
import pandas as pd
from scipy.integrate import odeint
from calcDEB_M import calcDEB_M
from DEBmodel import DEBmodel

def objectiveM (x):
    #define parameters
    ##yA, Km, v, m, g, ce, nX1, iX1
    pars = x
    
    #Import data
    d = pd.read_csv('/mnt/580CBE2464C5F83D/pracovni/data_statistika/SoilMBVariability/SoilMBVariabilityData/Marstorp1999.csv', sep=',')

    #initial conditions
    S_i = d.Sinit[0]
    e_i = 0.25*(d.Cmicinit[0]/d.DNAinit[0]*pars[7] - pars[6])
    X1_i = d.DNAinit[0]/(pars[5]*0.25*pars[7])
    
    y0 = np.array([S_i, e_i, X1_i, 0])

    #times
    t = d.Time

    #model simulations
    yhat_full = calcDEB_M(DEBmodel, pars, t, y0)

    #observations
    obs=np.concatenate((np.array([d.S]).reshape(len(d.Time),1),
                        np.array([d.CO212cumul]).reshape(len(d.Time),1),
                        np.array([d.Cmic12 + d.Cmic14]).reshape(len(d.Time),1),
                        #np.array([d.Cmic14]).reshape(len(d.Time),1),
                        np.array([d.DNA]).reshape(len(d.Time),1)),
                     axis=1)

    #weights
    weights=np.concatenate((np.nanstd(d.S).repeat(len(d.Time)).reshape(len(d.Time),1),
                            np.nanstd(d.CO212cumul).repeat(len(d.Time)).reshape(len(d.Time),1),
                            np.nanstd((d.Cmic12 + d.Cmic14)).repeat(len(d.Time)).reshape(len(d.Time),1),
                            #np.nanmean(d.Cmic14).repeat(len(d.Time)).reshape(len(d.Time),1),
                            np.nanstd(d.DNA).repeat(len(d.Time)).reshape(len(d.Time),1)),
                       axis=1)


    out=np.nansum(((yhat_full-obs)/weights)**2)

    return out
