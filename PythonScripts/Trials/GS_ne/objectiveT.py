import numpy as np
import pandas as pd
from scipy.integrate import odeint
from calcDEB_T import calcDEB_T
from DEBmodel import DEBmodel

def objectiveT (x):
    #define parameters
    ##yA, Km, v, m, g, ce, nX1, te, tX1
    x0 = np.array(x)
    pars = x0[[0,1,2,3,4,5,6, 8, 9]]
    
    #read data
    d = pd.read_csv('/mnt/580CBE2464C5F83D/pracovni/data_statistika/SoilMBVariability/SoilMBVariabilityData/Tsai1997.csv', sep=',')

    #initial conditions
    S_iHG = d.Sinit[0]
    S_iLG = d.Sinit[18]
        
    e_i = 0.25*((d.ATPinit[0]/d.Cmicinit[0])*pars[6] - pars[8])/(pars[7] - (d.ATPinit[0]/d.Cmicinit[0]))
    X1_i = d.Cmicinit[0]/(0.25*pars[6] + e_i)/pars[5]
    
    y0HG = np.array([S_iHG, e_i, X1_i, 0])
    y0LG = np.array([S_iLG, e_i, X1_i, 0])

    #times
    t = d.Time[0:18]

    #model simulations
    yhat_fullHG = calcDEB_T(DEBmodel, pars, t, y0HG)
    yhat_fullLG = calcDEB_T(DEBmodel, pars, t, y0LG)
    
    yhat_full = np.concatenate((yhat_fullHG, yhat_fullLG))

    #observations
    obs=np.concatenate((np.array([d.CO2cumul]).reshape(len(d.Time),1),
                        np.array([d.Cmic]).reshape(len(d.Time),1),
                        np.array([d.ATP]).reshape(len(d.Time),1)),
                     axis=1)

    #weights
    weightsHG=np.concatenate((np.nanstd(d.CO2cumul[0:18]).repeat(len(d.Time[0:18])).reshape(len(d.Time[0:18]),1),
                            np.nanstd((d.Cmic[0:18])).repeat(len(d.Time[0:18])).reshape(len(d.Time[0:18]),1),
                            np.nanstd((d.ATP[0:18])).repeat(len(d.Time[0:18])).reshape(len(d.Time[0:18]),1)),
                       axis=1)
    weightsLG=np.concatenate((np.nanstd(d.CO2cumul[18:len(d.Time)]).repeat(18).reshape(18,1),
                            np.nanstd((d.Cmic[18:len(d.Time)])).repeat(18).reshape(18,1),
                            np.nanstd((d.ATP[18:len(d.Time)])).repeat(18).reshape(18,1)),
                       axis=1)
    
    weights = np.concatenate((weightsHG, weightsLG))

    out=np.nansum(((yhat_full-obs)/weights)**2)

    return out
