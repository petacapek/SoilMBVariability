import numpy as np
import pandas as pd
from scipy.integrate import odeint
from calcDEB_T import calcDEB_T
from DEBmodel import DEBmodel

def FitT (x):
    #define parameters
    ##yA, Km, v, m, g, ce, nX1, te, tX1
    pars = x[[0,1,2,3,4,5,6, 8, 9]]
    
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
    
    #Convert to array
    yhat = np.concatenate([yhat_full[:, 0], yhat_full[:, 1], yhat_full[:, 2]])

    #observations
    obs=np.concatenate([d.CO2cumul, d.Cmic, d.ATP])

    #variables
    var = np.concatenate([['CO2']*len(d.Time), ['Flush']*len(d.Time), ['ATP']*len(d.Time)])
    
    #Study
    study = ['Tsai et al. (1997)']*len(d.Time)*3
    
    #Create the data frame
    data_out = {'Simulation':yhat, 'Observation': obs, 'Variable': var, 'Study': study}
    data_out = pd.DataFrame(data_out)

    return data_out
