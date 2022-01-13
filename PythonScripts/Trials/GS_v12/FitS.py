import numpy as np
import pandas as pd
from scipy.integrate import odeint
from calcDEB_S import calcDEB_S
from DEBmodelIso import DEBmodelIso

def FitS (x):
    #define parameters
    ##yA, Km, v1, v2, m, g, ce, nX1, eu_i = x[16]
    pars = x[0:8]
    
    #read data
    d = pd.read_csv('../Data/Santruckova2004.csv', sep=',')

    #initial conditions
    Sl_i = d.Sinit[0]
    #eu_i = x[16]  
    X1u_i = d.Cmicinit[0]/(pars[6]*0.25*pars[7])
        
    y0 = np.array([Sl_i, 0, 0, 0, X1u_i])
    #times
    t = d.Time
    #model simulations
    yhat_full = calcDEB_S(DEBmodelIso, pars, t, y0)
     #Convert to array
    yhat = np.concatenate([yhat_full[:, 0], yhat_full[:, 1], yhat_full[:, 2], yhat_full[:, 3]])
    #observations
    obs=np.concatenate([d.S, d.CO214cumul, d.kec, d.Cmic14])
    #variables
    var = np.concatenate([['S']*len(d.Time), ['CO2']*len(d.Time), ['kec']*len(d.Time), ['Flush']*len(d.Time)])
    
    #Study
    study = ['Santruckova et al. (2004)']*len(d.Time)*4
    
    #Create the data frame
    data_out = {'Simulation':yhat, 'Observation': obs, 'Variable': var, 'Study': study}
    data_out = pd.DataFrame(data_out)

    return data_out
