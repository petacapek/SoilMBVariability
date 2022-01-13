import numpy as np
import pandas as pd
from scipy.integrate import odeint
from calcDEB_M import calcDEB_M
from DEBmodel import DEBmodel

def FitM (x):
    #define parameters
    ##yA, Km, v1, v2, m, g, ce, nX1, iX1
    pars = x[0:9]
    
    #Import data
    d = pd.read_csv('../Data/Marstorp1999.csv', sep=',')

    #initial conditions
    S_i = d.Sinit[0]
    e_i = 0.25*(d.Cmicinit[0]/d.DNAinit[0]*pars[8] - pars[7])
    X1_i = d.DNAinit[0]/(pars[8]*0.25*pars[6])
    
    y0 = np.array([S_i, e_i, X1_i, 0])

    #times
    t = d.Time

    #model simulations
    yhat_full = calcDEB_M(DEBmodel, pars, t, y0)
    #Convert to array
    yhat = np.concatenate([yhat_full[:, 0], yhat_full[:, 1], yhat_full[:, 2], yhat_full[:, 3]])

    #observations
    obs=np.concatenate([d.S, d.CO212cumul, d.Cmic12 + d.Cmic14, d.DNA])

    #variables
    var = np.concatenate([['S']*len(d.Time), ['CO2']*len(d.Time), ['Flush']*len(d.Time), ['DNA']*len(d.Time)])
    
    #Study
    study = ['Marstorp and Witter (1999)']*len(d.Time)*4
    
    #Create the data frame
    data_out = {'Simulation':yhat, 'Observation': obs, 'Variable': var, 'Study': study}
    data_out = pd.DataFrame(data_out)

    return data_out
