import numpy as np
import pandas as pd
from scipy.integrate import odeint
from calcDEB_N import calcDEB_N
from DEBmodel import DEBmodel

def FitN_C (x):
    #define parameters
    ##yA, Km, v, m, g, ce, te, tX1, MwX1, Mwe
    pars = x[[0, 1, 2, 3, 4, 5, 8, 9, 10, 11]]
    
    #reading data for treatment A
    dall = pd.read_csv('/mnt/580CBE2464C5F83D/pracovni/data_statistika/SoilMBVariability/SoilMBVariabilityData/Nannipieri1977.csv', sep=',')
    d = dall[dall.Treatment == "C"]
    d = d.reset_index()

    #initial conditions
    S_i = d.Sinit[0]
            
    e_i = 0.25*(pars[7] - d.ATPinit[0]/d.Winit[0]*pars[8])/(d.ATPinit[0]/d.Winit[0]*(pars[9] - pars[6]))
    X1_i = d.Winit[0]/(pars[5]*(0.25*pars[8] + pars[9]*e_i))
    
    y0 = np.array([S_i, e_i, X1_i, 0])
    #times
    t = d.Time

    #model simulations
    yhat_full = calcDEB_N(DEBmodel, pars, t, y0)
    #Convert to array
    yhat = np.concatenate([yhat_full[:, 0], yhat_full[:, 1], yhat_full[:, 2]])
     
    #observations
    obs=np.concatenate([d.CO2,d.W, d.ATP])

    #variables
    var = np.concatenate([['CO2']*len(d.Time), ['W']*len(d.Time), ['ATP']*len(d.Time)])
    
    #Study
    study = ['Nannipieri et al. (1977)']*len(d.Time)*3
    
    #Create the data frame
    data_out = {'Simulation':yhat, 'Observation': obs, 'Variable': var, 'Study': study}
    data_out = pd.DataFrame(data_out)

    return data_out

