import numpy as np
import pandas as pd
from scipy.integrate import odeint
from calcDEB_Bl import calcDEB_Bl
from DEBmodel import DEBmodel

def FitBl (x):
    #define parameters
    ##yA, Km, v, m, g, ce, iX1, e_i = x[14]
    pars = x[[0, 1, 2, 3, 4, 5, 7]]
    
    #reading all data
    dall = pd.read_csv('/mnt/580CBE2464C5F83D/pracovni/data_statistika/SoilMBVariability/SoilMBVariabilityData/Blagodatskaya2014.csv', sep=',')
    
    #====================================Rhizosphere soil
    d = dall[dall.Treatment == "Rhizosphere"]
    d = d.reset_index()

    #initial conditions
    S_i = d.Sinit[0]
    e_i = x[14]
    X1_i = d.DNAinit[0]/(pars[6]*0.25*pars[5])
    
    y0 = np.array([S_i, e_i, X1_i, 0])

    #times
    t = d.Time

    #model simulations
    yhat_full = calcDEB_Bl(DEBmodel, pars, t, y0)
    
    #Convert to array
    yhat = np.concatenate([yhat_full[:, 0], yhat_full[:, 1]])

    #observations
    obs=np.concatenate([d.CO2, d.DNA])

    #variables
    var = np.concatenate([['CO2']*len(d.Time), ['DNA']*len(d.Time)])
    
    #Study
    study = ['Blagodatskaya et al. (2014)']*len(d.Time)*2
                       
    #====================================Bulk soil
    d2 = dall[dall.Treatment != "Rhizosphere"]
    d2 = d2.reset_index()

    #initial conditions
    S_i2 = d2.Sinit[0]
    X1_i2 = d2.DNAinit[0]/(pars[6]*0.25*pars[5])
    
    y02 = np.array([S_i2, e_i, X1_i2, 0])

    #times
    t2 = d2.Time

    #model simulations
    yhat_full2 = calcDEB_Bl(DEBmodel, pars, t2, y02)

    #Convert to array
    yhat2 = np.concatenate([yhat_full2[:, 0], yhat_full2[:, 1]])

    #observations
    obs2=np.concatenate([d2.CO2, d2.DNA])

    #variables
    var2 = np.concatenate([['CO2']*len(d2.Time), ['DNA']*len(d2.Time)])
    
    #Study
    study2 = ['Blagodatskaya et al. (2014)']*len(d2.Time)*2
    
    #Create the data frame
    data_out = {'Simulation':np.concatenate([yhat, yhat2]), 'Observation': np.concatenate([obs, obs2]), 'Variable': np.concatenate([var, var2]), 'Study': np.concatenate([study, study2])}
    data_out = pd.DataFrame(data_out)

    return data_out
