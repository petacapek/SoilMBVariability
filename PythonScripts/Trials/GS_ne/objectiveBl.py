import numpy as np
import pandas as pd
from scipy.integrate import odeint
from calcDEB_Bl import calcDEB_Bl
from DEBmodel import DEBmodel

def objectiveBl (x):
    #define parameters
    ##yA, Km, v, m, g, ce, iX1, e_i = x[15]
    pars = x[[0, 1, 2, 3, 4, 5, 8]]
    
    #reading all data
    dall = pd.read_csv('/mnt/580CBE2464C5F83D/pracovni/data_statistika/SoilMBVariability/SoilMBVariabilityData/Blagodatskaya2014.csv', sep=',')
    
    #====================================Rhizosphere soil
    d = dall[dall.Treatment == "Rhizosphere"]
    d = d.reset_index()

    #initial conditions
    S_i = d.Sinit[0]
    e_i = ((d.Ctot[0]/12.01*10000)*x[15]+x[14])/(((d.Ctot[0]/12.01*10000)*x[15]+x[14])+x[1])
    X1_i = d.DNAinit[0]/(pars[6]*0.25*pars[5])
    
    y0 = np.array([S_i, e_i, X1_i, 0])

    #times
    t = d.Time

    #model simulations
    yhat_full = calcDEB_Bl(DEBmodel, pars, t, y0)

    #observations
    obs=np.concatenate((np.array([d.CO2]).reshape(len(d.Time),1),
                        np.array([d.DNA]).reshape(len(d.Time),1)),
                     axis=1)

    #weights
    weights=np.concatenate((np.nanstd(d.CO2).repeat(len(d.Time)).reshape(len(d.Time),1),
                            np.nanstd(d.DNA).repeat(len(d.Time)).reshape(len(d.Time),1)),
                       axis=1)
                       
    #====================================Bulk soil
    d2 = dall[dall.Treatment != "Rhizosphere"]
    d2 = d2.reset_index()

    #initial conditions
    S_i2 = d2.Sinit[0]
    e_i2 = ((d2.Ctot[0]/12.01*10000)*x[15]+x[14])/(((d2.Ctot[0]/12.01*10000)*x[15]+x[14])+x[1])
    X1_i2 = d2.DNAinit[0]/(pars[6]*0.25*pars[5])
    
    y02 = np.array([S_i2, e_i, X1_i2, 0])

    #times
    t2 = d2.Time

    #model simulations
    yhat_full2 = calcDEB_Bl(DEBmodel, pars, t2, y02)

    #observations
    obs2=np.concatenate((np.array([d2.CO2]).reshape(len(d2.Time),1),
                        np.array([d2.DNA]).reshape(len(d2.Time),1)),
                     axis=1)

    #weights
    weights2=np.concatenate((np.nanstd(d2.CO2).repeat(len(d2.Time)).reshape(len(d2.Time),1),
                            np.nanstd(d2.DNA).repeat(len(d2.Time)).reshape(len(d2.Time),1)),
                       axis=1)
                       
    #============Merging
    yhat = np.concatenate([yhat_full, yhat_full2])
    obs_all = np.concatenate([obs, obs2])
    weights_all = np.concatenate([weights, weights2])

    out=np.nansum(((yhat-obs_all)/weights_all)**2)

    return out
