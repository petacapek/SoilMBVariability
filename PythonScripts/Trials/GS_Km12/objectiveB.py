import numpy as np
import pandas as pd
from scipy.integrate import odeint
from calcDEB_B import calcDEB_B
from DEBmodelIso import DEBmodelIso

def objectiveB (x):
    #define parameters
    ##yA, Km, v1, v2, m, g, ce, nX1, eu_i = x[17]
    pars = x[0:8]
    #read all data
    dAll = pd.read_csv('../Data/BremerKessel1990.csv', sep=',')
    #==============================================High carbon and nitrogen
    d = dAll[(dAll.Treatment=="HCHN")]
    d = d.reset_index()
    #initial conditions
    Sl_i = d.Sinit[0]
    #eu_i = x[17]
    X1u_i = d.Cmicinit[0]/(pars[7]*pars[6]*0.25)
        
    y0 = np.array([Sl_i, 0, 0, 0, X1u_i])
    #times
    t = d.Time
    #model simulations
    yhat_full = calcDEB_B(DEBmodelIso, pars, t, y0)
    #observations
    obs=np.concatenate((np.array([d.S]).reshape(len(d.Time),1),
                        np.array([d.CO214cumul]).reshape(len(d.Time),1),
                        np.array([d.kec]).reshape(len(d.Time),1),
                        np.array([d.Cmic14]).reshape(len(d.Time),1)), axis=1)
    #weights
    weights=np.concatenate((np.nanstd(d.S).repeat(len(d.Time)).reshape(len(d.Time),1),
                            np.nanstd(d.CO214cumul).repeat(len(d.Time)).reshape(len(d.Time),1),
                            np.nanstd(d.kec).repeat(len(d.Time)).reshape(len(d.Time),1),
                            np.nanstd((d.Cmic14)).repeat(len(d.Time)).reshape(len(d.Time),1)),
                       axis=1)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #==============================================High carbon
    d2 = dAll[(dAll.Treatment=="HC")]
    d2 = d2.reset_index()
    #initial conditions
    Sl_i2 = d2.Sinit[0]
    #eu_i = x[17]
    #X1u_i2 = d2.Cmicinit[0]/(pars[6]*pars[5]*0.25)
    
    y02 = np.array([Sl_i2, 0, 0, 0, X1u_i])
    #times
    t2 = d2.Time
    #model simulations
    yhat_full2 = calcDEB_B(DEBmodelIso, pars, t2, y02)
    #observations
    obs2=np.concatenate((np.array([d2.S]).reshape(len(d2.Time),1),
                        np.array([d2.CO214cumul]).reshape(len(d2.Time),1),
                        np.array([d2.kec]).reshape(len(d2.Time),1),
                        np.array([d2.Cmic14]).reshape(len(d2.Time),1)), axis=1)
    #weights
    weights2=np.concatenate((np.nanstd(d2.S).repeat(len(d2.Time)).reshape(len(d2.Time),1),
                            np.nanstd(d2.CO214cumul).repeat(len(d2.Time)).reshape(len(d2.Time),1),
                            np.nanstd(d2.kec).repeat(len(d2.Time)).reshape(len(d2.Time),1),
                            np.nanstd((d2.Cmic14)).repeat(len(d2.Time)).reshape(len(d2.Time),1)),
                       axis=1)                   
                       
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #==============================================Low carbon
    d3 = dAll[(dAll.Treatment=="LC")]
    d3 = d3.reset_index()
    #initial conditions
    Sl_i3 = d3.Sinit[0]
    #eu_i = x[17]
    #X1u_i3 = d3.Cmicinit[0]/(pars[6]*pars[5]*0.25)
    
    y03 = np.array([Sl_i3, 0, 0, 0, X1u_i])
    #times
    t3 = d3.Time
    #model simulations
    yhat_full3 = calcDEB_B(DEBmodelIso, pars, t3, y03)
    #observations
    obs3=np.concatenate((np.array([d3.S]).reshape(len(d3.Time),1),
                        np.array([d3.CO214cumul]).reshape(len(d3.Time),1),
                        np.array([d3.kec]).reshape(len(d3.Time),1),
                        np.array([d3.Cmic14]).reshape(len(d3.Time),1)), axis=1)
    #weights
    weights3=np.concatenate((np.nanstd(d3.S).repeat(len(d3.Time)).reshape(len(d3.Time),1),
                            np.nanstd(d3.CO214cumul).repeat(len(d3.Time)).reshape(len(d3.Time),1),
                            np.nanstd(d3.kec).repeat(len(d3.Time)).reshape(len(d3.Time),1),
                            np.nanstd((d3.Cmic14)).repeat(len(d3.Time)).reshape(len(d3.Time),1)),
                       axis=1)     
    #====================Merging
    yhat = np.concatenate([yhat_full, yhat_full2, yhat_full3]) 
    obs_all = np.concatenate([obs, obs2, obs3]) 
    weights_all = np.concatenate([weights, weights2, weights3])                      
    
    out=np.nansum(((yhat-obs_all)/weights_all)**2)
    return out
