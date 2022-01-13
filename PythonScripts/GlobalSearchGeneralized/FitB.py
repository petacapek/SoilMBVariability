import numpy as np
import pandas as pd
from scipy.integrate import odeint
from calcDEB_B import calcDEB_B
from DEBmodelIso import DEBmodelIso

def FitB (x):
    #define parameters
    ##yA, Km, v, m, g, ce, nX1, eu_i = x[17]
    pars = x[0:7]
    #read all data
    dAll = pd.read_csv('/mnt/580CBE2464C5F83D/pracovni/data_statistika/SoilMBVariability/SoilMBVariabilityData/BremerKessel1990.csv', sep=',')
    #==============================================High carbon and nitrogen
    d = dAll[(dAll.Treatment=="HCHN")]
    d = d.reset_index()
    #initial conditions
    Sl_i = d.Sinit[0]
    #eu_i = x[17]
    X1u_i = d.Cmicinit[0]/(pars[6]*pars[5]*0.25)
    
    y0 = np.array([Sl_i, 0, 0, 0, X1u_i])
    #times
    t = d.Time
    #model simulations
    yhat_full = calcDEB_B(DEBmodelIso, pars, t, y0)
    #Convert to array
    yhat = np.concatenate([yhat_full[:, 0], yhat_full[:, 1], yhat_full[:, 2], yhat_full[:, 3]])
    #observations
    obs=np.concatenate([d.S, d.CO214cumul, d.kec, d.Cmic14])
    #variables
    var = np.concatenate([['S']*len(d.Time), ['CO2']*len(d.Time), ['kec']*len(d.Time), ['MBC']*len(d.Time)])
    
    #Study
    study = ['Bremer and van Kessel (1990)']*len(d.Time)*4
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #==============================================High carbon
    d2 = dAll[(dAll.Treatment=="HC")]
    d2 = d2.reset_index()
    #initial conditions
    Sl_i2 = d2.Sinit[0]
    #eu_i = x[17]
    X1u_i2 = d2.Cmicinit[0]/(pars[6]*pars[5]*0.25)
    
    y02 = np.array([Sl_i2, 0, 0, 0, X1u_i2])
    #times
    t2 = d2.Time
    #model simulations
    yhat_full2 = calcDEB_B(DEBmodelIso, pars, t2, y02)
    #Convert to array
    yhat2 = np.concatenate([yhat_full2[:, 0], yhat_full2[:, 1], yhat_full2[:, 2], yhat_full2[:, 3]])
    #observations
    obs2=np.concatenate([d2.S, d2.CO214cumul, d2.kec, d2.Cmic14])
    #variables
    var2 = np.concatenate([['S']*len(d2.Time), ['CO2']*len(d2.Time), ['kec']*len(d2.Time), ['MBC']*len(d2.Time)])
    
    #Study
    study2 = ['Bremer and van Kessel (1990)']*len(d2.Time)*4              
                       
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #==============================================Low carbon
    d3 = dAll[(dAll.Treatment=="LC")]
    d3 = d3.reset_index()
    #initial conditions
    Sl_i3 = d3.Sinit[0]
    #eu_i = x[17]
    X1u_i3 = d3.Cmicinit[0]/(pars[6]*pars[5]*0.25)
    
    y03 = np.array([Sl_i3, 0, 0, 0, X1u_i3])
    #times
    t3 = d3.Time
    #model simulations
    yhat_full3 = calcDEB_B(DEBmodelIso, pars, t3, y03)
    #Convert to array
    yhat3 = np.concatenate([yhat_full3[:, 0], yhat_full3[:, 1], yhat_full3[:, 2], yhat_full3[:, 3]])
    #observations
    obs3=np.concatenate([d3.S, d3.CO214cumul, d3.kec, d3.Cmic14])
    #variables
    var3 = np.concatenate([['S']*len(d3.Time), ['CO2']*len(d3.Time), ['kec']*len(d3.Time), ['MBC']*len(d3.Time)])
    
    #Study
    study3 = ['Bremer and van Kessel (1990)']*len(d3.Time)*4    
    
    #Create the data frame
    data_out = {'Simulation':np.concatenate([yhat, yhat2, yhat3]), 'Observation': np.concatenate([obs, obs2, obs3]), 'Variable': np.concatenate([var, var2, var3]), 'Study': np.concatenate([study, study2, study3])}
    data_out = pd.DataFrame(data_out)
    
    return data_out
