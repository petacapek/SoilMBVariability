import numpy as np
import pandas as pd
from scipy.integrate import odeint

def calcDEB_Bl (model, pars, t, y0):
    #model parameters
    ##yA, Km, v, m, g, ce
    pars_model=pars[0:6]
    #conversion factors
    ##ce, iX1
    conversions=pars[5:7]

    #solve the model
    y=odeint(model,y0,t, args=(pars_model,))

    #calculate biomass (B) and total DNA
    B=(0.25*conversions[0] + conversions[0]*y[:, 1])*y[:, 2]
    DNA = 0.25*conversions[1]/(0.25 + y[:, 1])*B
    
    #Create data with predictions
    yhat = np.concatenate((#y[:, 0].reshape(len(d.Time),1),#glucose
                           y[:, 3].reshape(39,1),#CO2
                           #Flush.reshape(len(d.Time),1),
                           DNA.reshape(39,1)), axis=1)

    return yhat
