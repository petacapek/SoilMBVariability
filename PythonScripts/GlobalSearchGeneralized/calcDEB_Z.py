import numpy as np
import pandas as pd
from scipy.integrate import odeint

def calcDEB_Z (model, pars, t, y0):
    #model parameters
    ##yA, Km, v, m, g, ce
    pars_model=pars[0:6]
    #conversion factors
    ##ce, lX1, le
    conversions=pars[5:8]

    #solve the model
    y=odeint(model,y0,t, args=(pars_model,))

    #calculate biomass (B) and PLFA
    B=(conversions[0]/4 + conversions[0]*y[:, 1])*y[:, 2]
    PLFA = (conversions[2]/4 + conversions[1]*y[:, 1]/(0.25 + y[:, 1]))*B
    
    #Create data with predictions
    yhat = np.concatenate((y[:, 0].reshape(9,1),#glucose
                           y[:, 3].reshape(9,1),#CO2
                           PLFA.reshape(9,1)), axis=1)

    return yhat
