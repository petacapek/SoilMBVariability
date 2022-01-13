import numpy as np
import pandas as pd
from scipy.integrate import odeint

def calcDEB_T (model, pars, t, y0):
    #model parameters
    ##yA, Km, v, m, g, ce
    pars_model=pars[0:6]
    #conversion factors
    ##ce, nX1, te, tX1
    conversions=pars[5:9]

    #solve the model
    y=odeint(model,y0,t, args=(pars_model,))

    #calculate biomass (B), total ATP, and flush (Flush)
    B=(conversions[0]/4 + conversions[0]*y[:, 1])*y[:, 2]
    Flush = (conversions[1]/4 + y[:, 1])/(0.25 + y[:, 1])*B
    ATP = (conversions[3]/4 + conversions[2]*y[:, 1])/(0.25 + y[:, 1])*B
    
    #Create data with predictions
    yhat = np.concatenate((y[:, 3].reshape(18,1),#CO2
                           Flush.reshape(18,1),
                           ATP.reshape(18,1)), axis=1)

    return yhat
