import numpy as np
import pandas as pd
from scipy.integrate import odeint

def calcDEB_J (model, pars, t, y0):
    #model parameters
    ##yA, Km, v, m, g, ce
    pars_model=pars[0:6]
    #conversion factors
    ##ce, ne, nX1, te, tX1
    conversions=pars[5:10]

    #solve the model
    y=odeint(model,y0,t, args=(pars_model,))

    #calculate biomass (B), total ATP, and flush (Flush)
    B=(conversions[0]/4 + conversions[0]*y[:, 1])*y[:, 2]
    Flush = (conversions[2]/4 + conversions[1]*y[:, 1])/(0.25 + y[:, 1])*B
    ATP = (conversions[4]/4 + conversions[3]*y[:, 1])/(0.25 + y[:, 1])*B
    
    #Create data with predictions
    yhat = np.concatenate((y[:, 0].reshape(7,1),#CO2
                           Flush.reshape(7,1),
                           ATP.reshape(7,1)), axis=1)

    return yhat
