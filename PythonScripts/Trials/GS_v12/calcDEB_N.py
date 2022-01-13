import numpy as np
import pandas as pd
from scipy.integrate import odeint

def calcDEB_N (model, pars, t, y0):
    #model parameters
    ##yA, Km, v1, v2, m, g, ce
    pars_model=pars[0:7]
    #conversion factors
    ##ce, te, tX1, MwX1, Mwe
    conversions=pars[6:11]

    #solve the model
    y=odeint(model,y0,t, args=(pars_model,))

    #calculate biomass (B) and total ATP
    B=(conversions[0]/4 + conversions[0]*y[:, 1])*y[:, 2]
    ATP = (conversions[2]/4 + conversions[1]*y[:, 1])/(0.25 + y[:, 1])*B
    W = (conversions[3]*conversions[0]/4 + conversions[4]*conversions[0]*y[:, 1])*y[:, 2]
    
    #Create data with predictions
    yhat = np.concatenate((y[:, 3].reshape(8,1),#CO2
                           W.reshape(8,1),
                           ATP.reshape(8,1)), axis=1)

    return yhat
