import numpy as np
import pandas as pd
from scipy.integrate import odeint

def calcDEB_M (model, pars, t, y0):
    #model parameters
    ##yA, Km, v1, v2, m, g, ce
    pars_model=pars[0:7]
    #conversion factors
    ##ce, nX1, iX1
    conversions=pars[6:9]

    #solve the model
    y=odeint(model,y0,t, args=(pars_model,))

    #calculate biomass (B), total DNA, and flush (Flush)
    B=(conversions[0]/4 + conversions[0]*y[:, 1])*y[:, 2]
    Flush = (conversions[1]/4 + y[:, 1])/(0.25 + y[:, 1])*B
    DNA = 0.25*conversions[2]/(0.25 + y[:, 1])*B
    
    #Create data with predictions
    yhat = np.concatenate((y[:, 0].reshape(41,1),#glucose
                           y[:, 3].reshape(41,1),#CO2
                           Flush.reshape(41,1),
                           DNA.reshape(41,1)), axis=1)

    return yhat
