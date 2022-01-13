import numpy as np
import pandas as pd
from scipy.integrate import odeint

def calcDEB_B (model, pars, t, y0):
    #model parameters
    ##yA, Km, v, m, g, ce
    pars_model=pars[0:6]
    #conversion factors
    ##ce, ne, nX1
    conversions=pars[5:8]
    
    #solve the model
    y=odeint(model,y0,t, args=(pars_model,))
    #calculate labelled biomass (Bl) and kec factor
    Bl=0.25*conversions[0]*y[:, 2] + conversions[0]*y[:, 1]*(y[:, 2] + y[:, 4])
    kec = (0.25*conversions[2] + conversions[1]*y[:, 1])/(0.25 + y[:, 1])
                
    #Create data with predictions
    yhat = np.concatenate((y[:, 0].reshape(6,1), #Glucose
                           y[:, 3].reshape(6,1),#CO2
                           kec.reshape(6,1),
                           Bl.reshape(6,1)), axis=1)
    
    return yhat
