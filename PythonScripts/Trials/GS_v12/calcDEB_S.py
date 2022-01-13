import numpy as np
import pandas as pd
from scipy.integrate import odeint

def calcDEB_S (model, pars, t, y0):
    #model parameters
    ##yA, Km, v1, v2, m, g, ce
    pars_model=pars[0:7]
    #conversion factors
    ##ce, nX1
    conversions=pars[6:8]
    
    #solve the model
    y=odeint(model,y0,t, args=(pars_model,))
    #calculate labelled biomass (Bl), kec factor, and labelled chloroform flush (Flush14C) 
    Bl=0.25*conversions[0]*y[:, 2] + conversions[0]*y[:, 1]*(y[:, 2] + y[:, 4])
    kec = (0.25*conversions[1] + y[:, 1])/(0.25 + y[:, 1])
    Flush14C = Bl*kec
        
    #Create data with predictions
    yhat = np.concatenate((y[:, 0].reshape(6,1), #Glucose
                           y[:, 3].reshape(6,1),#CO2
                           kec.reshape(6,1),
                           Flush14C.reshape(6,1)), axis=1)
    
    return yhat
