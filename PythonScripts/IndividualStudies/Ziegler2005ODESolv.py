import numpy as np
import pandas as pd
from scipy.integrate import odeint

def Ziegler2005ODESolv (model, pars, t, y0):
    #model parameters
    ##Im, Km, yA, Em, m, g
    pars_model=pars[0:6]
    #conversion factors
    ##le, lX1
    conversions=pars[6:8]
    
    #solve the model
    y=odeint(model,y0,t, args=(pars_model,))
    #calculate total biomass (B), kec factor, DNA, and CFC Flush
    B  = y[:, 2] + y[:, 2]*y[:, 1]
    kl = (conversions[1] + conversions[0]*y[:, 1])/(1 + y[:, 1])
    PLFA = B*kl
                
    #Create data with predictions
    yhat = np.concatenate((y[:, 0].reshape(len(t), 1), #Glucose
                           y[:, 3].reshape(len(t), 1),#CO2
                           PLFA.reshape(len(t), 1)), 
                           axis = 1)
    
    return yhat
