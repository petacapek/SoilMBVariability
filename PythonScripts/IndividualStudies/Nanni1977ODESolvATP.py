import numpy as np
import pandas as pd
from scipy.integrate import odeint

def Nanni1977ODESolvATP (model, pars, t, y0):
    #model parameters
    ##Im, Km, yA, Em, m, g
    pars_model=pars[0:6]
    #conversion factors
    ##te, tX1
    conversions=pars[6:8]
    
    #solve the model
    y=odeint(model,y0,t, args=(pars_model,))
    #calculate total biomass (B), CFC Flush, and TAP concentration
    B  = y[:, 2] + y[:, 2]*y[:, 1]
    kATP = (conversions[1] + conversions[0]*y[:, 1])/(1 + y[:, 1])
    ATP = B*kATP
                
    #Create data with predictions
    yhat = np.concatenate((y[:, 3].reshape(len(t), 1),#CO2
                           #W.reshape(len(t), 1),
                           ATP.reshape(len(t), 1)), 
                           axis = 1)
    
    return yhat
