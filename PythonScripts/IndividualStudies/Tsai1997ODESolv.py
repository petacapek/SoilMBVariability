import numpy as np
import pandas as pd
from scipy.integrate import odeint

def Tsai1997ODESolv (model, pars, t, y0):
    #model parameters
    ##Im, Km, yA, Em, m, g
    pars_model=pars[0:6]
    #conversion factors
    ##ne, nX1, te, tX1
    conversions=pars[6:10]
    
    #solve the model
    y=odeint(model,y0,t, args=(pars_model,))
    #calculate total biomass (B), CFC Flush, and ATP concentration
    B  = y[:, 2] + y[:, 2]*y[:, 1]
    kec = (conversions[1] + conversions[0]*y[:, 1])/(1 + y[:, 1])
    kATP = (conversions[3] + conversions[2]*y[:, 1])/(1 + y[:, 1])
    Flush = B*kec
    ATP = B*kATP
                
    #Create data with predictions
    yhat = np.concatenate((y[:, 3].reshape(len(t), 1),#CO2
                           Flush.reshape(len(t), 1),
                           ATP.reshape(len(t), 1)), 
                           axis = 1)
    
    return yhat
