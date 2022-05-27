import numpy as np
import pandas as pd
from scipy.integrate import odeint

def Tsai1997ODESolvM (model, pars, t, y0):
    #model parameters
    ##Im, Km, yA, d
    pars_model=pars[0:4]
    #conversion factors
    ##kec, kATP
    conversions=pars[4:6]
    
    #solve the model
    y=odeint(model,y0,t, args=(pars_model,))
    #calculate total biomass (B), kec factor, DNA, and CFC Flush
    Flush = y[:, 1]*conversions[0]
    ATP = y[:, 1]*conversions[1]
                
    #Create data with predictions
    yhat = np.concatenate((y[:, 2].reshape(len(t), 1),#CO2
                           Flush.reshape(len(t), 1),
                           ATP.reshape(len(t), 1)), 
                           axis = 1)
    
    return yhat
