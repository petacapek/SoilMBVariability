import numpy as np
import pandas as pd
from scipy.integrate import odeint

def Tsai1997ODESolvP (model, pars, t, y0):
    #model parameters
    ##Im, Km, yA, d, m
    pars_model=pars[0:5]
    #conversion factors
    ##kec, kATP
    conversions=pars[5:7]
    
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
