import numpy as np
import pandas as pd
from scipy.integrate import odeint

def Ziegler2005ODESolvP (model, pars, t, y0):
    #model parameters
    ##Im, Km, yA, d, m
    pars_model=pars[0:5]
    #conversion factors
    ##kl
    conversions=pars[5]
    
    #solve the model
    y=odeint(model,y0,t, args=(pars_model,))
    #calculate total biomass (B), kec factor, DNA, and CFC Flush
    PLFA = y[:, 1]*conversions
                
    #Create data with predictions
    yhat = np.concatenate((y[:, 0].reshape(len(t), 1), #Glucose
                           y[:, 2].reshape(len(t), 1),#CO2
                           PLFA.reshape(len(t), 1)), 
                           axis = 1)
    
    return yhat
