import numpy as np
import pandas as pd
from scipy.integrate import odeint

def Ziegler2005ODESolv (model, pars, t, y0):
    #model parameters
    ##Im, Km, yA, Em, m, g
    pars_model=pars[0:6]
    #conversion factors
    ##le, lX1
    conversions=pars[6]
    
    #solve the model
    y=odeint(model,y0,t, args=(pars_model,))
    #calculate total biomass (B), kec factor, DNA, and CFC Flush
    Bl  = y[:, 2] + (y[:, 2] + y[:, 5])*y[:, 1]
    kl = conversions*y[:, 2]/Bl
    PLFA = Bl*kl
                
    #Create data with predictions
    yhat = np.concatenate((y[:, 0].reshape(len(t), 1), #Glucose
                           y[:, 3].reshape(len(t), 1),#CO2
                           PLFA.reshape(len(t), 1)), 
                           axis = 1)
    
    return yhat
