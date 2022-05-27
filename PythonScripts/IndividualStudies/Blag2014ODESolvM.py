import numpy as np
import pandas as pd
from scipy.integrate import odeint

def Blag2014ODESolvM (model, pars, t, y0):
    #model parameters
    ##Im, Km, yA, d
    pars_model=pars[0:4]
    #conversion factors
    ##ix1
    conversions=pars[4]
    
    #solve the model
    y=odeint(model,y0,t, args=(pars_model,))
    #calculate total biomass (B), kec factor, DNA, and CFC Flush
    DNA  = y[:, 1]*conversions
               
    #Create data with predictions
    yhat = np.concatenate((#y[:, 0].reshape(len(t), 1), #Glucose
                           y[:, 2].reshape(len(t), 1),#CO2
                           DNA.reshape(len(t), 1)), 
                           axis = 1)
    
    return yhat
