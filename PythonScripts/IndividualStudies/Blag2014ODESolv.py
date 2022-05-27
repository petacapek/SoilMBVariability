import numpy as np
import pandas as pd
from scipy.integrate import odeint

def Blag2014ODESolv (model, pars, t, y0):
    #model parameters
    ##Im, Km, yA, Em, m, g
    pars_model=pars[0:6]
    #conversion factors
    ##ix1
    conversions=pars[6]
    
    #solve the model
    y=odeint(model,y0,t, args=(pars_model,))
    #calculate total biomass (B), kec factor, DNA, and CFC Flush
    B  = y[:, 2] + y[:, 2]*y[:, 1]
    kDNA = conversions/(1 + y[:, 1])
    DNA = B*kDNA
                
    #Create data with predictions
    yhat = np.concatenate((#y[:, 0].reshape(len(t), 1), #Glucose
                           y[:, 3].reshape(len(t), 1),#CO2
                           DNA.reshape(len(t), 1)), 
                           axis = 1)
    
    return yhat
