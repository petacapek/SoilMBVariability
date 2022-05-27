import numpy as np
import pandas as pd
from scipy.integrate import odeint

def HasanODESolv (model, pars, t, y0):
    #model parameters
    ##Im, Km, yA, Em, m, g
    pars_model=pars[0:6]
    #conversion factors
    ##ne, nX1, ef
    conversions=pars[6:8]
    
    #solve the model
    y=odeint(model,y0,t, args=(pars_model,))
    #calculate variables the model is calibrated to
    B  = y[:, 2] + y[:, 2]*y[:, 1]
    #Bmeasured = B/conversions[2]
    kec = (conversions[1] + conversions[0]*y[:, 1])/(1 + y[:, 1])
    Flush = B*kec
    Walls = B*(1-kec)#/conversions[2]
                
    #Create data with predictions
    yhat = np.concatenate((y[:, 3].reshape(len(t), 1),#CO2
                           B.reshape(len(t), 1), 
                           #Flush.reshape(len(t), 1),
                           Walls.reshape(len(t), 1)), 
                           axis = 1)
    
    return yhat
