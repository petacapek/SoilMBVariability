import numpy as np
import pandas as pd
from scipy.integrate import odeint

def Bremer1990ODESolv (model, pars, t, y0):
    #model parameters
    ##Im, Km, yA, Em, m, g
    pars_model=pars[0:6]
    #conversion factors
    ##ne, nX1
    conversions=pars[6:8]
    
    #solve the model
    y=odeint(model,y0,t, args=(pars_model,))
    #calculate labelled biomass (Bl), kec factor, and 14CFC Flush
    Bl  = y[:, 2] + (y[:, 2] + y[:, 5])*y[:, 1]
    kec = (conversions[1]*y[:, 2] + conversions[0]*(y[:, 2] + y[:, 5])*y[:, 1])/Bl
                    
    #Create data with predictions
    yhat = np.concatenate((y[:, 0].reshape(len(t), 1), #Glucose
                           y[:, 3].reshape(len(t), 1),#labelled CO2
                           kec.reshape(len(t), 1), 
                           Bl.reshape(len(t), 1)), 
                           axis = 1)
    
    return yhat
