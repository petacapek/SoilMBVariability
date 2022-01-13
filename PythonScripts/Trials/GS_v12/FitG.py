import numpy as np
import pandas as pd

def FitG (x):
    #Read data
    d = pd.read_csv('../Data/Glanville2016.csv', sep=',')
    #Extract observations
    obs = d.kec_original
    #Calculate predictions 
    preds = (0.25*x[7] + x[14]*np.exp(-x[2]*d.Time/24) + x[14]*np.exp(-x[3]*d.Time/24))/(0.25 + x[14]*np.exp(-x[2]*d.Time/24) + x[14]*np.exp(-x[3]*d.Time/24))
    #variables
    var = ['kec']*len(d.Time)
    #Study
    study = ['Glanville et al. (2016)']*len(d.Time)
    #Create the data frame
    data_out = {'Simulation':preds, 'Observation': obs, 'Variable': var, 'Study': study}
    data_out = pd.DataFrame(data_out)
    return data_out
