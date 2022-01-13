import numpy as np
import pandas as pd

def FitWP (x):
    #Read data
    d = pd.read_csv('/mnt/580CBE2464C5F83D/pracovni/data_statistika/SoilMBVariability/SoilMBVariabilityData/Wardle1990.csv', sep=',')
    #Extract observations
    obs = d.kec
    #Calculate predictions 
    preds = (0.25*x[6] + x[12]*np.exp(-x[2]*d.Time))/(0.25 + x[12]*np.exp(-x[2]*d.Time))
    #Study
    study = np.concatenate([['Wardle and Parkinson (1990)']*len(d.Time)])
    #variables
    var = np.concatenate([['kec']*len(d.Time)])
    #Create the data frame
    data_out = {'Simulation':preds, 'Observation': obs, 'Variable': var, 'Study': study}
    data_out = pd.DataFrame(data_out)
    return data_out
