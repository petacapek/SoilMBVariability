import numpy as np
import pandas as pd

def FitG (x):
    #Read data
    d = pd.read_csv('/mnt/580CBE2464C5F83D/pracovni/data_statistika/SoilMBVariability/SoilMBVariabilityData/Glanville2016.csv', sep=',')
    #Extract observations
    obs = d.kec_original
    f = ((d.Ctot[0]/12.01*10000)*x[14]+x[13])/(((d.Ctot[0]/12.01*10000)*x[14]+x[13])+x[1])
    preds = (0.25*x[6] + f*np.exp(-x[2]*d.Time/24))/(0.25 + f*np.exp(-x[2]*d.Time/24))
    #variables
    var = ['kec']*len(d.Time)
    #Study
    study = ['Glanville et al. (2016)']*len(d.Time)
    #Create the data frame
    data_out = {'Simulation':preds, 'Observation': obs, 'Variable': var, 'Study': study}
    data_out = pd.DataFrame(data_out)
    return data_out
