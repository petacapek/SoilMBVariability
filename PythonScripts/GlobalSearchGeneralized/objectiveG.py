import numpy as np
import pandas as pd

def objectiveG (x):
	#Read data
	d = pd.read_csv('/mnt/580CBE2464C5F83D/pracovni/data_statistika/SoilMBVariability/SoilMBVariabilityData/Glanville2016.csv', sep=',')
	#Extract observations
	obs = d.kec_original
	#Calculate predictions 
	f = ((d.Ctot[0]/12.01*10000)*x[14]+x[13])/(((d.Ctot[0]/12.01*10000)*x[14]+x[13])+x[1])
	preds = (0.25*x[6] + f*np.exp(-x[2]*d.Time/24))/(0.25 + f*np.exp(-x[2]*d.Time/24))
	#Calculate weights
	weights = np.nanstd(d.kec_original).repeat(len(d.Time))
	#Calculate and return weighted residual sum of squares
	out = np.nansum(((preds - obs)/weights)**2)
	return out
