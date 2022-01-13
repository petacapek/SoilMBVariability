import numpy as np
import pandas as pd

def objectiveWP (x):
	#Read data
	d = pd.read_csv('/mnt/580CBE2464C5F83D/pracovni/data_statistika/SoilMBVariability/SoilMBVariabilityData/Wardle1990.csv', sep=',')
	#Extract observations
	obs = d.kec
	#Calculate predictions 
	preds = (0.25*x[6] + x[12]*np.exp(-x[2]*d.Time))/(0.25 + x[12]*np.exp(-x[2]*d.Time))
	#Calculate weights
	weights = np.concatenate([np.nanstd(d.kec[0:20]).repeat(20), np.nanstd(d.kec[20:40]).repeat(20)])
	#Calculate and return weighted residual sum of squares
	out = np.nansum(((preds - obs)/weights)**2)
	return out
