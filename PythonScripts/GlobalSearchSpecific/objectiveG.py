import numpy as np
import pandas as pd

def objectiveG (x):
	#Read data
	d = pd.read_csv('/mnt/580CBE2464C5F83D/pracovni/data_statistika/SoilMBVariability/SoilMBVariabilityData/Glanville2016.csv', sep=',')
	#Extract observations
	obs = d.kec_original
	#Calculate predictions 
	preds = (0.25*x[6] + x[13]*np.exp(-x[2]*d.Time/24))/(0.25 + x[13]*np.exp(-x[2]*d.Time/24))
	#Calculate weights
	weights = np.nanstd(d.kec_original).repeat(len(d.Time))
	#Calculate and return weighted residual sum of squares
	out = np.nansum(((preds - obs)/weights)**2)
	return out
