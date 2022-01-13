import numpy as np
import pandas as pd

def objectiveG (x):
	#Read data
	d = pd.read_csv('../Data/Glanville2016.csv', sep=',')
	#Extract observations
	obs = d.kec_original
	#Calculate predictions 
	preds = (0.25*x[7] + x[14]*np.exp(-x[2]*d.Time/24) + x[14]*np.exp(-x[3]*d.Time/24))/(0.25 + x[14]*np.exp(-x[2]*d.Time/24) + x[14]*np.exp(-x[3]*d.Time/24))
	#Calculate weights
	weights = np.nanstd(d.kec_original).repeat(len(d.Time))
	#Calculate and return weighted residual sum of squares
	out = np.nansum(((preds - obs)/weights)**2)
	return out
