import numpy as np
import pandas as pd

def objectiveTS (x):
	#Read data
	d = pd.read_csv('/mnt/580CBE2464C5F83D/pracovni/data_statistika/SoilMBVariability/SoilMBVariabilityData/Tessier1998.csv', sep=',')
	#Extract observations
	obs = d.kec
	#Calculate predictions 
	preds = (0.25*x[6] + (d.Ctot*x[14]+x[13])/(x[1] + (d.Ctot*x[14]+x[13])))/(0.25 + (d.Ctot*x[14]+x[13])/(x[1] + (d.Ctot*x[14]+x[13])))
	#Calculate weights
	weights = np.nanstd(d.kec).repeat(len(preds))
	#Calculate and return weighted residual sum of squares
	out = np.nansum(((preds - obs)/weights)**2)
	return out
