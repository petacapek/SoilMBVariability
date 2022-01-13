import numpy as np
import pandas as pd

def FitTS (x):
	#Read data
	d = pd.read_csv('/mnt/580CBE2464C5F83D/pracovni/data_statistika/SoilMBVariability/SoilMBVariabilityData/Tessier1998.csv', sep=',')
	#Extract observations
	obs = d.kec
	#Calculate predictions 
	preds = (0.25*x[6] + (d.Ctot*x[14]+x[13])/(x[1] + (d.Ctot*x[14]+x[13])))/(0.25 + (d.Ctot*x[14]+x[13])/(x[1] + (d.Ctot*x[14]+x[13])))
	#variables
    	var = ['kec']*len(preds)
    	#Study
    	study = ['Tessier et al. (1998)']*len(preds)
    	#Create the data frame
   	 data_out = {'Simulation':preds, 'Observation': obs, 'Variable': var, 'Study': study}
    	data_out = pd.DataFrame(data_out)
    	return data_out
