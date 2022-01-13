#Import libraries
import numpy as np
import pandas as pd
from scipy.integrate import odeint
from scipy.optimize import dual_annealing
from scipy.optimize import differential_evolution

#Import modules
from objectiveWP import objectiveWP
from objectiveM import objectiveM
from objectiveG import objectiveG
from objectiveB import objectiveB
#from objectiveZ import objectiveZ - do not estimate conversion to PLFA (only one study and two parameters) 
from objectiveBl import objectiveBl
from objectiveS import objectiveS
from objectiveJ import objectiveJ
from objectiveT import objectiveT
from objectiveN_A import objectiveN_A
from objectiveN_B import objectiveN_B
from objectiveN_C import objectiveN_C

#Define objective function
def objectiveGlobal (x):
    out = sum([objectiveWP(x), objectiveG(x), objectiveM(x), objectiveS(x), objectiveB(x), objectiveBl(x), 
     objectiveT(x), objectiveJ(x), objectiveN_A(x), objectiveN_B(x), objectiveN_C(x)])
    return out
    
#Define parameters
#Individual parameters
yA = (0.01, 1)   #[0]
Km = (0.1, 4000) #[1]  
v1 = (0.01, 30)   #[2]
v2 = (0.01, 30)   #[3]           
m = (1e-15, 0.1)  #[4]         
g = (0.08, 5)     #[5]     
ce = (0.3, 7)     #[6]
nX1 = (0.1, 0.45) #[7]         
iX1 = (0.001, 0.5) #[8]          
te = (0.00001, 0.2)#[9]           
tX1 = (0.00001, 0.2)#[10]           
MwX1 = (1, 100)  #[11]        
Mwe = (1, 100)   #[12]
eW = (0, 1)       #[13]
eG = (0, 1)       #[14]   
eBl = (0, 1)      #[15]


#Parameter space
parSpace = [yA, Km, v1, v2, m, g, ce, nX1, iX1, te, tX1, MwX1, Mwe, 
            eW, eG, eBl]
            
#Estimate parameters
parGlobalDE = differential_evolution(objectiveGlobal, parSpace, maxiter = 10000000, workers = -1)
print(parGlobalDE)
np.savetxt('parGlobalDE.csv', parGlobalDE.x.reshape(1,16), delimiter=",")
#parGlobalF = forest_minimize(objectiveGlobal, parSpace, n_calls = 10000, n_jobs = -1)
#print(parGlobalF)
#parGlobalF = np.array(parGlobalF.x)
#np.savetxt('parGlobalF.csv', parGlobalF.reshape(1,18), delimiter=",")

#Compare simulations with observations
from FitWP import FitWP
from FitM import FitM
from FitG import FitG
from FitB import FitB
#from objectiveZ import objectiveZ - do not estimate conversion to PLFA (only one study and two parameters) 
from FitBl import FitBl
from FitS import FitS
from FitJ import FitJ
from FitT import FitT
from FitN_A import FitN_A
from FitN_B import FitN_B
from FitN_C import FitN_C

#DE
globalResultsDE = pd.concat([FitWP(parGlobalDE.x), FitM(parGlobalDE.x), FitG(parGlobalDE.x), FitB(parGlobalDE.x),
                            FitBl(parGlobalDE.x), FitS(parGlobalDE.x), FitJ(parGlobalDE.x), FitT(parGlobalDE.x),
                            FitN_A(parGlobalDE.x), FitN_B(parGlobalDE.x), FitN_C(parGlobalDE.x)])
globalResultsDE.to_csv('globalResultsDE.csv')
#Random forests
#globalResultsF = pd.concat([FitWP(parGlobalF), FitM(parGlobalF), FitG(parGlobalF), FitB(parGlobalF),
#                           FitBl(parGlobalF), FitS(parGlobalF), FitJ(parGlobalF), FitT(parGlobalF),
#                           FitN_A(parGlobalF), FitN_B(parGlobalF), FitN_C(parGlobalF)])
#globalResultsF.to_csv('globalResultsF.csv')

#DA
parGlobalDA = dual_annealing(objectiveGlobal, parSpace, maxiter = 10000)
print(parGlobalDA)
np.savetxt('parGlobalDA.csv', parGlobalDA.x.reshape(1,16), delimiter=",")
globalResultsDA = pd.concat([FitWP(parGlobalDA.x), FitM(parGlobalDA.x), FitG(parGlobalDA.x), FitB(parGlobalDA.x),
                            FitBl(parGlobalDA.x), FitS(parGlobalDA.x), FitJ(parGlobalDA.x), FitT(parGlobalDA.x),
                            FitN_A(parGlobalDA.x), FitN_B(parGlobalDA.x), FitN_C(parGlobalDA.x)])
globalResultsDA.to_csv('globalResultsDA.csv')
