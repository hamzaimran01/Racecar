import numpy as np
from plots import plotTrackProj
from spline_tasks.readDataFcn import getTrack
from casadi import *
import matplotlib.pyplot as plt
from spline_tasks.model.tracks.readDataFcn import getTrack
import matplotlib as mpl

#X_measurement = np.loadtxt("final_data/1.4_data/states_final_alligned_1.4")
#X_measurement=np.loadtxt("final_data/0.8_training_data/states_final_alligned_0.8_training")

#X_measurement=np.loadtxt("final_data/0.8_training_data/states_final_alligned_0.8_training")
X_measurement = np.loadtxt("final_data/0.8_test_data/states_final_alligned_test_0.8")
#con = np.loadtxt("final_data/0.8_training_data/estimated_data/control_final_alligned_0.8_training.txt")
X_measurement=X_measurement[:]
""" 

track = "LMS_Track.txt"

[s0,xref,yref,psiref,kapparef] = getTrack(track)
kapparef_s=interpolant('kapparef_s','bspline',[s0],kapparef)



#plots(X_measurement,0,s0,xref,yref,psiref,keyword="estimated")
plt.show()
X_meas=[]

x_learned_0=0    
x_learned_1=0
x_learned_2=0
x_learned_3=0
x_learned_4=0
x_learned_5=0
#plotTrackProj(X_measurement , s0, filename='LMS_Track.txt', T_opt='None') 
plotTrackProj(1000 ,x_learned_0,x_learned_1,x_learned_2,x_learned_3,x_learned_4,x_learned_5, X_measurement, s0, track, T_opt='None')
plt.show()
 """
print(mpl.get_backend())