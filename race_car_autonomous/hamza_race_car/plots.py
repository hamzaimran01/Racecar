

import time
import sys,os
sys.path.append("../hamza_race_car")
import numpy as np
import sys
from coord_transform.spattoOrig import *
from casadi import *
import matplotlib.pyplot as plt
#from tracks.track_load import track_load
from model.tracks.readDataFcn import getTrack



# TODO: make this plot function

def plots(x_simulation,x_measurements,s0,xref,yref,psiref):
	#print(s0)
	#length=len(s0)
	#pathlength=s0[-1]
	#s0 = np.append(s0, [s0[length - 1] + s0[1:length]])
	#kapparef = np.append(kapparef, kapparef[1:length])
	#s0 = np.append([-s0[length - 2] + s0[length - 81 : length - 2]], s0)
	#kapparef = np.append(kapparef[length - 80 : length - 1], kapparef)
	#xref=np.append(xref,[xref[length-1] +xref[1:length]])
	#yref=np.append(yref,[yref[length-1]+yref[1:length]])
	#psiref=np.append(psiref,[psiref[length-1]+psiref[1:length]])
	
	xtrack,ytrack,psiref,vref = SpattoOrig(x_simulation[:,0],x_simulation[:,1],x_simulation[:,2],x_simulation[:,3],s0,xref,yref,psiref)
	#xmeas,ymeas,psimeas,vmeas = SpattoOrig(x_measurements[:,0],x_measurements[:,1],x_measurements[:,2],x_measurements[:,3],sref,xref,yref,psiref)
	
	#plt.plot(xmeas,ymeas,label='measurement')
	plt.plot(xref,yref,label ='track')
	plt.plot(xtrack,ytrack,label ='simulation time optimal')
	#plt.plot(xmeas,ymeas,label ='measurement data')
	plt.legend()	
	plt.show()
	#print(s0)



