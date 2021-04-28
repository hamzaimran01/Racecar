import time
import sys,os
sys.path.append("../THESIS")
import numpy as np
import sys
#from coord_transform.spattoOrig import transformProj2Orig
from coord_transform.spattoOrig import SpattoOrig,transformProj2Orig
from casadi import *
import matplotlib.pyplot as plt
from spline_tasks import readDataFcn

#from spline_tasks.model.tracks.track_load import track_load
from spline_tasks.model.tracks.readDataFcn import getTrack



# TODO: make this plot function

#def plots(x_simulation,x_measurements,s0,xref,yref,psiref,keyword):

def plotTrackProj(N_MHE,x_learned_0,x_learned_1,x_learned_2,x_learned_3,x_learned_4,x_learned_5,x_learned_6, x_measured, s0, filename, T_opt='None'):
	plt.rcParams.update({'font.size': 19})

	#[s0,xref,yref,psiref,kapparef] = track_load(filename)
	
	
	s_learned_0 = x_learned_0[0:N_MHE,0]
	n_learned_0 = x_learned_0[0:N_MHE,1]
	alpha_learned_0 = x_learned_0[0:N_MHE,2]
	v_learned_0 = x_learned_0[0:N_MHE,3]
	
	s_learned_1 = x_learned_1[0:N_MHE,0]
	n_learned_1 = x_learned_1[0:N_MHE,1]
	alpha_learned_1 = x_learned_1[0:N_MHE,2]
	v_learned_1 = x_learned_1[0:N_MHE,3]
	
	s_learned_2 = x_learned_2[0:N_MHE,0]
	n_learned_2 = x_learned_2[0:N_MHE,1]
	alpha_learned_2 = x_learned_2[0:N_MHE,2]
	v_learned_2 = x_learned_2[0:N_MHE,3]
	
	s_learned_3 = x_learned_3[0:N_MHE,0]
	n_learned_3 = x_learned_3[0:N_MHE,1]
	alpha_learned_3 = x_learned_3[0:N_MHE,2]
	v_learned_3 = x_learned_3[0:N_MHE,3]
	

	s_learned_4 = x_learned_4[0:N_MHE,0]
	n_learned_4 = x_learned_4[0:N_MHE,1]
	alpha_learned_4 = x_learned_4[0:N_MHE,2]
	v_learned_4 = x_learned_4[0:N_MHE,3]
 	
	
	s_learned_5 = x_learned_5[0:N_MHE,0]
	n_learned_5 = x_learned_5[0:N_MHE,1]
	alpha_learned_5 = x_learned_5[0:N_MHE,2]
	v_learned_5 = x_learned_5[0:N_MHE,3]
 	
 	

	s_learned_6 = x_learned_6[0:N_MHE,0]
	n_learned_6 = x_learned_6[0:N_MHE,1]
	alpha_learned_6 = x_learned_6[0:N_MHE,2]
	v_learned_6 = x_learned_6[0:N_MHE,3]
 	


	distance=0.12
    # transform data
	s_measured = x_measured[0:N_MHE,0]
	n_measured = x_measured[0:N_MHE,1]
	alpha_measured = x_measured[0:N_MHE,2]
	v_measured = x_measured[0:N_MHE,3]

	
	
	[x_learned_0, y_learned_0, _, _,idxmindist,xref,yref] = transformProj2Orig(s_learned_0, n_learned_0, alpha_learned_0, v_learned_0,filename)
	[x_learned_1, y_learned_1, _, _,idxmindist,xref,yref] = transformProj2Orig(s_learned_1, n_learned_1, alpha_learned_1, v_learned_1,filename)
	[x_learned_2, y_learned_2, _, _,idxmindist,xref,yref] = transformProj2Orig(s_learned_2, n_learned_2, alpha_learned_2, v_learned_2,filename)
	[x_learned_3, y_learned_3, _, _,idxmindist,xref,yref] = transformProj2Orig(s_learned_3, n_learned_3, alpha_learned_3, v_learned_3,filename)
	[x_learned_4, y_learned_4, _, _,idxmindist,xref,yref] = transformProj2Orig(s_learned_4, n_learned_4, alpha_learned_4, v_learned_4,filename)
	[x_learned_5, y_learned_5, _, _,idxmindist,xref,yref] = transformProj2Orig(s_learned_5, n_learned_5, alpha_learned_5, v_learned_5,filename)
	
	[x_learned_6, y_learned_6, _, _,idxmindist,xref,yref] = transformProj2Orig(s_learned_6, n_learned_6, alpha_learned_6, v_learned_6,filename)
	
	[x_measured, y_measured, _, _,idxmindist,xref,yref] = transformProj2Orig(s_measured, n_measured, alpha_measured, v_measured,filename)
    # plot racetrack map

	plt.figure()
	plt.ylim(bottom=-1.75,top=0.35)
	plt.xlim(left=-1.1,right=1.6)
	plt.ylabel('y[m]')
	plt.xlabel('x[m]')

	# Plot center line
	[Sref,Xref,Yref,Psiref,_]=getTrack(filename)

	Xboundleft = Xref-distance*np.sin(Psiref)
	Yboundleft = Yref+distance*np.cos(Psiref)
	Xboundright = Xref+distance*np.sin(Psiref)
	Yboundright = Yref-distance*np.cos(Psiref)
	plt.plot(Xboundleft,Yboundleft,color='k',linewidth=4)
	plt.plot(Xboundright,Yboundright,color='k',linewidth=4)
	plt.plot(x_learned_0,y_learned_0,'k',linestyle='--',label='first',linewidth=4)
	plt.plot(x_learned_1,y_learned_1,color='g',label='second segment',linewidth=4)
	plt.plot(x_learned_2,y_learned_2,color='r',label='third segment',linewidth=4)
	plt.plot(x_learned_3,y_learned_3,color='k',label='fourth segment',linewidth=4)
	plt.plot(x_learned_4,y_learned_4,color='m',label='fifth segment',linewidth=4)
	plt.plot(x_learned_5,y_learned_5,color='r',label='sixth segment',linewidth=4)
	plt.plot(x_learned_6,y_learned_6,color='r',label='seventh segment',linewidth=4)
	plt.plot(x_measured,y_measured,color='y',alpha=0.7,label='measurements',linewidth=4)
	#plt.plot(Xref,Yref,label='center line',linewidth=0.5)
	plt.legend()

	#plt.savefig('final_data/0.8_test_data/states_recorded.pdf')
	

	figure=plt.gcf()
	figure.set_size_inches(24, 18)

	plt.savefig('final_data/1.1_test_data/states_recorded.pdf',bbox_inches='tight',dpi=100)

	#for i in range(N_MHE):
		#print(s_measured[i])
	
	#plt.plot(xref,yref)
	
	
	
    
""" 
    #Setup plot
	plt.figure()
	plt.ylim(bottom=-1.75,top=0.35)
	plt.xlim(left=-1.1,right=1.6)
	plt.ylabel('y[m]')
	plt.xlabel('x[m]')

	# Plot center line
	[Sref,Xref,Yref,Psiref,_]=getTrack(filename)
	#plt.plot(Xref,Yref,'--',color='k')

	#plt.plot([0:250], x, 'r-')

	# Draw Trackboundaries
	Xboundleft=Xref-distance*np.sin(Psiref)
	Yboundleft=Yref+distance*np.cos(Psiref)
	Xboundright=Xref+distance*np.sin(Psiref)
	Yboundright=Yref-distance*np.cos(Psiref)
	plt.plot(Xboundleft,Yboundleft,color='k',linewidth=1)
	plt.plot(Xboundright,Yboundright,color='k',linewidth=1)

	xtrack,ytrack,psiref,vref = SpattoOrig(x_simulation[:,0],x_simulation[:,1],x_simulation[:,2],x_simulation[:,3],s0,xref,yref,psiref)
	#xmeas,ymeas,psimeas,vmeas = SpattoOrig(x_measurements[:,0],x_measurements[:,1],x_measurements[:,2],x_measurements[:,3],sref,xref,yref,psiref)

	#plt.plot(xmeas,ymeas,label='measurement')
	plt.plot(xref,yref)
	plt.plot(xtrack,ytrack,label =keyword)
	#plt.plot(xmeas,ymeas,label ='measurement data')
	plt.legend() """	
	#plt.show()
    
	#print(s0)
""" from spline_tasks.model.tracks.readDataFcn import getTrack	
track="LMS_Track.txt"
[s0, xref, yref, psiref, kapparef] = getTrack(track)	
N_MHE=1620
x_learned_0=np.loadtxt('final_data/1.1_test_data/estimated_data/simX_res.txt')
x_measured=X_measurement=np.loadtxt("final_data/1.1_test_data/states_final_alligned_1.1_test")
x_learned_1=0
x_learned_2=0
x_learned_3=0
x_learned_4=0
x_learned_5=0
plotTrackProj(N_MHE,x_learned_0,x_learned_1,x_learned_2,x_learned_3,x_learned_4,x_learned_5, x_measured, s0, track, T_opt='None')
plt.show() """