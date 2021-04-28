import sys

import numpy as np
import time,os
#from acados_template import AcadosOcp,AcadosOcpSolver,AcadosSimSolver
#from model.bycicle_model import bycicle_model
#from acados_settings import *
#from model.tracks.readDataFcn import getTrack
import matplotlib.pyplot as plt
from plots import *
from casadi_integrator import casadi_integrator_fun
from plotFcn import plotTrackProj
 



states_columns = np.loadtxt('data/final_data/states_final_alligned',dtype=float)
der_control_columns = np.loadtxt('data/final_data/control_final_alligned',dtype=float)



#states_columns=np.loadtxt('daniel_states',dtype=float)
#der_control_columns=np.loadtxt('daniel_controls',dtype=float)



#print(der_control_columns)

Tf = 1.0  # prediction horizon
N = 50  # number of discretization steps
T =3# maximum simulation time[s]
DT=Tf/N
Nsim=int(N*T/Tf)



def acados_integrator(control,measurements):
	#change here every time
	
	control = der_control_columns
	measurements = states_columns
	#print(der_control_columns)

	track="LMS_Track.txt"
	[s0, xref, yref, psiref, kapparef] = getTrack(track)
	
	

	constraint, model, acados_solver, ocp = acados_settings(Tf, N, track)
	

	# dimensions
	nx = model.x.size()[0]
	nu = model.u.size()[0]
	ny = nx + nu
	
	
       
	# initialize data structs
	simX = np.ndarray((Nsim, nx))
	simU = np.array(control) # der_control_columns, D_delta

	
	acados_integrator = AcadosSimSolver(ocp, json_file="acados_ocp.json")
	xcurrent = x0
	simX[0,:] = xcurrent

	#print(x0)
	for i in range(0,Nsim-1):

		#i_control=int(np.floor(i/2))
		#print(i_control)
		
		
               	
		
		acados_integrator.set("x", xcurrent)
		acados_integrator.set("u", simU[i,0:2])
		
		
		status = acados_integrator.solve()
		if status != 0:
			raise Exception('acados integrator returned status {}. Exiting.'.format(status))
		
		#import pdb; pdb.set_trace()
		# update state
		xcurrent = acados_integrator.get("x")		
		simX[i+1,:] = xcurrent 
		#print(xcurrent)
		
	np.savetxt('x_simulation', simX)
	return simX,measurements,s0,xref,yref,psiref





#change name here
name = "Casadi"#Acados,Casadi

#daniel starting point
#x0=np.array([-3.510884924689705588e-02, 1.074992450288252965e-01, 1.750991154203301037e-01, 2.042442022843528715e+00, 7.150026169364859241e-01, -6.086127799402280686e-02])

x0=np.array(states_columns[0,0:6])


if (name == "Casadi"):
	simX,measurements,s0,xref,yref,psiref = casadi_integrator_fun(states_columns,der_control_columns,DT)
	
elif(name == "Acados"):
	print("###########################################'starting point")
	simX,measurements,s0,xref,yref,psiref=acados_integrator(states_columns,der_control_columns)
	#print(simX)



plotTrackProj(simX,s0,'LMS_Track.txt', 'None')
plt.show()










	
	

