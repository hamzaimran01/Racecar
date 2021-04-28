import sys
import numpy as np
sys.path.append("../THESIS")
from casadi import *
from spline_tasks.readDataFcn import getTrack
from plots import *

#from model.bycicle_model import bycicle_model 

#control=np.loadtxt("daniel_controls")
#measurements=np.loadtxt("daniel_states")


def casadi_integrator_fun(measurements,control,DT):
	n_steps = 3

	track="LMS_Track.txt"
	Tf = 1.0  # prediction horizon
	N = 50  # number of discretization steps
	T=0.8
	#T =0.2 #5 maximum simulation time[s]
	DT=Tf/N
	Nsim=int(N*T/Tf)
	print(T)
	
	m = 0.043
	C1=0.5
	C2=15.5
	Cm1=0.28
	Cm2=0.05
	Cr0=0.011
	Cr2=0.006
	#mu_x = 0.8 # coefficient for maximum longitudinal force Fx = mu_x * F_gravity
	#mu_y = 0.5


	track = "LMS_Track.txt"
	[s0,xref,yref,psiref,kapparef] = getTrack(track)
	kapparef_s=interpolant('kapparef_s','bspline',[s0],kapparef)


	#dynamic system

	s = MX.sym('s')
	n = MX.sym('n')
	alpha = MX.sym('alpha')
	v   = MX.sym('v')
	D = MX.sym('D')
	delta = MX.sym('delta')
	x = vertcat(s,n,alpha,v,D,delta)

	# controls

	derD = MX.sym('derD')
	derDelta = MX.sym('derDelta')
	u = vertcat(derD,derDelta)

	Fxd = (Cm1 - Cm2 * v) * D - Cr2 * v * v - Cr0 * tanh(5 * v)
	sdota = ( v* np.cos( alpha+ C1* delta)) / (1- kapparef_s (s) * n)


	xdot = vertcat( sdota,\
	v * sin( alpha + C1 * delta),\
	v * C2 * delta - kapparef_s(s) * sdota,\
	Fxd/m* cos( C1* delta),\
	derD,\
	derDelta
	)


	#applying rk4 function scheme here
	f_rk4 = Function('f_rk4',[x,u],[xdot])

	X0 = MX.sym('X0',6,1)
	U0 = MX.sym('U0',2,1)
	DT = DT/n_steps
        
	k1 = f_rk4(X0, U0)
	k2 = f_rk4(X0 + DT/2 * k1, U0)
	k3 = f_rk4(X0 + DT/2 * k2, U0)
	k4 = f_rk4(X0 + DT * k3, U0)
	X_out = X0 + DT/6*(k1 +2*k2 +2*k3 +k4)		
	
	integrator_fun=Function('integrator_fun',[X0,U0],[X_out])

	simX = np.ndarray((Nsim,6))
	#x0=np.array([-3.510884924689705588e-02, 1.074992450288252965e-01, 1.750991154203301037e-01, 2.042442022843528715e+00, 7.150026169364859241e-01, -6.086127799402280686e-02])
	#x0=np.array([4.148474454789999877e-01, 3.200575584970000165e-02, 6.271748998659999685e+00, 4.544659226380000083e-01,4.621279507789999852e-01, -1.123296743100000022e-01])
	x0=np.array(measurements[0,0:6])
	simX[0,:] = x0
	#print(x0)
	simU=control
	
	for i in range(Nsim-1):
		if(x0[1]<0.14 and x0[1]>-0.14):
			
				
			for k in range(n_steps): # for multiple steps
			
				x0=integrator_fun(x0,simU[i,0:2])
				x0=x0.T
				simX[i+1,:]=x0
				print(x0)
			
							
		else:
			print("is it 0 here?")
			print(x0)	
			
	np.savetxt('x_simulation', simX)
	
	
	return simX,measurements,s0,xref,yref,psiref



