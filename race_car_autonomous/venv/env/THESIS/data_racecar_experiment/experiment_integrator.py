

import time
import sys,os

#path=os.getcwd()
print(os.getcwd())

#os.path.join(path,"acados/examples")

  
# Join various path components  

import sys
sys.path.append("/home/hamza/Desktop/THESIS/")


import numpy as np
from casadi import *
#import acados_template.acados_sim
from acados_template import AcadosSim, AcadosSimSolver

from acados.examples.acados_python.race_cars.bycicle_model import bycicle_model
from acados.examples.acados_python.race_cars.tracks.readDataFcn import getTrack
from acados.examples.acados_python.getting_started.common.export_pendulum_ode_model import export_pendulum_ode_model
import matplotlib.pyplot as plt

#states=np.loadtxt("data_racecar_experiment/states_columns")
#con=np.loadtxt("data_racecar_experiment/control_columns",dtype=float)
#####making the forward simulator for my u
#control



""" con=con[:,1:3]
con[:,0]=-con[:,0]
con[:,1]=-con[:,1]
 """
#con=con[:,1:3]
#con[:,0]=-con[:,0]
#con[:,1]=-con[:,1]

""" s_state= states[0,:]
#s_state=s_state.T
#s_state=states[1,:]
states_pred=[]

s_pred=[]
n_pred=[]
alpha_pred=[]
vel_pred=[]
#print(con[:,1])
for i in range(4,1033):
    s_pred=np.append(s_pred,s_state[0])
    n_pred=np.append(n_pred,s_state[1])
    alpha_pred=np.append(alpha_pred,s_state[2])
    vel_pred=np.append(vel_pred,s_state[3])
    
    #print(s_state)
    s_state=integrator_fun(s_state.T,[con[i,0],con[i,1].T])

    #print(con[i,:])
    
    #print(s_state)
    #print(con[i,:])


 """



""" 
xtrack,ytrack,psiref,vref=SpattoOrig(states[:,0],states[:,1],states[:,2],states[:,3],s0,xref,yref,psiref)
xpred,ypred,psirefpred,vrefpred=SpattoOrig(s_pred,n_pred,alpha_pred,vel_pred,s0,xref,yref,psiref)

#plt.plot(xtrack,ytrack,label='measurement')
plt.plot(xref,yref,label='track')
plt.plot(xpred,ypred,label='pred')
#plt.legend()
plt.show() """


#track details and settings for paramteres
Tf=1
N=50
track="LMS_Track.txt"
[Sref, _, _, _, _] = getTrack(track)

Tf = 1.0  # prediction horizon
N = 50  # number of discretization steps
T = 10.00  # maximum simulation time[s]
sref_N = 3  # reference for final reference progress

sim=AcadosSim()
#model,constraint=bycicle_model(track)
model,constraint=bycicle_model(track)
print(model)

sim.model=model

Tf = 0.1
# dimensions
nx = model.x.size()[0]
nu = model.u.size()[0]
ny = nx + nu
Nsim = int(T * N / Tf)



# set simulation time
sim.solver_options.T = Tf
# set options
sim.solver_options.num_stages = 4
sim.solver_options.num_steps = 3
sim.solver_options.newton_iter = 3 # for implicit integrator
sim.dims.nx=nx
sim.dims.nu=nu

print(model.x)




sim=AcadosSim()

acados_integrator = AcadosSimSolver(sim)






