

import time
import sys,os

#path=os.getcwd()
print(os.getcwd())

#os.path.join(path,"acados/examples")

  
# Join various path components  

import sys
sys.path.append("/home/hamza/Desktop/THESIS/")

import matplotlib.pyplot as plt
import numpy as np
from casadi import *
from model.track_load import *
from coord_transform.spattoOrig import SpattoOrig



#from final_data import control_der_final
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
control=np.loadtxt("data_racecar_experiment/final_data/daniel_controls")

#track details and settings for paramteres
Tf=1
N=50
track="LMS_Track.txt"
#Sref, _, _, _, _] = getTrack(track)

Tf = 1.0  # prediction horizon
N = 50  # number of discretization steps
T = 10.00  # maximum simulation time[s]
sref_N = 3  # reference for final reference progress
DT=Tf/N

N_sim=int(N*T/Tf)
####usinf the casadi integrator

m = 0.043
C1=0.5
C2=15.5
Cm1=0.58
Cm2=0.15
Cr0=0.029
Cr2=0.006
mu_x = 0.8 # coefficient for maximum longitudinal force Fx = mu_x * F_gravity
mu_y = 0.5


s   = MX.sym('s')
n   = MX.sym('n')
alpha = MX.sym('alpha')
v   = MX.sym('v')
D=MX.sym('D')
delta=MX.sym('delta')
x=vertcat(s,n,alpha,v,D,delta)

 # controls
derD = MX.sym('derD')
derDelta=MX.sym('derDelta')
u = vertcat(derD,derDelta)

Fxd=(Cm1-Cm2*v)*derD-Cr2*v*v-Cr0*tanh(5*v)
sdota=(v*np.cos(alpha+C1*derDelta))/(1-kapparef_s(s)*n)

xdot = vertcat(sdota,\
v*sin(alpha+C1*derDelta),\
v*C2*derDelta-kapparef_s(s)*sdota,\
Fxd/m*cos(C1*derDelta),\
D,\
delta
)



f_rk4=Function('f_rk4',[x,u],[xdot])

X0=MX.sym('X0',6,1)
U0=MX.sym('U0',2,1)

k1 = f_rk4(X0, U0)
k2 = f_rk4(X0 + DT/2 * k1, U0)
k3 = f_rk4(X0 + DT/2 * k2, U0)
k4 = f_rk4(X0 + DT * k3, U0)
X_out = X0 + DT/6*(k1 +2*k2 +2*k3 +k4)

integrator_fun=Function('integrator_fun',[X0,U0],[X_out])

x0=[-3.510884924689705588e-02, 1.074992450288252965e-01, 1.750991154203301037e-01, 2.042442022843528715e+00, 7.150026169364859241e-01, -6.086127799402280686e-02]


""" 
x0 =[0,\
    0,\
    0,\
    0]

 """
simX = np.ndarray((N_sim,6))
print(N_sim)
print(simX.shape)
simX[0,:] = x0


simU=control
#simU[:,0]=-1*simU[:,0]
#simU[:,1]=-1*simU[:,1]
for i in range(200):
    x0=integrator_fun(x0,simU[i,0:2])
    x0=x0.T
    simX[i+1,:]=x0
    print(x0)

plt.plot(xref,yref)
xsim,ysim,psiref,vref=SpattoOrig(simX[:,0],simX[:,1],simX[:,2],simX[:,3],s0,xref,yref,psiref)
plt.plot(xsim,ysim)
plt.show()