
import sys
sys.path.insert(0,"../THESIS")
import numpy as np
from casadi import *
import matplotlib.pyplot as plt
#from spattoOrig import *
from model.non_augmented import *
from model.non_augmented import f_rk4,integrator_fun,integrator_jac_x,integrator_jac_x_cont,integrator_jac_u_cont
from plots import *

from scipy.optimize import root
import control as cn
from coord_transform.spattoOrig import *

#N=100
Q=np.diag ([0.1, 10, 1, 2 ])
R=np.diag([0.5, 0.1])

X_initial=[0.1,\
    0.1,\
    0.1,\
    1] 
    

 
X_est =np.array([0.2,\
    0.2,\
    0.2,\
    1.5,])

""" 
X_initial=[0.2\
    0\
    0\
    1.5]

X_est =np.array([0.1\
    0.1\
    0.1\
    1])

 """

X_true =X_initial

P_model_cov=diag([1111])

W=0.000025
W_cov=W*diag([1,1,1,1])
P_model_cov_est=W
V=0.00025
V_cov=V*diag([1,1,1])

xbar_fun = integrator_fun(xbar,ubar)
print(xbar)
#print(f_rk4(xbarubar))
print(xbar_fun)

s_tilde_bar = xbar_fun[0] - \
    (DT*(xbar[3]*np.cos(xbar[2] \
    +C1*ubar[1])/(1-kapparef[0]*xbar[1])))


if np.linalg.norm(s_tilde_bar) > 1.0e-3:
    raise Exception('computed steady-state is not a steate-state!')


A_ss=integrator_jac_x_cont(xbar,ubar)
B_ss=integrator_jac_u_cont(xbar,ubar)
C = np.array([[1, 0, 0, 0 ],\
[0, 1, 0, 0 ], \
[0, 0, 1, 0 ]])
D=np.array([[1,0,0],[0,1,0],[0,1,0]])

delta_x_est=X_est-xbar
delta_u_est=[0,0]

#plotting initilization
s_plot=np.append(s_plot,X_true[0])
n_plot=np.append(n_plot,X_true[1])
alpha_plot=np.append(alpha_plot,delta_x_est[2])
velocity_plot=np.append(velocity_plot,delta_x_est[3])
s_plot_est=np.append(s_plot_est,delta_x_est[0])
n_plot_est=np.append(n_plot_est,delta_x_est[1])
alpha_plot_est=np.append(alpha_plot_est,delta_x_est[2])
velocity_plot_est=np.append(velocity_plot_est,delta_x_est[3])

A_ss=np.array(A_ss)
B_ss=np.array(B_ss)
#sys=cn.ss(A_ss,B_ss,C,D)    
S,L,G=cn.care(A_ss,B_ss,Q,R)


def plant_dyn(X_true,u,A,P_model_cov_est):
    P_model_cov_true=mtimes(A, mtimes(P_model_cov_est, \
        transpose(A))) + W_cov
    X_true=integrator_fun(X_true,u)+wk    
    #print(X_true)
    return (X_true,P_model_cov_true)


def observer_dyn(P_model_cov_true,X_est,u):
    P_model_cov_est=inv(inv(P_model_cov_true)+\
        mtimes(mtimes(transpose(C),inv(V_cov)),C))
    X_est=integrator_fun(X_est,u)
    
    Y_meas=X_true[0:3:]+vk
    Ck_gk=np.dot(C,X_est)
    measurement_residual=Y_meas-Ck_gk
    Kalman_gain = mtimes(mtimes(P_model_cov_est,transpose(C)),inv(V_cov))
    X_est=X_est+mtimes(Kalman_gain,measurement_residual)
    return (X_est,P_model_cov_est)


def ss(xbar_tilde):
    
    xbar_tilde[1] = xbar[1]
    xbar_tilde[2] = xbar[2]
    xbar_tilde[3] = xbar[3]
    xbar_tilde[0] = xbar_tilde[0] + (i-1)*DT*(xbar[3]*np.cos(xbar[2] \
    +C1*ubar[1])/(1-kapparef_s(s)*xbar[1]))
    
    return xbar_tilde

def LQR(A_ss,B_Ss,Q,R):
    SLG=cn.care(A_ss,B_ss,Q,R)
    delta_u=-mtimes(G,X_est[0:4]-xbar_tilde)
    u1=np.add(ubar,delta_u)
    return u1  

for i in range(N):

    vk=transpose(np.random.multivariate_normal(np.zeros(3),V_cov,1))
    wk=transpose(np.random.multivariate_normal(np.zeros(4),W_cov,1))

    xbar_tilde = np.zeros((4))
    xbar_tilde=ss(xbar_tilde)
    u=LQR(A_ss,B_ss,Q,R)
    A=integrator_jac_x(X_true,u)

    #plant x
    X_true,P_model_cov_true=plant_dyn(X_true,u,A,P_model_cov_est)
    print(X_true)
    #observer
    X_estP_model_cov_est=observer_dyn(P_model_cov_true,X_est,u)
    print(X_est)
    #plot plant and controller
    s_plot=np.append(s_plot,X_true[0])
    n_plot=np.append(n_plot,X_true[1])
    alpha_plot=np.append(alpha_plot,X_true[2])
    velocity_plot=np.append(velocity_plot,X_true[3])
    s_plot_est=np.append(s_plot_est,X_est[0])
    n_plot_est=np.append(n_plot_est,X_est[1])
    alpha_plot_est=np.append(alpha_plot_est,X_est[2])
    velocity_plot_est=np.append(velocity_plot_est,X_est[3])

plt.figure()
ax1 = plt.subplot2grid(shape=(2,6), loc=(0,0), colspan=2)
plt.ylabel('s_plot')
plt.plot(range(N+1),s_plot,'r',label='true s ')
plt.plot(range(N+1),s_plot_est,'g',label='estimated s')
plt.legend()
ax2 = plt.subplot2grid((2,6), (0,2), colspan=2)
plt.ylabel('n_plot')
plt.plot(range(N+1),n_plot,'r',label='true n ')
plt.plot(range(N+1),n_plot_est,'g',label='estimated n')
plt.legend()
ax3 = plt.subplot2grid((2,6), (0,4), colspan=2)
plt.ylabel('alpha_plot')
plt.plot(range(N+1),alpha_plot,'r',label='true alpha')
plt.plot(range(N+1),alpha_plot_est,'g',label='estimated alpha')
plt.legend()
ax4 = plt.subplot2grid((2,6), (1,1), colspan=2)
plt.ylabel('velocity_plot')
plt.plot(range(N+1),velocity_plot,'r',label='true velocity')
plt.plot(range(N+1),velocity_plot_est,'g',label='estimated velocity')
plt.legend()
plt.show()


plt.figure()
x_tracky_trackpsi_trackv_track=SpattoOrig(s_plot_est,n_plot_est,alpha_plot_est,velocity_plot_est,s0,xref,yref,psiref)
plt.plot(x_tracky_tracklabel='track in x and y coordinates')
plt.plot(xrefyreflabel='path taken by the car')
plt.legend()
plt.show()
