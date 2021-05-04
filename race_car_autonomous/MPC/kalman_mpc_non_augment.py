
import sys
sys.path.append("../THESIS")

from common_parameters.common import *
from plots import *
import numpy as np
from casadi import *
import matplotlib.pyplot as plt
from coord_transform.spattoOrig import *

from scipy.optimize import root
from MPC.mpc_module_ndisturbance import *

print(sys.path)

X_initial=[0.986912, -0.00104412, 0.0779257, 0.934145]

""" X_initial=[0.1,\
    0.1,\
    0.1,\
    1]
 """
X_est =np.array([0.2,\
    0.2,\
    0.2,\
    1.5])

X_true =X_initial

P_model_cov=diag([1,1,1,1])

W=0.000025
W_cov=W*diag([1,1,1,1])
P_model_cov_est=W
V=0.00025
V_cov=V*diag([1,1,1])

print(P_model_cov_est)

xbar_fun = integrator_fun(xbar,ubar)


s_tilde_bar = xbar_fun[0] - \
    (DT*(xbar[3]*np.cos(xbar[2] \
    +C1*ubar[1])/(1-kapparef[0]*xbar[1])))
print(s_tilde_bar)
if np.linalg.norm(s_tilde_bar) > 1.0e-3:
    raise Exception('computed steady-state is not a steate-state!')



delta_x_est=X_est-xbar
delta_u_est=[0,0]

#plotting initilization
s_plot=np.append(s_plot,X_true[0])
n_plot=np.append(n_plot,X_true[1])
alpha_plot=np.append(alpha_plot,X_true[2])
velocity_plot=np.append(velocity_plot,X_true[3])
s_plot_est=np.append(s_plot_est,X_est[0])
n_plot_est=np.append(n_plot_est,X_est[1])
alpha_plot_est=np.append(alpha_plot_est,X_est[2])
velocity_plot_est=np.append(velocity_plot_est,X_est[3])

def plant_dyn(X_true,u,A,P_model_cov_est):
    P_model_cov_true=mtimes(A, mtimes(P_model_cov_est, \
        transpose(A))) + W_cov
    X_true=integrator_fun(X_true,u)+wk
    #X_true=integrator_fun(X_true,u)    
    print(X_true)
    return (X_true,P_model_cov_true)


def observer_dyn(P_model_cov_true,X_est,u):
    P_model_cov_est=inv(inv(P_model_cov_true)+\
        mtimes(mtimes(transpose(C),inv(V_cov)),C))
    X_est=integrator_fun(X_est,u)
    
    Y_meas=X_true[0:3,:]+vk
    # Y_meas=X_true[0:3,:]
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




t_current=0    


for i in range(N):
    
    vk=transpose(np.random.multivariate_normal(np.zeros(3),V_cov,1))
    wk=transpose(np.random.multivariate_normal(np.zeros(4),W_cov,1))
    
    
    u,current_cost,status,x_reference,ubar,result = mpc_module_ndisturbance(i,X_est,N_main)
    print("u value is")
    print(u)
    A = integrator_jac_x(X_true,u)

    #plant x
    X_true,P_model_cov_true = plant_dyn(X_true,u,A,P_model_cov_est)
    #observer
    X_est,P_model_cov_est=observer_dyn(P_model_cov_true,X_est,u)

    #plot plant and controller
    s_plot=np.append(s_plot,X_true[0])
    n_plot=np.append(n_plot,X_true[1])
    alpha_plot=np.append(alpha_plot,X_true[2])
    velocity_plot=np.append(velocity_plot,X_true[3])
    s_plot_est=np.append(s_plot_est,X_est[0])
    n_plot_est=np.append(n_plot_est,X_est[1])
    alpha_plot_est=np.append(alpha_plot_est,X_est[2])
    velocity_plot_est=np.append(velocity_plot_est,X_est[3])




ax1 = plt.subplot2grid(shape=(2,6), loc=(0,0), colspan=2)
plt.ylabel('s_plot')
plt.plot(range(N+1),s_plot,'r')
plt.plot(range(N+1),s_plot_est,'g')

ax2 = plt.subplot2grid((2,6), (0,2), colspan=2)
plt.ylabel('n_plot')
plt.plot(range(N+1),n_plot,'r')
plt.plot(range(N+1),n_plot_est,'g')

ax3 = plt.subplot2grid((2,6), (0,4), colspan=2)
plt.ylabel('alpha_plot')
plt.plot(range(N+1),alpha_plot,'r')
plt.plot(range(N+1),alpha_plot_est,'g')

ax4 = plt.subplot2grid((2,6), (1,1), colspan=2)
plt.ylabel('velocity_plot')
plt.plot(range(N+1),velocity_plot,'r')
plt.plot(range(N+1),velocity_plot_est,'g') 



plt.show()
#print(velocity_plot_est)
x_track,y_track,psi_track,v_track=SpattoOrig(s_plot_est,n_plot_est,alpha_plot_est,velocity_plot_est,s0,xref,yref,psiref)
plt.figure()

plt.plot(x_track,y_track)
plt.plot(xref,yref)
plt.show()

