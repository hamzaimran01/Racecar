import sys
sys.path.append("../THESIS")
import numpy as np
from casadi import *
import matplotlib.pyplot as plt
#from spattoOrig import *
from scipy.optimize import root
import control.matlab as cn
from coord_transform.spattoOrig import *
from common_parameters.common import *
from model.augmented import *
from model.augmented import integrator_fun_a
from model.augmented import integrator_jac_x_cont_a
from model.augmented import integrator_jac_u_cont_a
from model.augmented import integrator_jac_x_a
from plots import *


X_initial=[0.4,0.4,0.4,1.5,\
0.25,\
0.25,\
0.25,\
]


X_initial=[0.4,0.4,0.4,1,\
0.4,\
0.4,\
0.4,\
]


#U0=[K*np.cos(omega*DT),0.1]

P_model_cov=diag([1,1,1,1])#3rd diagonal term changed to 1 from 0 because otherwise its inverse twill turn into infinity
#P_model_cov=diag([1,1,1,1,1,1,1])
P_model_cov_true=diag([1,1,1,1,1,1,1])

W=0.00025
W_cov=W*diag([1,1,1,1,10,10,10])
P_model_cov_est=W
V=0.000025
V_cov=V*diag([1,1,1])

X_est =[0.1,0.1,0.1,2,-1,1,0]
X_values=X_initial
X_true =X_initial

xbar=np.array([0,\
    0,\
    0,\
    1,
    0,
    0,
    0\
    ])


ubar=np.array([0.08138923, 0.09206817])

xbar_fun=integrator_fun_a(xbar,ubar)

xbar_fun[0]=xbar_fun[0]-DT*(xbar[3]*np.cos(xbar[2] \
    +C1*ubar[1])/(1-kapparef[0]*xbar[1]))-DT*xbar[4]


s_plot=np.append(s_plot,X_true[0])
n_plot=np.append(n_plot,X_true[1])
alpha_plot=np.append(alpha_plot,X_true[2])
velocity_plot=np.append(velocity_plot,X_true[3])
s_plot_est=np.append(s_plot_est,X_est[0])
n_plot_est=np.append(n_plot_est,X_est[1])
alpha_plot_est=np.append(alpha_plot_est,X_est[2])
velocity_plot_est=np.append(velocity_plot_est,X_est[3])


A_ss=integrator_jac_x_cont_a(xbar,ubar)
A_ss=A_ss[0:4,0:4]

A_ss=np.array(A_ss)
B_ss=integrator_jac_u_cont_a(xbar,ubar)
B_ss=B_ss[0:4,0:2]
B_ss=np.array(B_ss)
Q=np.diag ([0.1, 10, 1, 2 ])
R=np.diag([0.5 ,0.1])
S,L,G=cn.care(A_ss,B_ss,Q,R)


def plant_dyn(X_true,u,A,P_model_cov_est):
    X_true=integrator_fun_a(X_true,u)
    P_model_cov_true=mtimes(A, mtimes(P_model_cov_est, \
        transpose(A))) + W_cov

    """ 
    if(i<80):
        X_true[4]=X_initial[4]
        X_true[5]=X_initial[5]
        X_true[6]=X_initial[6]
   
    elif(i>80):
        X_true[4]=-X_initial[4]
        X_true[5]=-X_initial[5]
        X_true[6]=-X_initial[6] """
   
    X_true[4:7]=X_initial[4:7]    
    print(X_true)    
    return (X_true,P_model_cov_true)


def observer_dyn(P_model_cov_true,X_est,u):
    P_model_cov_est=inv(inv(P_model_cov_true)+\
        mtimes(mtimes(transpose(C),inv(V_cov)),C))
    X_est=integrator_fun_a(X_est,u)
    
    Y_meas=X_true[0:3,:]+vk
    Ck_gk=np.dot(C,X_est)
    measurement_residual=Y_meas-Ck_gk
    Kalman_gain = mtimes(mtimes(P_model_cov_est,transpose(C)),inv(V_cov))
    X_est=X_est+mtimes(Kalman_gain,measurement_residual)
    print(Kalman_gain)
    return (X_est,P_model_cov_est)    

def ss(xbar_tilde):
    
    xbar_tilde[1] = xbar[1]
    xbar_tilde[2] = xbar[2]
    xbar_tilde[3] = xbar[3]
    xbar_tilde[0] = xbar_tilde[0] + (i-1)*DT*(xbar[3]*np.cos(xbar[2] \
     +C1*ubar[1])/(1-kapparef_s(X_true[0])*xbar[1]))
    
    return xbar_tilde

def LQR(A_ss,B_Ss,Q,R):
    S,L,G=cn.care(A_ss,B_ss,Q,R)
    delta_u=-mtimes(G,X_est[0:4]-xbar_tilde)
    u1=np.add(ubar,delta_u)
    
    return u1
         
for i in range(N):

    vk=transpose(np.random.multivariate_normal(np.zeros(3),V_cov,1))
    wk=transpose(np.random.multivariate_normal(np.zeros(7),W_cov,1))
    xbar_tilde = np.zeros((4,))
    xbar_tilde=ss(xbar_tilde)

    u1=LQR(A_ss,B_ss,Q,R)
    A=integrator_jac_x_a(X_true,u1) 
    X_true,P_model_cov_true=plant_dyn(X_true,u1,A,P_model_cov_est)

    s_plot=np.append(s_plot,X_true[0])
    n_plot=np.append(n_plot,X_true[1])
    alpha_plot=np.append(alpha_plot,X_true[2])
    velocity_plot=np.append(velocity_plot,X_true[3])
    model_error_plot_true_1=np.append(model_error_plot_true_1,X_true[4])
    model_error_plot_true_2=np.append(model_error_plot_true_2,X_true[5])
    model_error_plot_true_3=np.append(model_error_plot_true_3,X_true[6])
    

    X_est,P_model_cov_est=observer_dyn(P_model_cov_true,X_est,u1)    
  

    s_plot_est=np.append(s_plot_est,X_est[0])
    n_plot_est=np.append(n_plot_est,X_est[1])
    alpha_plot_est=np.append(alpha_plot_est,X_est[2])
    velocity_plot_est=np.append(velocity_plot_est,X_est[3])
    model_error_plot_est_1=np.append(model_error_plot_est_1,X_est[4])
    model_error_plot_est_2=np.append(model_error_plot_est_2,X_est[5])
    model_error_plot_est_3=np.append(model_error_plot_est_3,X_est[6])
    
    
ax1 = plt.subplot2grid(shape=(4,4), loc=(0,0),rowspan=2 ,colspan=2)
plt.ylabel('s')
plt.plot(range(N+1),s_plot,'r',label='true')
plt.plot(range(N+1),s_plot_est,'g',label='est')
plt.tight_layout()
plt.legend()

ax2 = ax2 = plt.subplot2grid(shape=(4,4), loc=(0,2),rowspan=2 , colspan=2)
plt.ylabel('n')
plt.plot(range(N+1),n_plot,'r',label='true')
plt.plot(range(N+1),n_plot_est,'g',label='est')
plt.tight_layout()
plt.legend()

ax3 = plt.subplot2grid(shape=(4,4), loc=(2,0),rowspan=2 , colspan=2)
plt.ylabel(r'$\alpha$')
plt.plot(range(N+1),alpha_plot,'r',label='true')
plt.plot(range(N+1),alpha_plot_est,'g',label='est')
plt.tight_layout()
plt.legend()

ax4 = plt.subplot2grid(shape=(4,4), loc=(2,2),rowspan=2 , colspan=2)
plt.ylabel('v')
plt.plot(range(N+1),velocity_plot,'r',label='true')
plt.plot(range(N+1),velocity_plot_est,'g',label='est')
plt.tight_layout()
plt.legend()

plt.tight_layout()
plt.savefig('LQR/pics/EKF_constant.pdf')
plt.show()

plt.figure(200)



ax1 = plt.subplot2grid(shape=(4,4), loc=(0,0), rowspan=2,colspan=2)
# ax1.set(xlim=(0,N), ylim=(-0.5, 0.5))
plt.xlabel('N')
plt.ylabel(r'$d_s$')
plt.plot(range(N),model_error_plot_true_1,'r',label='true')
plt.plot(range(N),model_error_plot_est_1,'g',label='est')
plt.tight_layout()
plt.legend()


ax2 = plt.subplot2grid((4,4), (0,2), rowspan=2,colspan=2)
# ax2.set(xlim=(0,N), ylim=(-0.5, 0.5))
plt.xlabel('N')
plt.ylabel(r'$d_n$')
plt.plot(range(N),model_error_plot_true_2,'r',label='true')
plt.plot(range(N),model_error_plot_est_2,'g',label='est')
plt.tight_layout()
plt.legend()

ax3 = plt.subplot2grid((4,4), (2,1), rowspan=2,colspan=2)
# ax2.set(xlim=(0,N), ylim=(-0.5, 0.5))
plt.xlabel('N')
plt.ylabel(r'$d_\alpha$')
plt.plot(range(N),model_error_plot_true_3,'r',label='true')
plt.plot(range(N),model_error_plot_est_3,'g',label='est')
plt.tight_layout()
plt.legend()

plt.savefig('LQR/pics/disturbance_n_constant.pdf')
plt.show()

ax1 = plt.subplot2grid(shape=(4,4), loc=(0,0), rowspan=2,colspan=2) 
plt.xlabel('s_reference')
plt.ylabel(r'$d_s$')

plt.plot(s0[0:N],model_error_plot_est_1,label='est')
plt.plot(s0[0:N],model_error_plot_true_1,label='true')
plt.tight_layout()
plt.legend()

ax2 = plt.subplot2grid((4,4), (0,2), rowspan=2,colspan=2)
plt.xlabel('s_reference')
plt.ylabel(r'$d_n$')
plt.plot(s0[0:N],model_error_plot_est_2,label='est')
plt.plot(s0[0:N],model_error_plot_true_2,label='true')
plt.tight_layout()
plt.legend()


ax3 = plt.subplot2grid((4,4), (2,1), rowspan=2,colspan=2)
plt.xlabel('s_reference')
plt.ylabel(r'$d_\alpha$')
plt.plot(s0[0:N],model_error_plot_est_3,label='est')
plt.plot(s0[0:N],model_error_plot_true_3,label='true')
plt.tight_layout()
plt.legend()
plt.savefig('LQR/pics/error_constant_s_reference.pdf')
plt.show()


""" if(X_initial[4]==0):
    np.savetxt("configurations/augmented_0",np.transpose([model_error_plot_est_1,model_error_plot_est_2,model_error_plot_est_3]))
elif(X_initial[4]==0.1):
    np.savetxt("configurations/augmented_0.1",np.transpose([model_error_plot_est_1,model_error_plot_est_2,model_error_plot_est_3]))
elif(X_initial[4]==0.2):
    np.savetxt("configurations/augmented_0.2",np.transpose([model_error_plot_est_1,model_error_plot_est_2,model_error_plot_est_3]))
elif(X_initial[4]==0.3):
    np.savetxt("configurations/augmented_0.3",np.transpose([model_error_plot_est_1,model_error_plot_est_2,model_error_plot_est_3]))
elif(X_initial[4]==0.5):
    np.savetxt("configurations/augmented_0.5",np.transpose([model_error_plot_est_1,model_error_plot_est_2,model_error_plot_est_3]))
elif(X_initial[4]==1.0):
    np.savetxt("configurations/augmented_1.0",np.transpose([model_error_plot_est_1,model_error_plot_est_2,model_error_plot_est_3]))
 """

plt.figure()
plt.xlabel('x')
plt.ylabel('y')

px,py,psi_i,v_i =SpattoOrig(s_plot_est,n_plot_est,alpha_plot_est,velocity_plot_est,s0[0:N],xref,yref,psiref)
plt.plot(px,py,label='Estimator output of the car')
px,py,psi_i,v_i =SpattoOrig(s_plot,n_plot,alpha_plot,velocity_plot,s0[0:N],xref,yref,psiref)
plt.plot(px,py,label='Plant output of the car')

plt.plot(xref,yref,label='track')
plt.legend()
plt.tight_layout()
plt.savefig('LQR/pics/track_LQR.pdf')
plt.show()
