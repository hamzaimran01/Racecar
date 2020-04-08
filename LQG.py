import numpy as np
from casadi import *
import matplotlib.pyplot as plt
import control
from slycot import sb02md
from slycot import sb02mt

from scipy.optimize import root




# load track parameters
line_file=np.loadtxt('circle_07.txt')


#[s0,xref,yref,psiref,kapparef]=line_file[0]
s0=line_file[:,0]
xref=line_file[:,1]
yref=line_file[:,2]
psiref=line_file[:,3]
kapparef=line_file[:,4]

length=len(s0)
pathlength=s0[-1]
#copy loop to beginning and end
s0=np.append(s0,[s0[length-1] + s0[1:length]])

kapparef=np.append(kapparef,kapparef[1:length])
s0 = np.append([-s0[length-2] + s0[length-81:length-2]],s0)
kapparef = np.append(kapparef[length-80:length-1],kapparef)
kapparef_s=interpolant('kapparef_s','bspline',[s0],kapparef)





## Race car parameters
m = 0.043
C1=0.5
C2=15.5
Cm1=0.58
Cm2=0.15
Cr0=0.029
Cr2=0.006
mu_x = 0.8 # coefficient for maximum longitudinal force Fx = mu_x * F_gravity
mu_y = 0.5 # coefficient for maximum lateral force Fy = mu_y * F_gravity
K=0.2
omega=0.2

#Plot initilization true and est
s_plot=[]
n_plot=[]
alpha_plot=[]
velocity_plot=[]
s_plot_est=[]
n_plot_est=[]
alpha_plot_est=[]
velocity_plot_est=[]

model_error_plot_true_1=[]
model_error_plot_true_2=[]
model_error_plot_true_3=[]
model_error_plot_est_1=[]
model_error_plot_est_2=[]
model_error_plot_est_3=[]

#augmented matrix initiliayation
g_1=MX.sym('g_1')
g_2=MX.sym('g_2')
g_3=MX.sym('g_3')
g=vertcat(g_1,g_2,g_3)

#dealing with 4 states for now
s   = MX.sym('s')
n   = MX.sym('n')
alpha = MX.sym('alpha')
v   = MX.sym('v')
x = vertcat(s,n,alpha,v)
#D = MX.sym('D')
#delta=MX.sym('delta')
#x = vertcat(s,n,alpha,v,g)

# controls
derD = MX.sym('derD')
derDelta=MX.sym('derDelta')
u = vertcat(derD,derDelta)

#parameters needed for the ode
Fxd=(Cm1-Cm2*v)*derD-Cr2*v*v-Cr0*tanh(5*v)
sdota=(v*np.cos(alpha+C1*derDelta))/(1-kapparef_s(s)*n)

# model in race track coordinates
xdot = vertcat(sdota,\
v*sin(alpha+C1*derDelta),\
v*C2*derDelta-kapparef_s(s)*sdota,\
Fxd/m*cos(C1*derDelta)
) 


N=250#preditction horizon(nr.ofsteps)
#dT sampling time
T=5
DT=T/N


#function rk4 
f_rk4=Function('f_rk4',[x,u],[xdot])
X0 = MX.sym('X0', 4, 1)
U0 = MX.sym('U0', 2, 1)
k1 = f_rk4(X0, U0)
k2 = f_rk4(X0 + DT/2 * k1, U0)
k3 = f_rk4(X0 + DT/2 * k2, U0)
k4 = f_rk4(X0 + DT * k3, U0)
X_out = X0 + DT/6*(k1 +2*k2 +2*k3 +k4)
integrator_fun = Function('integrator_fun', [X0, U0], [X_out])
integrator_jac_x = Function('integrator_jac_x', [X0, U0], [jacobian(X_out,X0)])
integrator_jac_u=Function('integrator_jac_u',[X0,U0],[jacobian(X_out,U0)])
#steady_state_fun=Function('function_steady_state',[X0,U0],[xdot])




X_initial=[0.1,\
    0.1,\
    0.1,\
    1]

X_est =np.array([0.1,\
    0.1,\
    0.1,\
    0.6])
X_values=X_initial
X_true =X_initial



P_model_cov=diag([1,1,1,1])

W=0.0025
W_cov=W*diag([1,1,1,1])
P_model_cov_est=W
V=0.0025
V_cov=V*diag([1,1,1])

Q=np.diag([1,1,1,1])
R=np.diag([1,1])



s_plot=np.append(s_plot,X_true[0])
n_plot=np.append(n_plot,X_true[1])
alpha_plot=np.append(alpha_plot,X_true[2])
velocity_plot=np.append(velocity_plot,X_true[3])
s_plot_est=np.append(s_plot_est,X_est[0])
n_plot_est=np.append(n_plot_est,X_est[1])
alpha_plot_est=np.append(alpha_plot_est,X_est[2])
velocity_plot_est=np.append(velocity_plot_est,X_est[3])



#operating points for a non linear system


xbar=np.array([0,\
    0,\
    0,\
    0.5])
ubar=np.array([0.05962736, 0.09206817])
xbar_fun=integrator_fun(xbar,ubar)
#print(xbar_fun)
#print(xbar)
uref=[0.1,0.1]
delta_u_est=uref-ubar
#jacobian @ steady state A_ss,B_ss,C,D

A_ss=integrator_jac_x(xbar,ubar)
B_ss=integrator_jac_u(xbar,ubar)



C = np.array([[1, 0, 0 ,0 ],\
[0, 1, 0 ,0 ], \
[0, 0, 1 ,0 ]])
D=np.array([[0,0],[0,0],[0,0]])
#
delta_x_est=0.1*(X_est-xbar)


#Ax+Bu
linearized_sym=mtimes(A_ss,X0)+mtimes(B_ss,U0)
f_rk4_l=Function('f_rk4_linearized',[X0,U0],[linearized_sym])
#print(mtimes(A_ss,xbar)+mtimes(B_ss,ubar))

k1_l = f_rk4_l(X0, U0)
k2_l = f_rk4_l(X0 + DT/2 * k1_l, U0)
k3_l = f_rk4_l(X0 + DT/2 * k2_l, U0)
k4_l = f_rk4_l(X0 + DT * k3_l, U0)
X_out_l = X0 + DT/6*(k1_l +2*k2_l +2*k3_l +k4_l)
integrator_fun_l = Function('integrator_fun', [X0, U0], [X_out_l])



for i in range(N):

    vk=transpose(np.random.multivariate_normal(np.zeros(3),V_cov,1))
    wk=transpose(np.random.multivariate_normal(np.zeros(4),W_cov,1))
    
    
    
    A=integrator_jac_x(X_true,delta_u_est)
   
    
    
    k1 = f_rk4(X_true, delta_u_est)
    k2 = f_rk4(X_true + DT/2 * k1,delta_u_est )
    k3 = f_rk4(X_true + DT/2 * k2, delta_u_est)
    k4 = f_rk4(X_true + DT * k3, delta_u_est)
    X_true=X_true+DT/6*(k1 +2*k2 +2*k3 +k4)
    X_true=X_true+2*wk
   
    
    s_plot=np.append(s_plot,X_true[0])
    n_plot=np.append(n_plot,X_true[1])
    alpha_plot=np.append(alpha_plot,X_true[2])
    velocity_plot=np.append(velocity_plot,X_true[3])
  

    
    
    P_model_cov_true=mtimes(A, mtimes(P_model_cov_est, \
        transpose(A))) + W_cov
    
      
    P_model_cov_est=inv(inv(P_model_cov_true)+\
        mtimes(mtimes(transpose(C),inv(V_cov)),C))



    #AX+BU integrator
    delta_x_est=integrator_fun_l(delta_x_est,delta_u_est)
    
    Y_meas=X_true[0:3,:]+vk
    Ck_gk=np.dot(C,delta_x_est)
    #print(Y_meas)
    measurement_residual=Y_meas-Ck_gk
    
    
    Kalman_gain = mtimes(mtimes(P_model_cov_est,transpose(C)),inv(V_cov))
    delta_x_est=delta_x_est+mtimes(Kalman_gain, measurement_residual)  
     
    #print(delta_x_est)
    

    #creating a state space
   
    sys=control.ss(A_ss,B_ss,C,D,DT)    
    G,S,E=control.lqr(A_ss,B_ss,Q,R)
    
    delta_u_est=-np.dot(G,delta_x_est)
    print(delta_u_est) 
    #delta_u_est=mtimes(G,delta_x_est)
    
    #print(delta_x_est)

    s_plot_est=np.append(s_plot_est,delta_x_est[0])
    n_plot_est=np.append(n_plot_est,delta_x_est[1])
    alpha_plot_est=np.append(alpha_plot_est,delta_x_est[2])
    velocity_plot_est=np.append(velocity_plot_est,delta_x_est[3])
    


power=10**(-3)
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

#plt.figure(200)




plt.show()