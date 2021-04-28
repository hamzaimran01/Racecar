import sys
sys.path.append("../THESIS")
import numpy as np
from casadi import *
import matplotlib.pyplot as plt
from model.augmented import f_rk4 as f_rk4
from model.augmented import integrator_jac as integrator_jac
#from augmented import integrator as integrator_fun
from common import *
from plots import *


N=250#preditction horizon(nr.ofsteps)
#dT sampling time
T=20
DT=T/N

# load track parameters
""" line_file=np.loadtxt('tracks/line.txt')
s0=line_file[:0]
xref=line_file[:1]
yref=line_file[:2]
psiref=line_file[:3]
kapparef=line_file[:4]

length=len(s0)
pathlength=s0[-1]
#copy loop to beginning and end
s0=np.append(s0[s0[length-1] + s0[1:length]])
kapparef=np.append(kapparefkapparef[1:length])
s0 = np.append([-s0[length-2] + s0[length-81:length-2]]s0)
kapparef = np.append(kapparef[length-80:length-1]kapparef) """




## Race car parameters
""" m = 0.043
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
 """
""" s_plot=[]
n_plot=[]
alpha_plot=[]
velocity_plot=[]
s_plot_est=[]
n_plot_est=[]
alpha_plot_est=[]
velocity_plot_est=[]
model_error_plot_est_1=[]
model_error_plot_est_2=[]
model_error_plot_est_3=[]


model_error_plot_true_1=[]
model_error_plot_true_2=[]
model_error_plot_true_3=[]

 """



""" g_1=MX.sym('g_1')
g_2=MX.sym('g_2')
g_3=MX.sym('g_3')
g=vertcat(g_1g_2g_3)



#dealing with 4 states for now
s   = MX.sym('s')
n   = MX.sym('n')
alpha = MX.sym('alpha')
v   = MX.sym('v')
#D = MX.sym('D')
#delta=MX.sym('delta')
x = vertcat(snalphavg)


 # controls
derD = MX.sym('derD')
derDelta=MX.sym('derDelta')
u = vertcat(derDderDelta)

#parameters needed for the ode
Fxd=(Cm1-Cm2*v)*derD-Cr2*v*v-Cr0*tanh(5*v)
sdota=(v*np.cos(alpha+C1*derDelta))/(1-kapparef[0]*n)


# model in race track coordinates
xdot = vertcat(sdota+g_1\
v*sin(alpha+C1*derDelta)+g_2\
v*C2*derDelta-kapparef[0]*sdota+g_3\
Fxd/m*cos(C1*derDelta)\
0\
0\
0)






#################################################function for rk4###########################################################
#to do function for jacobian!


f_rk4=Function('f_rk4'[xu][xdot])

X0 = MX.sym('X0' 7 1)
U0 = MX.sym('U0' 2 1)
k1 = f_rk4(X0 U0)
k2 = f_rk4(X0 + DT/2 * k1 U0)
k3 = f_rk4(X0 + DT/2 * k2 U0)
k4 = f_rk4(X0 + DT * k3 U0)
X_out = X0 + DT/6*(k1 +2*k2 +2*k3 +k4)
integrator_fun = Function('integrator_fun' [X0 U0] [X_out])
integrator_jac = Function('integrator_jac' [X0 U0] [jacobian(X_outX0)])
xk_jac=jacobian(xdotx)
xk_jac_val=Function('f_jac_val'[xu][xk_jac])




"""

#initial conditionsU0 impemented directly in the controls

X_initial=[0.10.112\
1\
1\
1\
]
#U0=[K*np.cos(omega*DT)0.1]

#P_model_cov=diag([1111])#3rd diagonal term changed to 1 from 0 because otherwise its inverse twill turn into infinity
P_model_cov=diag([1111111])
# P_model_cov_true=diag([1111111])


W=0.0025
W_cov=W*diag([1111101010])
P_model_cov_est=W
V=0.25
V_cov=V*diag([111])

X_est =[0.10.10.11-550]
X_values=X_initial
X_true =X_initial


s_plot=np.append(s_plotX_true[0])
n_plot=np.append(n_plotX_true[1])
alpha_plot=np.append(alpha_plotX_true[2])
velocity_plot=np.append(velocity_plotX_true[3])
s_plot_est=np.append(s_plot_estX_est[0])
n_plot_est=np.append(n_plot_estX_est[1])
alpha_plot_est=np.append(alpha_plot_estX_est[2])
velocity_plot_est=np.append(velocity_plot_estX_est[3])
    

                  
for i in range(N):

    vk=transpose(np.random.multivariate_normal(np.zeros(3)V_cov1))
    wk=transpose(np.random.multivariate_normal(np.zeros(7)W_cov1))
    # X_est=X_est
    
    # A=xk_jac_val(X_true[0.1*K*np.cos(i*omega)0.1])
    A=integrator_jac(X_true[0.1*K*np.cos(i*omega)0.1])
    #print(A)
    k1 = f_rk4(X_true [0.1*K*np.cos(i*omega)0.1])
    k2 = f_rk4(X_true + DT/2 * k1[0.1*K*np.cos(omega*i)0.1] )
    k3 = f_rk4(X_true + DT/2 * k2 [0.1*K*np.cos(i*omega)0.1])
    k4 = f_rk4(X_true + DT * k3 [0.1*K*np.cos(i*omega)0.1])
    X_true=X_true+DT/6*(k1 +2*k2 +2*k3 +k4)
    X_true=X_true+2*wk
    # do not apply process noise to disturbance states
    X_true[4:7]=X_initial[4:7] 
    #print(X_true)
    
    s_plot=np.append(s_plotX_true[0])
    n_plot=np.append(n_plotX_true[1])
    alpha_plot=np.append(alpha_plotX_true[2])
    velocity_plot=np.append(velocity_plotX_true[3])
    model_error_plot_true_1=np.append(model_error_plot_true_1X_true[4])
    model_error_plot_true_2=np.append(model_error_plot_true_2X_true[5])
    model_error_plot_true_3=np.append(model_error_plot_true_3X_true[6])
    

    
    # P_model_cov_true=A*P_model_cov_true*transpose(A)+W_cov
    P_model_cov_true=mtimes(A mtimes(P_model_cov_est \
        transpose(A))) + W_cov
    #ok!
    #innovation steps    
    #sn and alpha are the measured states!
    
    C = np.array([[1 0 0 0 0 0 0] [0 1 0 0 0 0 0] \
        [0 0 1 0 0 0 0]])

    #C=yk_jac_val(X_true[K*np.cos(i*omega)0.1])
    #print(Vk.size())  
    P_model_cov_est=inv(inv(P_model_cov_true)+\
        mtimes(mtimes(transpose(C)inv(V_cov))C))
    
    k1_y = f_rk4(X_est [0.1*K*np.cos(i*omega)0.1])
    k2_y = f_rk4(X_est + DT/2 * k1_y[0.1*K*np.cos(omega*i)0.1] )
    k3_y = f_rk4(X_est + DT/2 * k2_y [0.1*K*np.cos(i*omega)0.1])
    k4_y = f_rk4(X_est + DT * k3_y [0.1*K*np.cos(i*omega)0.1])
    X_est=X_est+DT/6*(k1_y +2*k2_y +2*k3_y +k4_y)
    #print(X_est)
    
    #check please
    Y_meas=X_true[0:3:]+vk
    Ck_gk=np.dot(CX_est)
    #print(Ck_gk)
    
    #print(vk.shape)
    measurement_residual=Y_meas-Ck_gk
    #print(measurement_residual)
    
    Kalman_gain = mtimes(mtimes(P_model_cov_esttranspose(C))inv(V_cov))
    X_est=X_est+mtimes(Kalman_gain measurement_residual) 
    print(Kalman_gain)
    
    
    s_plot_est=np.append(s_plot_estX_est[0])
    n_plot_est=np.append(n_plot_estX_est[1])
    alpha_plot_est=np.append(alpha_plot_estX_est[2])
    velocity_plot_est=np.append(velocity_plot_estX_est[3])
    model_error_plot_est_1=np.append(model_error_plot_est_1X_est[4])
    model_error_plot_est_2=np.append(model_error_plot_est_2X_est[5])
    model_error_plot_est_3=np.append(model_error_plot_est_3X_est[6])
    
 


power=10**(-3)
ax1 = plt.subplot2grid(shape=(26) loc=(00) colspan=2)
plt.ylabel('s_plot')
plt.plot(range(N+1)s_plot'r')
plt.plot(range(N+1)s_plot_est'g')

ax2 = plt.subplot2grid((26) (02) colspan=2)
plt.ylabel('n_plot')
plt.plot(range(N+1)n_plot'r')
plt.plot(range(N+1)n_plot_est'g')

ax3 = plt.subplot2grid((26) (04) colspan=2)
plt.ylabel('alpha_plot')
plt.plot(range(N+1)alpha_plot'r')
plt.plot(range(N+1)alpha_plot_est'g')

ax4 = plt.subplot2grid((26) (11) colspan=2)
plt.ylabel('velocity_plot')
plt.plot(range(N+1)velocity_plot'r')
plt.plot(range(N+1)velocity_plot_est'g')

plt.figure(200)



ax1 = plt.subplot2grid(shape=(26) loc=(00) colspan=2)
# ax1.set(xlim=(0N) ylim=(-0.5 0.5))
plt.ylabel('model_error_plot_1')

plt.plot(range(N)model_error_plot_true_1'r')
plt.plot(range(N)model_error_plot_est_1'g')

ax2 = plt.subplot2grid((26) (02) colspan=2)
# ax2.set(xlim=(0N) ylim=(-0.5 0.5))
plt.ylabel('model_error_plot_2')

plt.plot(range(N)model_error_plot_true_2'r')
plt.plot(range(N)model_error_plot_est_2'g')


ax3 = plt.subplot2grid((26) (04) colspan=2)
# ax2.set(xlim=(0N) ylim=(-0.5 0.5))
plt.ylabel('model_error_plot_3')
plt.ylabel('10^-1')
plt.plot(range(N)model_error_plot_true_3'r')
plt.plot(range(N)model_error_plot_est_3'g')
plt.show()
