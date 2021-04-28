
import sys
#sys.path.insert(0"../common_parameters")
#print(sys.path)
#from common_parameters.common import *
#from common_parameters.MPC_common import *
from plots import *
import numpy as np
from casadi import *
import matplotlib.pyplot as plt
import control
from slycot import sb02md
from slycot import sb02mt
import control
#import control.matlab as cn
from spattoOrig import *

from scipy.optimize import root
from MPC.mpc_module import *





   





# load track parameters
line_file=np.loadtxt("tracks/LMS_Track.txt")
#[s0xrefyrefpsirefkapparef]=line_file[0]
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
kapparef = np.append(kapparef[length-80:length-1]kapparef)
kapparef_s=interpolant('kapparef_s''bspline'[s0]kapparef)

""" 
#augmented matrix initiliayation
g_1=MX.sym('g_1')
g_2=MX.sym('g_2')
g_3=MX.sym('g_3')
g=vertcat(g_1g_2g_3)
 """
#dealing with 4 states for now
""" s   = MX.sym('s')
n   = MX.sym('n')
alpha = MX.sym('alpha')
v   = MX.sym('v')
x = vertcat(snalphav) """
#D = MX.sym('D')
#delta=MX.sym('delta')
#x = vertcat(snalphavg)

# controls
""" derD = MX.sym('derD')
derDelta=MX.sym('derDelta')
u = vertcat(derDderDelta)
 """
#parameters needed for the ode
""" Fxd=(Cm1-Cm2*v)*derD-Cr2*v*v-Cr0*tanh(5*v)
sdota=(v*np.cos(alpha+C1*derDelta))/(1-kapparef_s(s)*n)
 """
# model in race track coordinates
""" xdot = vertcat(sdota\
v*sin(alpha+C1*derDelta)\
v*C2*derDelta-kapparef_s(s)*sdota\
Fxd/m*cos(C1*derDelta)
)  """

#function rk4 
""" f_rk4=Function('f_rk4'[xu][xdot])
X0 = MX.sym('X0' 4 1)
U0 = MX.sym('U0' 2 1)
k1 = f_rk4(X0 U0)
k2 = f_rk4(X0 + DT/2 * k1 U0)
k3 = f_rk4(X0 + DT/2 * k2 U0)
k4 = f_rk4(X0 + DT * k3 U0)
X_out = X0 + DT/6*(k1 +2*k2 +2*k3 +k4)
integrator_fun = Function('integrator_fun' [X0 U0] [X_out])
integrator_jac_x = Function('integrator_jac_x' [X0 U0] [jacobian(X_outX0)])
integrator_jac_u=Function('integrator_jac_u'[X0U0][jacobian(X_outU0)]) """

""" integrator_jac_x_cont=Function('dx_continous'[xu][jacobian(xdotx)])
integrator_jac_u_cont=Function('du_continous'[xu][jacobian(xdotu)]) """
#steady_state_fun=Function('function_steady_state'[X0U0][xdot])

X_initial=[0.1\
    0.1\
    0.1\
    1]

X_est =np.array([0.2\
    0.2\
    0.2\
    1.5])

X_true =X_initial

P_model_cov=diag([1111])

W=0.000025
W_cov=W*diag([1111])
P_model_cov_est=W
V=0.00025
V_cov=V*diag([111])

xbar_fun = integrator_fun(xbarubar)


s_tilde_bar = xbar_fun[0] - \
    (DT*(xbar[3]*np.cos(xbar[2] \
    +C1*ubar[1])/(1-kapparef[0]*xbar[1])))
print(s_tilde_bar)
if np.linalg.norm(s_tilde_bar) > 1.0e-3:
    raise Exception('computed steady-state is not a steate-state!')

#A_ss=integrator_jac_x_cont(xbarubar)
#B_ss=integrator_jac_u_cont(xbarubar)

C = np.array([[1 0 0 0 ]\
[0 1 0 0 ] \
[0 0 1 0 ]])
D=np.array([[00][00][00]])

delta_x_est=X_est-xbar
delta_u_est=[00]

#plotting initilization
s_plot=np.append(s_plotX_true[0])
n_plot=np.append(n_plotX_true[1])
alpha_plot=np.append(alpha_plotX_true[2])
velocity_plot=np.append(velocity_plotX_true[3])
s_plot_est=np.append(s_plot_estX_est[0])
n_plot_est=np.append(n_plot_estX_est[1])
alpha_plot_est=np.append(alpha_plot_estX_est[2])
velocity_plot_est=np.append(velocity_plot_estX_est[3])

def plant_dyn(X_trueuAP_model_cov_est):
    P_model_cov_true=mtimes(A mtimes(P_model_cov_est \
        transpose(A))) + W_cov
    X_true=integrator_fun(X_true,u)+wk
    #X_true=integrator_fun(X_trueu)    
    return (X_trueP_model_cov_true)


def observer_dyn(P_model_cov_trueX_estu):
    P_model_cov_est=inv(inv(P_model_cov_true)+\
        mtimes(mtimes(transpose(C)inv(V_cov))C))
    X_est=integrator_fun(X_estu)
    
    Y_meas=X_true[0:3:]+vk
    # Y_meas=X_true[0:3:]
    Ck_gk=np.dot(CX_est)
    measurement_residual=Y_meas-Ck_gk
    Kalman_gain = mtimes(mtimes(P_model_cov_est,transpose(C))inv(V_cov))
    X_est=X_est+mtimes(Kalman_gainmeasurement_residual)
    return (X_estP_model_cov_est)


def ss(xbar_tilde):
    
    xbar_tilde[1] = xbar[1]
    xbar_tilde[2] = xbar[2]
    xbar_tilde[3] = xbar[3]
    xbar_tilde[0] = xbar_tilde[0] + (i-1)*DT*(xbar[3]*np.cos(xbar[2] \
    +C1*ubar[1])/(1-kapparef_s(s)*xbar[1]))
    return xbar_tilde




t_current=0    


for i in range(N):
    
    vk=transpose(np.random.multivariate_normal(np.zeros(3)V_cov1))
    wk=transpose(np.random.multivariate_normal(np.zeros(4)W_cov1))

    
    u=mpc_module(iX_est)
    A=integrator_jac_x(X_trueu)

    #plant x
    X_trueP_model_cov_true=plant_dyn(X_trueuAP_model_cov_est)
    #observer
    X_estP_model_cov_est=observer_dyn(P_model_cov_trueX_estu)

    #plot plant and controller
    s_plot=np.append(s_plotX_true[0])
    n_plot=np.append(n_plotX_true[1])
    alpha_plot=np.append(alpha_plotX_true[2])
    velocity_plot=np.append(velocity_plotX_true[3])
    s_plot_est=np.append(s_plot_estX_est[0])
    n_plot_est=np.append(n_plot_estX_est[1])
    alpha_plot_est=np.append(alpha_plot_estX_est[2])
    velocity_plot_est=np.append(velocity_plot_estX_est[3])




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



plt.show()
#print(velocity_plot_est)
x_tracky_trackpsi_trackv_track=SpattoOrig(s_plot_estn_plot_estalpha_plot_estvelocity_plot_ests0xrefyrefpsiref)
plt.figure()

plt.plot(x_tracky_track)
plt.plot(xrefyref)
plt.show()

