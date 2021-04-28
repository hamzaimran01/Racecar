
import sys
sys.path.append("../THESIS")
import numpy as np
from casadi import *
import matplotlib.pyplot as plt
from common_parameters.common import *
from model.track_load import *






length=len(s0)
pathlength=s0[-1]

C = np.array([[1, 0, 0,0, 0, 0, 0 ]\
[0, 1, 0, 0, 0, 0, 0 ] \
[0, 0, 1, 0, 0, 0, 0 ]])


s   = MX.sym('s')
n   = MX.sym('n')
alpha = MX.sym('alpha')
v   = MX.sym('v')
#D = MX.sym('D')
#delta=MX.sym('delta')
g_1=MX.sym('g_1')
g_2=MX.sym('g_2')
g_3=MX.sym('g_3')
g=vertcat(g_1,g_2,g_3)
x_a = vertcat(s,n,alpha,v,g)
 # controls
derD = MX.sym('derD')
derDelta=MX.sym('derDelta')
u_a = vertcat(derD,derDelta)
 
#parameters needed for the ode
Fxd=(Cm1-Cm2*v)*derD-Cr2*v*v-Cr0*tanh(5*v)
sdota=(v*np.cos(alpha+C1*derDelta))/(1-kapparef_s(s)*n)



# model in race track coordinates
xdot_a = vertcat(sdota+g_1\
v*sin(alpha+C1*derDelta)+g_2\
v*C2*derDelta-kapparef_s(s)*sdota+g_3\
Fxd/m*cos(C1*derDelta)\
0\
0\
0)


f_rk4_a=Function('f_rk4_a'[x_au_a][xdot_a])

X0_a = MX.sym('X0' 7 1)
U0_a = MX.sym('U0' 2 1)
k1_a = f_rk4_a(X0_a U0_a)
k2_a = f_rk4_a(X0_a + DT/2 * k1_a U0_a)
k3_a = f_rk4_a(X0_a + DT/2 * k2_a U0_a)
k4_a = f_rk4_a(X0_a + DT * k3_a U0_a)
X_out_a = X0_a + DT/6*(k1_a +2*k2_a +2*k3_a +k4_a)

integrator_fun_a = Function('integrator_fun_a' [X0_a U0_a] [X_out_a])
integrator_jac_x_cont_a=Function('integrator_jac_x_cont_a'[x_au_a][jacobian(xdot_ax_a)])
integrator_jac_u_cont_a=Function('integrator_jac_u_cont_a'[x_au_a][jacobian(xdot_au_a)])
integrator_jac_x_a = Function('integrator_jac_x_a' [X0_a U0_a] [jacobian(X_out_aX0_a)])


X0_mpc=MX.sym('X0_mpc'71)
U0_mpc=MX.sym('U0_mpc'21)


d_1_p=[]
d_2_p=[]
d_3_p=[]


data_disturbance=np.loadtxt('configurations/disturbance_varying_sin')




d_1_p=data_disturbance
d_2_p=data_disturbance
d_3_p=data_disturbance


for i in range(N_main):

    """ d_iter=np.random.multivariate_normal([000]d_cov1).T
    #print(d_iter.shape)
    d_1_p=np.append(d_1_pd_iter[0])
    d_2_p=np.append(d_2_pd_iter[1])
    d_3_p=np.append(d_3_pd_iter[2])
    d_1_p=np.append(d_1_pd_iter[0])
    d_2_p=np.append(d_2_pd_iter[1])
    d_3_p=np.append(d_3_pd_iter[2])
    d_1_p=np.append(d_1_pd_iter[0])
    d_2_p=np.append(d_2_pd_iter[1])
    d_3_p=np.append(d_3_pd_iter[2])
    d_1_p=np.append(d_1_pd_iter[0])
    d_2_p=np.append(d_2_pd_iter[1])
    d_3_p=np.append(d_3_pd_iter[2])
    d_1_p=np.append(d_1_pd_iter[0])
    d_2_p=np.append(d_2_pd_iter[1])
    d_3_p=np.append(d_3_pd_iter[2])
    d_1_p=np.append(d_1_pd_iter[0])
    d_2_p=np.append(d_2_pd_iter[1])
    d_3_p=np.append(d_3_pd_iter[2])
    d_1_p=np.append(d_1_pd_iter[0])
    d_2_p=np.append(d_2_pd_iter[1])
    d_3_p=np.append(d_3_pd_iter[2])
    d_1_p=np.append(d_1_pd_iter[0])
    d_2_p=np.append(d_2_pd_iter[1])
    d_3_p=np.append(d_3_pd_iter[2])
 """
print(d_1_p)
#print(d_1_p.shape)

#print(d_3_p)
P_model_cov=diag([1 1 1 1 10 10 10])
W_cov=diag([1 1 1 1 1 1 1])



augmented_file_d=np.loadtxt('configurations/d_estimated_values')
d_1=augmented_file_d[:0]
d_2=augmented_file_d[:1]
d_3=augmented_file_d[:2]


#X_initial=[0.986912 -0.00104412 0.0779257 0.934145 0.2 0.2 0.2]
X_initial=[0.1 0.1 0.1 0.3 0.3 0.3 0.3]

#X_initial=[0.986912 -0.00104412 0.0779257 0.934145 d_1_p[0] d_2_p[0] d_3_p[0]]

""" 
X_est=X_est =np.array([0.2\
    0.1\
    0.1\
    1.0d_1[0]d_2[0]d_3[0]]) """

X_est=X_est =np.array([0.4\
    0.3\
    0.3\
    0.5-0.20-0.20-0.20])



C = np.array([[1 0 0 0 0 0 0] [0 1 0 0 0 0 0] \
        [0 0 1 0 0 0 0]])




g_mpc=vertcat(g_1,g_2,g_3)
x_mpc=vertcat(s,n,alpha,v)

xdot_mpc=vertcat(sdota+g_1\
v*sin(alpha+C1*derDelta)+g_2\
v*C2*derDelta-kapparef_s(s)*sdota+g_3\
Fxd/m*cos(C1*derDelta)\
)
f_rk4_mpc=Function('f_rk4_mpc'[x_mpc,u_a,g_mpc][xdot_mpc])


X0_mpc=MX.sym('X0_mpc'41)
U0_mpc=MX.sym('U0_mpc'21)

k1_mpc = f_rk4_mpc(X0_mpc U0_mpcg_mpc)
k2_mpc = f_rk4_mpc(X0_mpc + DT/2 * k1_mpc U0_mpcg_mpc)
k3_mpc = f_rk4_mpc(X0_mpc + DT/2 * k2_mpc U0_mpcg_mpc)
k4_mpc = f_rk4_mpc(X0_mpc + DT * k3_mpc U0_mpcg_mpc)
X_out_mpc = X0_mpc + DT/6*(k1_mpc +2*k2_mpc +2*k3_mpc +k4_mpc)

integrator_fun_mpc=Function('integrator_fun_mpc'[X0_mpcU0_mpcg_mpc][X_out_mpc])

    