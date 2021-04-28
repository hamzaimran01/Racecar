import sys
sys.path.append("../THESIS")
import numpy as np
from casadi import *
import matplotlib.pyplot as plt
from common_parameters.common import *
from model.track_load import *

s   = MX.sym('s')
n   = MX.sym('n')
alpha = MX.sym('alpha')
v   = MX.sym('v')
x = vertcat(s,n,alpha,v)
#u_states
derD = MX.sym('derD')
derDelta=MX.sym('derDelta')
u = vertcat(derD,derDelta)
#parameters for ode
Fxd=(Cm1-Cm2*v)*derD-Cr2*v*v-Cr0*tanh(5*v)
sdota=(v*np.cos(alpha+C1*derDelta))/(1-kapparef_s(s)*n)


# model in race track coordinates
DT=0.02
##non augmented
xdot = vertcat(sdota,\
v*sin(alpha+C1*derDelta),\
v*C2*derDelta-kapparef_s(s)*sdota,\
Fxd/m*cos(C1*derDelta),\
)
f_rk4 = Function('f_rk4',[x,u],[xdot])
X0 = MX.sym('X0', 4, 1)
U0 = MX.sym('U0', 2, 1)
k1 = f_rk4(X0, U0)
k2 = f_rk4(X0 + DT/2 * k1, U0)
k3 = f_rk4(X0 + DT/2 * k2, U0)
k4 = f_rk4(X0 + DT * k3, U0)
X_out = X0 + DT/6*(k1 +2*k2 +2*k3 +k4)
integrator_fun=Function('F_integrator',[X0,U0],[X_out])
integrator_jac_x_cont=Function('dx_continous',[x,u],[jacobian(xdot,x)])
integrator_jac_u_cont=Function('du_continous',[x,u],[jacobian(xdot,u)])
integrator_jac_x = Function('integrator_jac_x', [X0, U0], [jacobian(X_out,X0)])

C = np.array([[1, 0, 0 ,0 ],\
[0, 1, 0 ,0 ], \
[0, 0, 1 ,0 ]])
D=np.array([[0,0],[0,0],[0,0]])

