import sys
sys.path.append("../THESIS")
import numpy as np
from casadi import *

from scipy.optimize import root
#####only for kapparef[0]=1.42857 implemented
#operating points for a non linear system
#for mpc

T=4
N=60

""" 
T=3
N=30
"""

#for lqg
""" 
T=10
N=590
 """

 #for lqr and lqg
""" 
T=8
N=500
 """ 

m =0.043
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
#MPC#########
obj=0
mpciter=0
DT=T/N 
n_states=4
n_controls=2
""" Q=np.diag ([0.1, 100, 1, 2])
R=np.diag([0.05 ,0.1]) """
xbar=np.array([0,\
    0,\
    0,\
    1])

ubar = np.array([0.08138923, 0.091973])

""" 
Q=np.diag ([0.1, 10, 1, 2 ])
R=np.diag([0.5 ,0.1]) """

##non augmented
#N_main=128
N_main= 499












