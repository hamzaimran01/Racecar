import sys
sys.path.append("../THESIS")

import numpy as np
from casadi import *
import matplotlib.pyplot as plt
import control.matlab as cn

#from spattoOrig import *
from common_parameters.common import *
#from common_parameters.common import *
from scipy.optimize import root

#####mpc paramters#####
U=MX.sym('U',n_controls,N)#2*N
P=MX.sym('P',n_states+N*(n_states+n_controls))#initial condition and reference
#P=MX.sym('P',n_states+N*(n_states+n_controls)+3)#initial condition and reference + disturbance for every N


X=MX.sym('X',n_states,(N+1))#4*(N+1)
X_aug=MX.sym('X_aug',7,(N+1))#7*N+1
z_k=MX.sym('z_k',n_states,(N+1))



p=np.zeros(n_states+N*(n_controls+n_states))
#p=np.zeros(n_states+N*(n_controls+n_states)+3)#initial condition and reference + disturbance for every N

xx=np.zeros((n_states,N+1))
xx1=np.zeros((N+1,4,N+1))
G=[]#constraint vectordef ss(xbar_tilde):
# constraints
args={}
#continuity constraints
lbg=np.zeros(4*(N+1))
ubg=np.zeros(4*(N+1))
lbg[0:4*(N+1):1]=0
ubg[0:4*(N+1):1]=0
#initiliazation of x
lbx=np.zeros(4*(N+1))
ubx=np.zeros(4*(N+1))
lbx=np.append(lbx,np.zeros(2*(N)),axis=0)
ubx=np.append(ubx,np.zeros(2*(N)),axis=0)

lbx[0:4*(N+1):4]=-inf#s
ubx[0:4*(N+1):4]=+inf

lbx[1:4*(N+1):4]=-8#n
ubx[1:4*(N+1):4]=+8


lbx[2:4*(N+1):4]=-10#alpha
ubx[2:4*(N+1):4]=+10
lbx[3:4*(N+1):4]=-10#v
ubx[3:4*(N+1):4]=+10

#coontrol boundaries
lbx[4*(N+1):4*(N+1)+1+2*(N+1):2]=-1#D
ubx[4*(N+1):4*(N+1)+1+2*(N+1):2]=1
lbx[4*(N+1)+1:4*(N+1)+1+2*(N+1):2]=-0.4
ubx[4*(N+1)+1:4*(N+1)+1+2*(N+1):2]=+0.4
args['lbx']=lbx
args['ubx']=ubx    
args['lbg']=lbg
args['ubg']=ubg
print(lbx)

Q=np.diag ([0.1, 100, 1, 5])
R=np.diag([0.05 ,0.1]) 
#Q=np.diag ([0.1, 100, 1, 2])
#R=np.diag([0.05 ,0.1]) 