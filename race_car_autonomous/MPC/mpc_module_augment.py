""" import sys
sys.path.append("../THESIS")
from casadi import *
import numpy as np
import math
import matplotlib.pyplot as plt
#from model.non_augmented import *
#from model.non_augmented import f_rk4 as f_rk4
#from model.non_augmented import integrator_fun as integrator_fun 
from common_parameters.common import *
#from model.augmented import integrator_fun_a
from model.non_augmented import *
from model.non_augmented import integrator_fun
#from model.non_augmented import xbar,ubar
from common_parameters.MPC_common import *
#from model.augmented import *
#from model.augmented import f_rk4 as f_rk4
#from model.augmented import integrator_fun as integrator_fun
#from model.augmented import *
 """

""" U=MX.sym('U',n_controls,N)#2*N
P=MX.sym('P',n_states+N*(n_states+n_controls))#initial condition and reference 
X=MX.sym('X',n_states,(N+1))#4*(N+1)
p=np.zeros(n_states+N*(n_controls+n_states))
xx=np.zeros((n_states,N+1))
xx1=np.zeros((N+1,4,N+1))
G=[]

args={}
#continuity constraints
lbg=np.zeros(4*(N+1))
ubg=np.zeros(4*(N+1))
#initiliazation of x
lbx=np.zeros(4*(N+1))
ubx=np.zeros(4*(N+1))
lbx=np.append(lbx,np.zeros(2*(N)),axis=0)
ubx=np.append(ubx,np.zeros(2*(N)),axis=0)
lbg[0:4*(N+1):1]=0
ubg[0:4*(N+1):1]=0
lbx[0:4*(N+1):4]=-inf#s
ubx[0:4*(N+1):4]=+inf
lbx[1:4*(N+1):4]=-0.12#n
ubx[1:4*(N+1):4]=+0.12
lbx[2:4*(N+1):4]=-10#alpha
ubx[2:4*(N+1):4]=+10
lbx[3:4*(N+1):4]=-10#v
ubx[3:4*(N+1):4]=+10

#coontrol boundaries
lbx[4*(N+1):4*(N+1)+1+2*(N+1):2]=-1.0#D
ubx[4*(N+1):4*(N+1)+1+2*(N+1):2]=1.0
lbx[4*(N+1)+1:4*(N+1)+1+2*(N+1):2]=-0.40
ubx[4*(N+1)+1:4*(N+1)+1+2*(N+1):2]=+0.40
args['lbx']=lbx
args['ubx']=ubx    
args['lbg']=lbg
args['ubg']=ubg

"""

import sys
sys.path.append("../THESIS")
from casadi import *
import numpy as np
import math
import matplotlib.pyplot as plt
from model.non_augmented import *
from model.non_augmented import f_rk4 as f_rk4
from model.non_augmented import integrator_fun as integrator_fun 
from common_parameters.common import *
from common_parameters.MPC_common import *
#from model.augmented import *
#from model.augmented import f_rk4 as f_rk4
#from model.augmented import integrator_fun as integrator_fun


xbar=np.array([0,\
    0,\
    0,\
    1])

ubar = np.array([0.08138923, 0.091973])
 
t=0
t0=0    
t_current=0
u_cl=ubar

t_predict=0
U0_1=np.zeros((N,2))

def shift(DT,t0,x0,u,f):
        st=x0
        con=u[0,:].T#select first row of control,and than integrate it,1 time rk4
        st=integrator_fun(st,con)
        x0=st.full()
        U0_1 = np.append(u[1:np.size(u,0),:],u[-1,:])#trim the first entry and use the last entry
        t0=t0+DT
        
        return t0,x0,U0_1

def ss(xbar_tilde,t_predict):
        
        xbar_tilde[1] = xbar[1]
        xbar_tilde[2] = xbar[2]
        xbar_tilde[3] = xbar[3]
        xbar_tilde[0] = xbar_tilde[0] + t_predict*(xbar[3]*np.cos(xbar[2] \
        +C1*ubar[1])/(1-kapparef_s(s)*xbar[1]))   
        return xbar_tilde   



args['p']=p
initial=X[:,0]
G=vertcat(initial-P[0:4])
con=U[:,0]

obj=obj+mtimes(mtimes((initial-P[4:8]).T,Q),(initial-P[4:8]))+mtimes(mtimes((con-P[9:11]).T,R),(con-P[9:11]))
initial_next=X[:,1]
initial_next_pred=integrator_fun(initial,con)
G=vertcat(G,initial_next-initial_next_pred)

for k in range (1,N):
        #donot have to change anything over here
    initial=X[:,k]
    con=U[:,k]
    obj=obj+mtimes(mtimes((initial-P[6*k-2:6*k+2]).T,Q),(initial-P[(6*k-2):(6*k+2)]))\
        +mtimes(mtimes((con-P[(6*k)+2:(6*k)+4]).T,R),(con-P[(6*k)+2:(6*k)+4]))
    initial_next=X[:,k+1]
        ##i think the shape is one big than it should be
    initial_next_pred=integrator_fun(initial,con)
    G=vertcat(G,initial_next-initial_next_pred)#21*4
opt_variables=vertcat(reshape(X,n_states*(N+1),1),reshape(U,n_controls*(N),1))
nlp_prob={'f':obj,'x':opt_variables,'g':G,'p':P}    

def mpc_module_augment(i,X_est):
    
    global U0_1
    global u_cl
    global t0  
    global t
    
    
    solver=nlpsol('solver','ipopt',nlp_prob)
    
    if (i==0):
    #xx1=np.zeros((N+1,4,N+1))#for storing prediction equal to length N for every time step
        #x0=X_est[0:4]#relṕpace it by X_est
        xx[:,0]=X_est
        #xx[:,0]=x0

    x0=X_est
    args['x0']=x0
    x0=[x0[0],x0[1],x0[2],x0[3]]
    args['p'][0:4]=x0
    
   
    #initilization for u and x trajectories for an initial state
    #100*2
    x0=np.array(x0)
    X0_1=np.array(repmat(x0,1,N+1))
    X0_1=X0_1.T#101*4


    #simulation loop starts here,but loop is 1

    x0=np.array(x0)
    t_current=i*DT
    #print("the value of t_current is")
    #print(t_current)
   
    for k in range(1,N+1):
        t_predict=t_current+(k-1)*DT
        xbar_tilde = np.zeros((4,))
        xbar_tilde=ss(xbar_tilde,t_predict)
        args['p'][(6*k)-2:(6*k)+2]=xbar_tilde
        args['p'][(6*k)+2:(6*k)+4]=ubar

   
   
    args['x0']=np.append(reshape((X0_1.T),4*(N+1),1),reshape(transpose(U0_1),2*(N),1)).T    
    sol=solver(x0=args['x0'],lbx=args['lbx'],ubx=args['ubx'],lbg=args['lbg'],ubg=args['ubg'],p=args['p'])
    u=reshape((sol['x'][4*(N+1):]).T,2,N).T
    sol['x']=np.array(sol['x']) 
    #xx1[:,0:4,i+1]=reshape((sol['x'][0:4*(N+1)]).T,4,N+1).T#101*4 for each iteration
    u_cl=np.vstack((u_cl,u[0,:]))#first controlt taken
    t0,x0,U0_1=shift(DT,t0,x0,u,integrator_fun)

    xx[0:4,i+1:i+2]=x0
    t=np.append(t,t0)
    #print(u_cl)
    
    return u_cl[-1]  


