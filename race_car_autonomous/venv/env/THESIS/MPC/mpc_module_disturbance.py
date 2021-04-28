import sys
sys.path.append("../THESIS")
from casadi import *
import numpy as np
import math
import matplotlib.pyplot as plt
from model.non_augmented import *
from model.non_augmented import f_rk4 as f_rk4
from model.non_augmented import integrator_fun as integrator_fun 
from model.augmented import integrator_fun_a,integrator_fun_mpc
from common_parameters.common import *
from common_parameters.MPC_common import *
from model.augmented import d_1,d_2,d_3,d_1_p,d_2_p,d_3_p,u_a

##implement in mpc_common
#Q=np.diag ([0.1, 100, 1, 2])
#R=np.diag([0.05 ,0.1]) 
from collections import namedtuple

xbar=np.array([0,\
    0,\
    0,\
    1])

ubar = np.array([0.08138923, 0.091973])

t=0
t0=0    
t_current=1
u_cl=[0,0]

t_predict=0
U0_1=np.zeros((N,2))
P=MX.sym('P',n_states+N*(n_states+n_controls)+3)
p=np.zeros(n_states+N*(n_controls+n_states)+3)


def shift(DT,t0,x0,u,integrator_fun_mpc,d_mpc):
        initial=x0
        con=u[0,:].T#select first row of control,and than integrate it,1 time rk4
        
        initial=integrator_fun_mpc(initial,con,d_mpc)
        x0=initial.full()
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
con=U[:,0]
initial=X[:,0]

G=vertcat(initial-P[0:4])
for k in range (0,N):
    ####check it
    if(k==0):
        disturbance=P[n_states+N*(n_controls+n_states):n_states+N*(n_controls+n_states)+3]
        #con=U[:,0]
        obj=obj+mtimes(mtimes((initial-P[6*1-2:6*1+2]).T,Q),(initial-P[(6*1-2):(6*1+2)]))\
        +mtimes(mtimes((con-P[(6*1)+2:(6*1)+4]).T,R),(con-P[(6*1)+2:(6*1)+4]))
            #+\
            #mtimes(mtimes((zk).T,1000000),(zk))
        initial_next=X[:,1]
        initial_next_predict=integrator_fun_mpc(initial,con,disturbance)
        G=vertcat(G,initial_next-initial_next_predict)
    elif(k>0):
        disturbance=P[n_states+N*(n_controls+n_states):n_states+N*(n_controls+n_states)+3]    
        initial=X[:,k]
        con=U[:,k]
        #zk=z_k[:,k]
        obj=obj+mtimes(mtimes((initial-P[6*k+6-2:6*k+2+6]).T,Q),(initial-P[(6*k-2+6):(6*k+2+6)]))\
            +mtimes(mtimes((con-P[(6*k)+2+6:(6*k)+4+6]).T,R),(con-P[(6*k)+2+6:(6*k)+4+6]))
            #+\
            #mtimes(mtimes((zk).T,1000000),(zk))
    
        initial_next=X[:,k+1]
        initial_next_predict=integrator_fun_mpc(initial,con,disturbance)
        
        G=vertcat(G,initial_next-initial_next_predict)
  
#G=vertcat(G,initial_next-initial_next_predict)
#print(G)
#print(P.shape)

""" initial=X[:,k]
con=U[:,k]
obj=obj+mtimes(mtimes((initial-P[6*k-2:6*k+2]).T,Q),(initial-P[(6*k-2):(6*k+2)]))\
    +mtimes(mtimes((con-P[(6*k)+2:(6*k)+4]).T,R),(con-P[(6*k)+2:(6*k)+4]))
initial_next=X[:,k+1]
initial_next_pred=integrator_fun(initial,con)
G=vertcat(G,initial_next-initial_next_pred)#21*4
"""


options = { 'ipopt': {'bound_relax_factor':1e-4,
            'constr_viol_tol': 1e-4
            
                }}



####changed here
opt_variables=vertcat(reshape(X,n_states*(N+1),1),reshape(U,n_controls*(N),1))
nlp_prob={'f':obj,'x':opt_variables,'g':G,'p':P}    
#ipopt_bound.bound_relax_factor

solver=nlpsol('solver','ipopt',nlp_prob,options)

def cost_fun(x_optimum,xr,ur,u_cl):
    u_opt=u_cl[-1]
    xr=xr.T
    ur=ur.T
    print(x_optimum.shape)
    x_opt=transpose(x_optimum[[0],:])
    u_opt=np.expand_dims(u_opt,axis=1)
    result=np.dot(np.dot(transpose(x_opt-xr),Q),(x_opt-xr))+np.dot(np.dot(transpose(u_opt-ur),R),(u_opt-ur))
    print(result)
    return result
    


def mpc_module_disturbance(i,X_est,N_main):
    
    global U0_1
    global u_cl
    global t0  
    global t
    
    
    
    #solver=nlpsol('solver','ipopt',nlp_prob)
    
    if (i==0):
        xx[:,0]=X_est[0:4]
        
    
    
    x0=X_est[0:4]#initial condition
    args['x0']=x0
    x0=[x0[0],x0[1],x0[2],x0[3]]

    args['p'][0:4]=x0#initial condition 
    
    d_mpc=np.array((X_est[4],X_est[5],X_est[6]))
    args['p'][n_states+N*(n_controls+n_states):n_states+N*(n_controls+n_states)+3]=d_mpc    
    
    x0=np.array(x0)
    X0_1=np.array(repmat(x0,1,N+1))
    X0_1=X0_1.T#101*4
    #simulation loop starts here,but loop is 1
    
    #x0=np.array(x0)
    t_current=i*DT+1

    
    t_predict=t_current-DT
    xbar_tilde = np.zeros((4,))
    xbar_tilde=ss(xbar_tilde,t_predict)
    args['p'][4:8]=xbar_tilde
    args['p'][8:10]=ubar
   
    
   
    x_ref=xbar_tilde
    u_ref=ubar
    for k in range(1,N):
        t_predict=t_current+(k-1)*DT
        xbar_tilde = np.zeros((4,))
        xbar_tilde=ss(xbar_tilde,t_predict)
        args['p'][(6*k)+4:(6*k)+8]=xbar_tilde
        args['p'][(6*k)+8:(6*k)+10]=ubar



    


          
    print("value of disturbance is")
    print(d_mpc)
    print("last 3 values should be disturbance added parameterized")
    print(args['p']) 
      
    
    args['x0']=np.append(reshape((X0_1.T),4*(N+1),1),reshape(transpose(U0_1),2*(N),1)).T    
    sol=solver(x0=args['x0'],lbx=args['lbx'],ubx=args['ubx'],lbg=args['lbg'],ubg=args['ubg'],p=args['p'])
    u=reshape((sol['x'][4*(N+1):]).T,2,N).T
    x_optimum=reshape((sol['x'][0:4*(N+1)]).T,4,N+1).T

    sol['x']=np.array(sol['x']) 
    sol['f']=np.array(sol['f'])
    #xx1[:,0:4,i+1]=reshape((sol['x'][0:4*(N+1)]).T,4,N+1).T#101*4 for each iteration
    u_cl=np.vstack((u_cl,u[0,:]))#first controlt taken

    current_cost=sol['f']
    
    

    
       
    t0,x0,U0_1=shift(DT,t0,x0,u,integrator_fun_mpc,d_mpc)

    xx[0:4,i+1:i+2]=x0

    xr=np.array(args['p'][4:8])
    xr=np.expand_dims(xr,axis=1)
    xr=xr.T
    ur=np.array(args['p'][8:10])
    ur=np.expand_dims(ur,axis=1)
    ur=ur.T
    status=solver.stats()['return_status']
    t=np.append(t,t0)

    result=cost_fun(x_optimum,xr,ur,u_cl)



    return (u_cl[-1],current_cost,status,x_ref,ubar,result)
        
 
