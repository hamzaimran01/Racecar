import sys
sys.path.append("../THESIS")
from casadi import *
import numpy as np
import math
import matplotlib.pyplot as plt
from common_parameters.common import *
from common_parameters.MPC_common import *
from model.non_augmented import *



Q=np.diag ([0.1, 100, 1, 2])
R=np.diag([0.05 ,0.1])

initial=X[:,0]
G=vertcat(initial-P[0:4])
con=U[:,0]
#g is filled up over here

obj=obj+mtimes(mtimes((initial-P[4:8]).T,Q),(initial-P[4:8]))+mtimes(mtimes((con-P[9:11]).T,R),(con-P[9:11]))
initial_next=X[:,1]

initial_next_pred=integrator_fun(initial,con)
print(initial_next)
G=vertcat(G,initial_next-initial_next_pred)


for k in range (1,N):

    initial=X[:,k]
    con=U[:,k]
    obj=obj+mtimes(mtimes((initial-P[6*k-2:6*k+2]).T,Q),(initial-P[(6*k-2):(6*k+2)]))\
        +mtimes(mtimes((con-P[(6*k)+2:(6*k)+4]).T,R),(con-P[(6*k)+2:(6*k)+4]))
    initial_next=X[:,k+1]
    ##i think the shape is one big than it should be
    initial_next_pred=integrator_fun(initial,con)
    G=vertcat(G,initial_next-initial_next_pred)#21*4






#4*21,2*20=124,x and u optimization variables,1 column vector
opt_variables=vertcat(reshape(X,n_states*(N+1),1),reshape(U,n_controls*(N),1))
nlp_prob={'f':obj,'x':opt_variables,'g':G,'p':P}
#settings nlp solver
#opts={'max_iter':2000,'print_level':0,'print_time':0,'acceptable_tol':1e-8}
solver=nlpsol('solver','ipopt',nlp_prob)



#simulation loop
t=0#series of values
t0=0#operating point
#xx1=np.zeros((N+1,4,N+1))#for storing prediction equal to length N for every time step
u_cl=ubar##check this out!!
#p=np.zeros(n_states+N*(n_controls+n_states))
#x0=np.zeros((4,1))

    
#x0=X_est
x0=np.array([0,0.2,0.2,2])

xx[:,0]=x0

args['x0']=x0
args['p']=p
args['p'][0:4]=x0

#initilization for u and x trajectories for an initial state
U0_1=np.zeros((N,2))#100*2
X0_1=np.array(repmat(x0,1,N+1))

X0_1=X0_1.T#101*4
count=0



def shift(DT,t0,x0,u,f):
    st=x0
    con=u[0,:].T#select first row of control,and than integrate it,1 time rk4
    st=integrator_fun(st,con)
    x0=st.full()

    t0=t0+DT
    

    U0_1 = np.append(u[1:np.size(u,0),:],u[-1,:])#trim the first entry and use the last entry
    
    return t0,x0,U0_1

def ss(xbar_tilde):
    
    xbar_tilde[1] = xbar[1]
    xbar_tilde[2] = xbar[2]
    xbar_tilde[3] = xbar[3]
    xbar_tilde[0] = xbar_tilde[0] + t_predict*(xbar[3]*np.cos(xbar[2] \
    +C1*ubar[1])/(1-kapparef_s(s)*xbar[1]))
    return xbar_tilde   

count=0
t_current=0
for i in range (1,N):

    x0=np.array(x0)
    
    args['p'][0:4]=x0.ravel()
    count=0    
    t_current=i*DT
    print(t_current)
    for k in range(1,N+1):
        t_predict=t_current+(k-1)*DT
        xbar_tilde = np.zeros((4,))
        xbar_tilde=ss(xbar_tilde)
        
        args['p'][(6*k)-2:(6*k)+2]=xbar_tilde
        args['p'][(6*k)+2:(6*k)+4]=ubar
        count=count+1
        
       
    #args['p'][118:122]=xbar_tilde
    #args['p'][122:124]=ubar    
    args['x0']=np.append(reshape((X0_1.T),4*(N+1),1),reshape(transpose(U0_1),2*(N),1)).T#124
    
    #args has lbg,ubg,lbx,ubx,,p
    sol=solver(x0=args['x0'],lbx=args['lbx'],ubx=args['ubx'],lbg=args['lbg'],ubg=args['ubg'],p=args['p'])
    u=reshape((sol['x'][4*(N+1):]).T,2,N).T
    sol['x']=np.array(sol['x']) 
    xx1[:,0:4,i+1]=reshape((sol['x'][0:4*(N+1)]).T,4,N+1).T#101*4 for each iteration
    u_cl=np.vstack((u_cl,u[0,:]))#first controlt taken
   
    t0,x0,U0_1=shift(DT,t0,x0,u,integrator)
  
    xx[0:4,i:i+1]=x0
    t=np.append(t,t0)
  
    X0_1 =reshape(sol['x'][0:4*(N+1)].T,4,N+1).T
    X0_1 = np.append(X0_1[1:,:],X0_1[-1,:])#new initilization 1 increment
   
    
#t=np.append(t,DT)    

#t=t.flatten()
#plt.legend((s_values,n_values,alpha,v))
print(xx)
plt.step(t,u_cl[:,0],label="derD")
plt.step(t,u_cl[:,1],label="derDelta")
plt.legend()
plt.show()


s_plot=plt.plot(t,xx[0,:-1],label="s_values")
n_plot=plt.plot(t,xx[1,:-1],label="n_values")
alpha_plot=plt.plot(t,xx[2,:-1],label="alpha_values")
vel_plot=plt.plot(t,xx[3,:-1],label="v_values")
plt.legend()
 
plt.show() 

print(u_cl)





    











