from casadi import *
import numpy as np
import sys
sys.path.append("../THESIS")
from common_parameters.common import *
#from model.augmented import g,x_a,u_a,xdot_aU0,X_out_a
from model.augmented import *
from model.non_augmented import integrator_fun
from model.non_augmented import C


#print(DT)
#N_main=110
print(DT)

line_file=np.loadtxt('tracks/LMS_Track.txt')
s0=line_file[:,0]
xref=line_file[:,1]
yref=line_file[:,2]
psiref=line_file[:,3]
kapparef=line_file[:,4]

length=len(s0)
# copy loop to beginning and end
s0=np.append(s0,[s0[length-1] + s0[1:length]])
kapparef=np.append(kapparef,kapparef[1:length])
s0 = np.append([-s0[length-2] + s0[length-81:length-2]],s0)
kapparef = np.append(kapparef[length-80:length-1],kapparef)

# compute spline interpolations
kapparef_s=interpolant('kapparef_s','bspline',[s0],kapparef)


est_true_data=np.loadtxt('configurations/est-true_mpc')
#model_error=np.loadtxt('configurations/model-error-est-true')
u_data=np.loadtxt('configurations/u_data')
d_data=np.loadtxt('configurations/model-error-est-true') 
y_measurements=np.loadtxt('configurations/y_measured_vk')



s_true=est_true_data[:,0]
n_true=est_true_data[:,1]
alpha_true=est_true_data[:,2]
velocity_true=est_true_data[:,3]

s_est=est_true_data[:,4]
n_est=est_true_data[:,5]
alpha_est=est_true_data[:,6]
velocity_est=est_true_data[:,7] 


D_values=u_data[:,0]
Der_values=u_data[:,1]

d_1_init=d_data[0:N_main,3]
d_2_init=d_data[0:N_main,4]
d_3_init=d_data[0:N_main,5]




N_MHE=N_main-1#94

X_true=np.zeros((N_MHE+1,4))
X_est=np.zeros((N_MHE+1,4))
#y_measurements=y_measurements[0:N_MHE+1,:]

print(y_measurements.shape)
print(X_true.shape)
for i in range (N_MHE+1):
  X_true[i,:]=s_true[i],n_true[i],alpha_true[i],velocity_true[i]
  X_est[i,:]=s_est[i],n_est[i],alpha_est[i],velocity_est[i]

#N_MHE=np.size(X_true,0)

x=vertcat(s,n,alpha,v)
g=vertcat(g_1,g_2,g_3)
derD = MX.sym('derD')
derDelta=MX.sym('derDelta')
V=DM([[0.000000025,0,0], [0,0.000000025,0],[0,0,0.000000025]])
#V=DM([[1,0,0], [0,1,0],[0,0,1]])
#V=inv(V)
#decision variables

X=MX.sym('X',n_states,(N_MHE+1))
G_aug=MX.sym('G_aug',3,(N_MHE))
#measurementes
#P=MX.sym('P',3,(N_MHE+1))
P=MX.sym('P',3*(N_MHE+1)+4)#measurements + state for starting
P=P.T
G=[]
print(P.shape)

#objective function formulation

st=X[:,0]
#0:3, 3:6, 6:9

y_tile=P[0:3]
h_x=st[0:3]
#obj=obj+mtimes(mtimes((y_tile-h_x).T,V),(y_tile-h_x))
obj=obj+mtimes(mtimes((y_tile.T-h_x).T,V),(y_tile.T-h_x))#1*3 * 3*3 *3*1

print(P.shape)
for k in range(1,N_MHE+1):
  st=X[:,k]
  
  y_tile=P[:,3*k:3*k+3]
  h_x=st[0:3]
  #obj=obj+mtimes(mtimes((y_tile-h_x).T,V),(y_tile-h_x))
  obj=obj+mtimes(mtimes((y_tile.T-h_x).T,V),(y_tile.T-h_x))#1*3 * 3*3 *3*1


#multiple shooting
st=X[:,0]
G=vertcat(st-P[3*(N_MHE+1):3*(N_MHE+1)+4].T)
#print(G)
for k in range(0,N_MHE):  

   
  st=X[:,k]
  gt=G_aug[:,k]
  x=vertcat(s,n,alpha,v)
  #g=vertcat(g_1,g_2,g_3)
  Fxd=(Cm1-Cm2*v)*D_values[k]-Cr2*v*v-Cr0*tanh(5*v)
  sdota=(v*np.cos(alpha+C1*Der_values[k]))/(1-kapparef_s(s)*n)
 
  xdot_a = vertcat(sdota+g_1,\
  v*sin(alpha+C1*Der_values[k])+g_2,\
  v*C2*Der_values[k]-kapparef_s(s)*sdota+g_3,\
  Fxd/m*cos(C1*Der_values[k]),\
  )
  f_rk4_p=Function('f_rk4_p',[x,g],[xdot_a])
  # model chaning every time wrt to u
  
  X0=MX.sym('X0',4,1)
  G0=MX.sym('G0',3,1)
  k1_p=f_rk4_p(X0,G0)
  k2_p=f_rk4_p(X0+DT/2*k1_p,G0)
  k3_p=f_rk4_p(X0+DT/2*k2_p,G0)
  k4_p=f_rk4_p(X0+DT*k3_p,G0)
  X_out=X0+DT/6*(k1_p+2*k2_p+2*k3_p+k4_p)
  integrator_p=Function('integrator_p',[X0,G0],[X_out])  
  
  st_next=X[:,k+1]
  st_next_euler=integrator_p(st,gt)  
  G=vertcat(G,st_next-st_next_euler)
  
  #print(G)  
  #add one for outside the loop with intergration
#G=vertcat(G,st_next-st_next_euler)  
#print(G)
#G=vertcat(G,st_next-st_next_euler)
OPT_variables = vertcat(reshape(X,4*(N_MHE+1),1),reshape(G_aug,3*(N_MHE),1))#x +d

nlp_mhe={'f':obj,'x':OPT_variables,'g':G,'p':P}
solver=nlpsol('solver','ipopt',nlp_mhe)


#PARAMETER Y _MEASUREMENTS
args={}
y_m=y_measurements.T
y_m=y_m.flatten('F')

p=np.zeros((3*(N_MHE+1)+4))
p[0:3*(N_MHE+2)]=y_m
p[3*(N_MHE+1):3*(N_MHE+1)+1]=0.4
p[3*(N_MHE+1)+1:3*(N_MHE+1)+2]=0.3
p[3*(N_MHE+1)+2:3*(N_MHE+1)+3]=0.3
p[3*(N_MHE+1)+3:3*(N_MHE+1)+4]=0.5

args['p']=p
#args['p']=y_m#3*150
print("p shape is")
print(args['p'])

lbg=np.zeros((4*(N_MHE+1)))
ubg=np.zeros((4*(N_MHE+1)))
ubg[0:4*(N):1]=0
lbg[0:4*(N):1]=0
args['lbg']=lbg
args['ubg']=ubg

ubx=np.zeros((4*(N_MHE+1)))
lbx=np.zeros((4*(N_MHE+1)))
#s
lbx[0:4*(N_MHE+1):4]=-inf
ubx[0:4*(N_MHE+1):4]=+inf
#n
lbx[1:4*(N_MHE+1):4]=-3
ubx[1:4*(N_MHE+1):4]=+3
#alpha
lbx[2:4*(N_MHE+1):4]=-10
ubx[2:4*(N_MHE+1):4]=+10
#v
lbx[3:4*(N_MHE+1):4]=-10
ubx[3:4*(N_MHE+1):4]=+10

lbx=np.append(lbx,np.zeros((3*N_MHE)))
ubx=np.append(ubx,np.zeros((3*N_MHE)))
#d
lbx[4*N_MHE+1:]=-10
ubx[4*N_MHE+1:]=10

args['lbx']=lbx
args['ubx']=ubx





###simulation####
#initilization of state decision variables
X0_1=np.zeros(((N_MHE+1,4)))#150*4
G0_1=np.zeros((N_MHE,3))#149*3



X0_1[:,0]=X_est[:,0]
X0_1[:,1]=X_est[:,1]
X0_1[:,2]=X_est[:,2]
X0_1[:,3]=X_est[:,3]

G0_1[:,0]=d_1_init[:-1]
G0_1[:,1]=d_2_init[:-1]
G0_1[:,2]=d_3_init[:-1]


#print(d_1_init[:-1].shape)

#print(G0_1)

#initilization value of the optimization problem
args['x0']=np.append(reshape((X0_1.T),4*(N_MHE+1),1),reshape((G0_1.T),3*(N_MHE),1))
#print(args['x0'])

sol=solver(x0=args['x0'],p=args['p'],lbg=args['lbg'],ubg=args['ubg'],lbx=args['lbx'],ubx=args['ubx'])
X_sol=reshape(sol['x'][0:4*(N_MHE+1)].T,4,N_MHE+1).T
d_sol=reshape(sol['x'][4*(N_MHE+1):].T,3,N_MHE).T
d_sol[0,:]=[0.1,0.1,0.1]



s_lls=(X_sol[:,0]).full()
s_lls=s_lls.flatten()

n_lls=(X_sol[:,1]).full()
n_lls=n_lls.flatten()


alpha_lls=(X_sol[:,2]).full()
alpha_lls=alpha_lls.flatten()


vel_lls=(X_sol[:,3]).full()
vel_lls=vel_lls.flatten()

d1_lls=(d_sol[:,0]).full()
d1_lls=d1_lls.flatten()

d2_lls=(d_sol[:,1]).full()
d2_lls=d2_lls.flatten()

d3_lls=(d_sol[:,2]).full()
d3_lls=d3_lls.flatten()


#s=s.flatten()

plt.figure(1)

plt.xlabel('number of total points ')
plt.ylabel('s_estimated')
plt.plot(range(N_MHE+1),s_lls,label='s_mhe')
#plt.plot(range(N_MHE+1),s_est[:-1],label='ground truth')
plt.plot(range(N_MHE),y_measurements[:-2,0],label='y_measurement_s')
plt.legend()
plt.show()

plt.figure(2)

plt.xlabel('number of total points ')
plt.ylabel('n_estimated')
plt.plot(range(N_MHE+1),n_lls,label='n_mhe')
#plt.plot(range(N_MHE+1),y_measurements[:,1],label='n_measured')
plt.plot(range(N_MHE+2),y_measurements[:,1],label='y_measurement_n')
plt.legend()
plt.show()
plt.figure(3)

plt.xlabel('number of total points ')
plt.ylabel('alpha_estimated')
plt.plot(range(N_MHE+1),alpha_lls,label='alpha_mhe')
#plt.plot(range(N_MHE+1),y_measurements[:,2],label='mesured alpha')
plt.plot(range(N_MHE+2),y_measurements[:,2],label='y_measurement_alpha')
plt.legend()

plt.show()


plt.figure(4)

plt.xlabel('s reference ')
plt.ylabel('vel_estimated')
plt.plot(range(N_MHE+1),vel_lls,label='vel_estimated')
plt.plot(range(N_MHE+1),velocity_true[0:N_main],label='true velocity')
plt.show()


plt.figure(5)

plt.xlabel('s reference ')
plt.ylabel(r'$d_s$')
plt.plot(s0[0:N_main-1],d1_lls,label='offline')
plt.plot(s0[0:N_main],d_1_init,label='true')
plt.legend()
plt.savefig('plots_save/disturbed/offline/disturbance_parameters/d_s.pdf')
plt.show()


plt.figure(6)

plt.xlabel('s reference ')
plt.ylabel(r'$d_n$')
plt.plot(s0[0:N_main-1],d2_lls,label='offline')
plt.plot(s0[0:N_main],d_2_init,label='true')
plt.legend()
plt.savefig('plots_save/disturbed/offline/disturbance_parameters/d_n.pdf')
plt.show()


plt.figure(7)

plt.xlabel('s reference ')
plt.ylabel(r'$d_\alpha$')
plt.plot(s0[0:N_main-1],d3_lls,label='offline')
plt.plot(s0[0:N_main],d_3_init,label='true')
plt.legend()
plt.savefig('plots_save/disturbed/offline/disturbance_parameters/d_alpha.pdf')
plt.show() 


""" plt.figure(8)
plt.xlabel('numer of total points')
plt.ylabel('s')
plt.plot(range(N_MHE+1),s_true[:-1],label='plant s')
plt.plot(range(N_MHE+1),s_est[:-1],label='estimated s')
plt.plot(range(N_MHE+1),y_measurements[:,0],label='y_measurement')
plt.legend()
plt.show()

plt.figure(9)
plt.xlabel('numer of total points')
plt.ylabel('n')
plt.plot(range(N_MHE+1),n_true[:-1],label='plant n')
plt.plot(range(N_MHE+1),n_est[:-1],label='estimated n')
plt.plot(range(N_MHE+1),y_measurements[:,1],label='y_measurement')
plt.legend()
plt.show()

plt.figure(10)
plt.xlabel('numer of total points')
plt.ylabel('alpha')
plt.plot(range(N_MHE+1),alpha_true[:-1],label='true_alpha')
plt.plot(range(N_MHE+1),alpha_est[:-1],label='estimated_alpha')
plt.plot(range(N_MHE+1),y_measurements[:,2],label='y_measurement')
plt.legend()
"""
plt.show()

np.savetxt('configurations/d_estimated_values',d_sol)
np.savetxt('configurations/x_estimated_values',X_sol)

   

