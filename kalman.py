import numpy as np
from casadi import *
import matplotlib.pyplot as plt


# load track parameters
line_file=np.loadtxt('line.txt')


#[s0,xref,yref,psiref,kapparef]=line_file[0]
s0=line_file[:,0]
xref=line_file[:,1]
yref=line_file[:,2]
psiref=line_file[:,3]
kapparef=line_file[:,4]

length=len(s0)
pathlength=s0[-1]
#copy loop to beginning and end
s0=np.append(s0,[s0[length-1] + s0[1:length]])
kapparef=np.append(kapparef,kapparef[1:length])
s0 = np.append([-s0[length-2] + s0[length-81:length-2]],s0)
kapparef = np.append(kapparef[length-80:length-1],kapparef)




## Race car parameters
m = 0.043
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

s_plot=[]
n_plot=[]
alpha_plot=[]
velocity_plot=[]
s_plot_est=[]
n_plot_est=[]
alpha_plot_est=[]
velocity_plot_est=[]
model_error_plot_est=[]


g_1=MX.sym('g_1')
g_2=MX.sym('g_2')
g_3=MX.sym('g_3')
g_4=MX.sym('g_4')
g=vertcat(g_1,g_2,g_3,g_4)



#dealing with 4 states for now
s   = MX.sym('s')
n   = MX.sym('n')
alpha = MX.sym('alpha')
v   = MX.sym('v')
#D = MX.sym('D')
#delta=MX.sym('delta')
x = vertcat(s,n,alpha,v,g)


 # controls
derD = MX.sym('derD')
derDelta=MX.sym('derDelta')
u = vertcat(derD,derDelta)

#parameters needed for the ode
Fxd=(Cm1-Cm2*v)*derD-Cr2*v*v-Cr0*tanh(5*v)
sdota=(v*np.cos(alpha+C1*derDelta))/(1-kapparef[0]*n)


# model in race track coordinates
xdot = vertcat(sdota,\
v*sin(alpha+C1*derDelta),\
v*C2*derDelta-kapparef[0]*sdota,\
Fxd/m*cos(C1*derDelta),\
g)



ydot=vertcat(vertcat(sdota,\
v*sin(alpha+C1*derDelta),\
v*C2*derDelta-kapparef[0]*sdota),\
g)

N=100#preditction horizon(nr.ofsteps)
M=1#1 rk4 step
#dT sampling time
T=1
DT=T/N/M
#################################################function for rk4###########################################################
#to do function for jacobian!

f_rk4=Function('f_rk4',[x,u],[xdot])
xk_jac=jacobian(xdot,x)
xk_jac_val=Function('f_jac_val',[x,u],[xk_jac])

yk_jac=jacobian(ydot,x)
#print(yk_jac.size())
yk_jac_val=Function('yk_jac_val',[x,u],[yk_jac])
#print(yk_jac_val)



mean,sigma=0,0.5

#initial conditions,U0 impemented directly in the controls
#X0=[0.1,0.1,0.1,6]#chose the value of 6 because it is giving better results for the covariance matrix
X0=[0.1,0.1,0.1,6,\
np.random.normal(mean,sigma,1),\
np.random.normal(mean,sigma,1),\
np.random.normal(mean,sigma,1),\
np.random.normal(mean,sigma,1)]
U0=[K*np.cos(omega*DT),0]

#P_model_cov=diag([1,1,1,1])#3rd diagonal term changed to 1 from 0 because otherwise its inverse twill turn into infinity
P_model_cov=diag([1,1,1,1,1,1,1,1])

Wk=np.random.normal(mean,sigma,1)
Vk=np.random.normal(mean,sigma,1)
#print(Wk)
#Wk=Wk*diag([1,1,1,1])
Wk=Wk*diag([1,1,1,1,1,1,1,1])
#Vk=Vk*diag([1,1,1])
Vk=Vk*diag([1,1,1,1,1,1,1])
#print(Wk)
s_plot=np.append(s_plot,X0[0])
n_plot=np.append(n_plot,X0[1])
alpha_plot=np.append(alpha_plot,X0[2])
velocity_plot=np.append(velocity_plot,X0[3])
s_plot_est=np.append(s_plot_est,X0[0])
n_plot_est=np.append(n_plot_est,X0[1])
alpha_plot_est=np.append(alpha_plot_est,X0[2])
velocity_plot_est=np.append(velocity_plot_est,X0[3])
    
    
for i in range(N):
    #f_jac_val(0.1)    
    #save for plotting
    
    
    A=xk_jac_val(X0,[K*np.cos(i*omega),0.1])
    #print(A)
    k1 = f_rk4(X0, [K*np.cos(i*omega),0.1])
    k2 = f_rk4(X0 + DT/2 * k1,[K*np.cos(omega*i),0.1] )
    k3 = f_rk4(X0 + DT/2 * k2, [K*np.cos(i*omega),0.1])
    k4 = f_rk4(X0 + DT * k3, [K*np.cos(i*omega),0.1])
    X0=X0+DT/6*(k1 +2*k2 +2*k3 +k4)
    #print("old values are")
    #print(X0)
    s_plot=np.append(s_plot,X0[0])
    n_plot=np.append(n_plot,X0[1])
    alpha_plot=np.append(alpha_plot,X0[2])
    velocity_plot=np.append(velocity_plot,X0[3])
    
    P_model_cov=A*P_model_cov*transpose(A)+Wk
    #print(P_model_cov)

    #innovation steps    
    #s,n and alpha are the measured states!
    
    C=yk_jac_val(X0,[K*np.cos(i*omega),0.1])
    #print(C)
    #print(Vk.size())  
    P_model_cov=inv(inv(P_model_cov)+mtimes(mtimes(transpose(C),inv(Vk)),C))
   # print(P_model_cov)
    
    #wk=np.random.normal(mean,sigma,1)*[1,1,1]
    wk=np.random.normal(mean,sigma,1)*[1,1,1,1,1,1,1]
    gk=vertcat(X0[0:3,:],X0[4:8,:])
    #print(gk)
    #print(wk)
    Y0=gk+wk
    measurement_residual=Y0-gk
    #print(measurement_residual)
    X0=X0+mtimes(mtimes(mtimes(P_model_cov,transpose(C)),inv(Vk)),measurement_residual)
    #print("new values are")
    #print(X0)
    #print(mtimes(mtimes(mtimes(P_model_cov,transpose(C)),inv(Vk)),measurement_residual))
    s_plot_est=np.append(s_plot_est,X0[0])
    n_plot_est=np.append(n_plot_est,X0[1])
    alpha_plot_est=np.append(alpha_plot_est,X0[2])
    velocity_plot_est=np.append(velocity_plot_est,X0[3])
    model_error_plot_est=np.append(model_error_plot_est,X0[6])
 
   
 

""" plt.figure()        
plt.subplots(221)
#print(range(N))
#print(s_plot)
plt.ylabel('s_plot')
plt.plot(range(N+1),s_plot,'r--')
plt.plot(range(N+1),s_plot_est,'g--')

plt.subplot(222)
plt.ylabel('n_plot')
plt.plot(range(N+1),n_plot,'r--')
plt.plot(range(N+1),n_plot_est,'g--')
plt.subplot(223)
plt.ylabel('alpha_plot')
plt.plot(range(N+1),alpha_plot,'r--')
plt.plot(range(N+1),alpha_plot_est,'g--')
plt.subplot(224)
plt.ylabel('velocity_plot')
plt.plot(range(N+1),velocity_plot,'r--')
plt.plot(range(N+1),velocity_plot_est,'g--')


#plt.plot      
plt.ylabel('model_error_plot')
plt.plot(range(N),model_error_plot_est,'r--')
#plt.plot(range(N+1),velocity_plot_est,'g--')
plt.show()
 """  
ax1 = plt.subplot2grid(shape=(2,6), loc=(0,0), colspan=2)
plt.ylabel('s_plot')
plt.plot(range(N+1),s_plot,'r--')
plt.plot(range(N+1),s_plot_est,'g--')

ax2 = plt.subplot2grid((2,6), (0,2), colspan=2)
plt.ylabel('n_plot')
plt.plot(range(N+1),n_plot,'r--')
plt.plot(range(N+1),n_plot_est,'g--')

ax3 = plt.subplot2grid((2,6), (0,4), colspan=2)
plt.ylabel('alpha_plot')
plt.plot(range(N+1),alpha_plot,'r--')
plt.plot(range(N+1),alpha_plot_est,'g--')

ax4 = plt.subplot2grid((2,6), (1,1), colspan=2)
plt.ylabel('velocity_plot')
plt.plot(range(N+1),velocity_plot,'r--')
plt.plot(range(N+1),velocity_plot_est,'g--')

ax5 = plt.subplot2grid((2,6), (1,3), colspan=2)
plt.ylabel('model_error_plot')
plt.plot(range(N),model_error_plot_est,'r--')
#plt.plot(range(N+1),velocity_plot_est,'g--')

plt.show()

 



  