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
s_plot=[]
n_plot=[]
alpha_plot=[]
velocity_plot=[]
K=0.2
omega=0.2
        


#dealing with 4 states for now
s   = MX.sym('s')
n   = MX.sym('n')
alpha = MX.sym('alpha')
v   = MX.sym('v')
#D = MX.sym('D')
#delta=MX.sym('delta')
x = vertcat(s,n,alpha,v)


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
Fxd/m*cos(C1*derDelta))


N=250#preditction horizon(nr.ofsteps)
M=1#1 rk4 step
#dT sampling time
T=1
DT=T/N/M
#################################################function for rk4###########################################################
#to do function for jacobian!

f_rk4=Function('f_rk4',[x,u],[xdot])
xk_jac=jacobian(xdot,x)
xk_jac_val=Function('f_jac_val',[x,u],[xk_jac])




#f_jac=Function('f_jac',)
#initial conditions,U0 impemented directly in the controls
X0=[0.1,0.1,0.1,0.1]
U0=[K*np.cos(omega*DT),0]

#print(xk_jac_val.size())

P_model_cov=diag([0.005,0.1,0.001,0.001])#3rd diagonal term changed to 1 from 0 because otherwise its inverse twill turn into infinity


for i in range(N):
    #f_jac_val(0.1)    
    #save for plotting
    s_plot=np.append(s_plot,X0[0])
    n_plot=np.append(n_plot,X0[1])
    alpha_plot=np.append(alpha_plot,X0[2])
    velocity_plot=np.append(velocity_plot,X0[3])
    
    
    A=xk_jac_val(X0,[K*np.cos(i*omega),0.1])
    k1 = f_rk4(X0, [K*np.cos(i*omega),0.1])
    k2 = f_rk4(X0 + DT/2 * k1,[K*np.cos(omega*i),0.1] )
    k3 = f_rk4(X0 + DT/2 * k2, [K*np.cos(i*omega),0.1])
    k4 = f_rk4(X0 + DT * k3, [K*np.cos(i*omega),0.1])
    X0=X0+DT/6*(k1 +2*k2 +2*k3 +k4)

    P_model_cov=A*P_model_cov
    print(A)

        
""" plt.figure()        
plt.subplot(221)
#print(range(N))
#print(s_plot)
 plt.ylabel('s_plot')
plt.plot(range(N),s_plot,'r--')
plt.subplot(222)
plt.ylabel('n_plot')
plt.plot(range(N),n_plot,'g--')
plt.subplot(223)
plt.ylabel('alpha_plot')
plt.plot(range(N),alpha_plot,'b--')
plt.subplot(224)
plt.ylabel('velocity_plot')
plt.plot(range(N),velocity_plot,'y--')
plt.show()
 """         








