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
model_error_plot_est_1=[]
model_error_plot_est_2=[]
model_error_plot_est_3=[]


model_error_plot_true_1=[]
model_error_plot_true_2=[]
model_error_plot_true_3=[]





g_1=MX.sym('g_1')
g_2=MX.sym('g_2')
g_3=MX.sym('g_3')
g=vertcat(g_1,g_2,g_3)



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
xdot = vertcat(sdota+g_1,\
v*sin(alpha+C1*derDelta)+g_2,\
v*C2*derDelta-kapparef[0]*sdota+g_3,\
Fxd/m*cos(C1*derDelta)+g_3,\
0,\
0,\
0)





N=250#preditction horizon(nr.ofsteps)
M=1#1 rk4 step
#dT sampling time
T=4
DT=T/N/M
#################################################function for rk4###########################################################
#to do function for jacobian!

f_rk4=Function('f_rk4',[x,u],[xdot])
xk_jac=jacobian(xdot,x)
xk_jac_val=Function('f_jac_val',[x,u],[xk_jac])




d_value=0.025

#initial conditions,U0 impemented directly in the controls

X_initial=[0.1,0.1,0.1,1,\
0,\
0,\
0,\
]
#U0=[K*np.cos(omega*DT),0.1]

#P_model_cov=diag([1,1,1,1])#3rd diagonal term changed to 1 from 0 because otherwise its inverse twill turn into infinity
#P_model_cov_intial=diag([1,1,1,1,1,1,1])
#P_model_cov_pred=P_model_cov_intial
P_model_cov_est=diag([1,1,1,1,1,1,1])
P_model_cov=diag([1,1,1,1,1,1,1])

W=0.0025
W_cov=W*diag([1,1,1,1,1,1,1])
V=0.0025
V_cov=V*diag([1,1,1])


X_meas=[1,1,1,1,0,0,0]
#X =X_initial
X_est=X_initial


s_plot=np.append(s_plot,X_est[0])
n_plot=np.append(n_plot,X_est[1])
alpha_plot=np.append(alpha_plot,X_est[2])
velocity_plot=np.append(velocity_plot,X_est[3])

s_plot_est=np.append(s_plot_est,X_meas[0])
n_plot_est=np.append(n_plot_est,X_meas[1])
alpha_plot_est=np.append(alpha_plot_est,X_meas[2])
velocity_plot_est=np.append(velocity_plot_est,X_meas[3])
    

                   
for i in range(N):

    vk=transpose(np.random.multivariate_normal(np.zeros(3),V_cov,1))
    wk=transpose(np.random.multivariate_normal(np.zeros(7),W_cov,1))
    # X_est=X_est
    
    A=xk_jac_val(X_est,[K*np.cos(i*omega),0.1])
    #print(A)
    k1 = f_rk4(X_est, [K*np.cos(i*omega),0.1])
    k2 = f_rk4(X_est + DT/2 * k1,[K*np.cos(omega*i),0.1] )
    k3 = f_rk4(X_est + DT/2 * k2, [K*np.cos(i*omega),0.1])
    k4 = f_rk4(X_est + DT * k3, [K*np.cos(i*omega),0.1])
    X_est=X_est+DT/6*(k1 +2*k2 +2*k3 +k4)+wk
    X_predict=X_est
    
    
    #print(X_est)
    
    s_plot=np.append(s_plot,X_predict[0])
    n_plot=np.append(n_plot,X_predict[1])
    alpha_plot=np.append(alpha_plot,X_predict[2])
    velocity_plot=np.append(velocity_plot,X_predict[3])
    model_error_plot_true_1=np.append(model_error_plot_true_1,X_predict[4])
    model_error_plot_true_2=np.append(model_error_plot_true_2,X_predict[5])
    model_error_plot_true_3=np.append(model_error_plot_true_3,X_predict[6])
    

    
    P_model_cov_predict=A*P_model_cov_est*transpose(A)+W_cov
        
    
    C = np.array([[1, 0, 0 ,0 ,0 ,0 ,0], [0, 1, 0 ,0 ,0 ,0 ,0], [0, 0, 1 ,0 ,0 ,0 ,0]])
    
    P_model_cov_est=inv(inv(P_model_cov_predict)+mtimes(mtimes(transpose(C),inv(V_cov)),C))
   
    #print(P_model_cov_est)
    
    k1_y = f_rk4(X_meas, [K*np.cos(i*omega),0.1])
    k2_y = f_rk4(X_meas + DT/2 * k1_y,[K*np.cos(omega*i),0.1] )
    k3_y = f_rk4(X_meas + DT/2 * k2_y, [K*np.cos(i*omega),0.1])
    k4_y = f_rk4(X_meas + DT * k3_y, [K*np.cos(i*omega),0.1])
    X_meas=X_meas+DT/6*(k1_y +2*k2_y +2*k3_y +k4_y)
    #print(X_meas)
    
    #check please
    Y_meas=np.dot(C,X_meas)
    
    Ck_gk=np.dot(C,X_predict)
    
    #print(vk.shape)
    measurement_residual=Y_meas-Ck_gk
    #print(measurement_residual)
    print(P_model_cov_est)
    X_est=X_est+mtimes(mtimes(mtimes(P_model_cov_est,transpose(C)),inv(V_cov)),measurement_residual) 
    #print(X_est)
    #print(X_est)
    #X_pred=X_est
    
    s_plot_est=np.append(s_plot_est,X_est[0])
    n_plot_est=np.append(n_plot_est,X_est[1])
    alpha_plot_est=np.append(alpha_plot_est,X_est[2])
    velocity_plot_est=np.append(velocity_plot_est,X_est[3])
    model_error_plot_est_1=np.append(model_error_plot_est_1,X_est[4])
    model_error_plot_est_2=np.append(model_error_plot_est_2,X_est[5])
    model_error_plot_est_3=np.append(model_error_plot_est_3,X_est[6]) 
     
 



ax1 = plt.subplot2grid(shape=(2,6), loc=(0,0), colspan=2)
plt.ylabel('s_plot')
plt.plot(range(N+1),s_plot,'r')
plt.plot(range(N+1),s_plot_est,'g')

ax2 = plt.subplot2grid((2,6), (0,2), colspan=2)
plt.ylabel('n_plot')
plt.plot(range(N+1),n_plot,'r')
plt.plot(range(N+1),n_plot_est,'g')

ax3 = plt.subplot2grid((2,6), (0,4), colspan=2)
plt.ylabel('alpha_plot')
plt.plot(range(N+1),alpha_plot,'r')
plt.plot(range(N+1),alpha_plot_est,'g')

ax4 = plt.subplot2grid((2,6), (1,1), colspan=2)
plt.ylabel('velocity_plot')
plt.plot(range(N+1),velocity_plot,'r')
plt.plot(range(N+1),velocity_plot_est,'g')


plt.show()  

ax1 = plt.subplot2grid(shape=(2,6), loc=(0,0), colspan=2)
plt.ylabel('model_error_plot_1')
plt.plot(range(N),model_error_plot_true_1,'r')
plt.plot(range(N),model_error_plot_est_1,'g')

ax2 = plt.subplot2grid((2,6), (0,2), colspan=2)
plt.ylabel('model_error_plot_2')
plt.plot(range(N),model_error_plot_true_2,'r')
plt.plot(range(N),model_error_plot_est_2,'g')


ax3 = plt.subplot2grid((2,6), (0,4), colspan=2)
plt.ylabel('model_error_plot_3')
plt.plot(range(N),model_error_plot_true_3,'r')
plt.plot(range(N),model_error_plot_est_3,'g')
plt.show()