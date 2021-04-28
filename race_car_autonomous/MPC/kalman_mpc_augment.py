
import sys
sys.path.append("../THESIS")

import imp
from common_parameters.common import *
from plots import *
import numpy as np
from casadi import *
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pgf import PdfPages
from coord_transform.spattoOrig import *
from scipy.optimize import root
from MPC.mpc_module_disturbance import *
#from MPC.mpc_module_ndisturbance import *
import model.non_augmented as mod_nonaug
import model.augmented as mod_aug
import time





module_name='MPC.mpc_module_ndisturbance'
if module_name  in sys.modules:
    mpc=mpc_module_ndisturbance
    cost_flag='non_disturbed'
    print("imported non disturbed model mpc")

module_name='MPC.mpc_module_disturbance'
if module_name in sys.modules:
    mpc=mpc_module_disturbance
    cost_flag='disturbed'
    print("imported disturbed model mpc")




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


t_current=0 
C=mod_aug.C
integrator_fun_a=mod_aug.integrator_fun_a
X_initial=mod_aug.X_initial
X_est=mod_aug.X_est
integrator_fun=mod_nonaug.integrator_fun

u_D_plot=[]
u_delta_plot=[]
X_true =X_initial
#P_model_cov=diag([1,1,1,1,10,10,10])
#P_model_cov=diag([1,1,1,1,5,5,5])
W=0.000025
W_cov=W*diag([1, 1, 1, 1, 10, 10, 10])
P_model_cov_est=W
V=0.00025
V_cov=V*diag([1,1,1])


xbar_fun = integrator_fun(xbar,ubar)
s_tilde_bar = xbar_fun[0] - \
    (DT*(xbar[3]*np.cos(xbar[2] \
    +C1*ubar[1])/(1-kapparef[0]*xbar[1])))

if np.linalg.norm(s_tilde_bar) > 1.0e-3:
    raise Exception('computed steady-state ifs not a steate-state!')

#plotting initilization
s_plot=np.append(s_plot,X_true[0])
n_plot=np.append(n_plot,X_true[1])
alpha_plot=np.append(alpha_plot,X_true[2])
velocity_plot=np.append(velocity_plot,X_true[3])
s_plot_est=np.append(s_plot_est,X_est[0])
n_plot_est=np.append(n_plot_est,X_est[1])
alpha_plot_est=np.append(alpha_plot_est,X_est[2])
velocity_plot_est=np.append(velocity_plot_est,X_est[3])

s_plot_meas_y=[]
n_plot_meas_y=[]
alpha_plot_meas_y=[]



def plant_dyn_vary(X_true,u,A,P_model_cov_est,i):
    P_model_cov_true=mtimes(A, mtimes(P_model_cov_est, \
        transpose(A))) + W_cov
    X_true=integrator_fun_a(X_true,u)+wk
    #print(X_true)
    

    X_true[4]=d_1_p[i]
    X_true[5]=d_2_p[i]
    X_true[6]=d_3_p[i]
   
    return (X_true,P_model_cov_true)


def plant_dyn_const(X_true,u,A,P_model_cov_est,i):
    P_model_cov_true=mtimes(A, mtimes(P_model_cov_est, \
        transpose(A))) + W_cov
    X_true=integrator_fun_a(X_true,u)+wk

    
    if(i<125):
        X_true[4]=-0.25
        X_true[5]=-0.25
        X_true[6]=-0.25

    if(i>=125 and i<250):
        X_true[4]=0.25
        X_true[5]=0.25
        X_true[6]=0.25

    if(i>=250 and i<375):
        X_true[4]=-0.25
        X_true[5]=-0.25
        X_true[6]=-0.25


    if(i>=375 and i<500):
        X_true[4]=0.25
        X_true[5]=0.25
        X_true[6]=0.25

    


    """ X_true[4]=0.25
    X_true[5]=0.25
    X_true[6]=0.25
    """    #print("the value of N_main is")
    #print(i)


    return (X_true,P_model_cov_true)


def observer_dyn(P_model_cov_true,X_est,u):
   
    P_model_cov_est=inv(inv(P_model_cov_true)+\
        mtimes(mtimes(transpose(C),inv(V_cov)),C))
    X_est=mod_aug.integrator_fun_a(X_est,u)
    
    #reduced the value of vk
    Y_meas=X_true[0:3,:]+vk
    
    Ck_gk=np.dot(C,X_est)
    measurement_residual=Y_meas-Ck_gk
    Kalman_gain = mtimes(mtimes(P_model_cov_est,transpose(C)),inv(V_cov))
    X_est=X_est+mtimes(Kalman_gain,measurement_residual)
     
    return (X_est,P_model_cov_est,Y_meas)


def observer_dyn_offline(P_model_cov_true,X_est,u):
   
    P_model_cov_est=inv(inv(P_model_cov_true)+\
        mtimes(mtimes(transpose(C),inv(V_cov)),C))
    X_est=mod_aug.integrator_fun_a(X_est,u)
    
    #reduced the value of vk
    Y_meas=X_true[0:3,:]+vk
    
    Ck_gk=np.dot(C,X_est)
    measurement_residual=Y_meas-Ck_gk
    Kalman_gain = mtimes(mtimes(P_model_cov_est,transpose(C)),inv(V_cov))
    X_est=X_est+mtimes(Kalman_gain,measurement_residual)
    
    
    X_est[4]=d_1[i]
    X_est[5]=d_2[i]
    X_est[6]=d_3[i]
   
       
    return (X_est,P_model_cov_est,Y_meas)


tic=time.perf_counter()
#print(tic)
new_cost=0
cost_list=[]
cost=0
x_reference=np.zeros((N_main,4))
u_reference=np.zeros((N_main,2))
status=None
count=0

plant_model='None'
observer_mode='None'
observer_model='None'

for i in range(N_main):
    
    
    vk=transpose(np.random.multivariate_normal(np.zeros(3),V_cov,1))
    wk=transpose(np.random.multivariate_normal(np.zeros(7),W_cov,1))

    u,current_cost,status,x_r,ubar,result=mpc(i,X_est,N_main)
    

    
    x_reference[i,:]=x_r
    u_reference[i,:]=ubar

    new_cost=new_cost+result
    #cost_increase=cost_list+result

    cost_list=np.append(cost_list,new_cost)
    #total_cost=total_cost+current_cost
    if status == 'Infeasible_Problem_Detected':
        count=count+1
        
        pass
    #np.array(u)
    A=mod_aug.integrator_jac_x_a(X_true,u)
    u_D_plot=np.append(u_D_plot,u[0])
    u_delta_plot=np.append(u_delta_plot,u[1])

    plant_model=plant_dyn_vary#options plan_dyn_vary and plant_dyn_const
    observer_model=observer_dyn_offline#options observer_dyn and observer_dyn_offline


    #plant x
    X_true,P_model_cov_true=plant_model(X_true,u,A,P_model_cov_est,i)
    #observer   
    X_est,P_model_cov_est,Y_meas=observer_model(P_model_cov_true,X_est,u)
    print("xtrue and xest values are")
    print(X_true)
    print(X_est)
    #plot plant and controller
    s_plot_meas_y=np.append(s_plot_meas_y,Y_meas[0])
    n_plot_meas_y=np.append(n_plot_meas_y,Y_meas[1])
    alpha_plot_meas_y=np.append(alpha_plot_meas_y,Y_meas[2])

    s_plot=np.append(s_plot,X_true[0])
    n_plot=np.append(n_plot,X_true[1])
    alpha_plot=np.append(alpha_plot,X_true[2])
    velocity_plot=np.append(velocity_plot,X_true[3])
    model_error_plot_true_1=np.append(model_error_plot_true_1,X_true[4])
    model_error_plot_true_2=np.append(model_error_plot_true_2,X_true[5])
    model_error_plot_true_3=np.append(model_error_plot_true_3,X_true[6])
    
    s_plot_est=np.append(s_plot_est,X_est[0])
    n_plot_est=np.append(n_plot_est,X_est[1])
    alpha_plot_est=np.append(alpha_plot_est,X_est[2])
    velocity_plot_est=np.append(velocity_plot_est,X_est[3])
    model_error_plot_est_1=np.append(model_error_plot_est_1,X_est[4])
    model_error_plot_est_2=np.append(model_error_plot_est_2,X_est[5])
    model_error_plot_est_3=np.append(model_error_plot_est_3,X_est[6])
    #print("value of x_reference is")
    #print(x_reference)    

print("count value is")
print(count)    
#toc=time.perf_counter() 

#np.savetxt('configurations/x_est_cost_function',cost_list)
#np.savetxt('configurations/x_cost_function',cost_list)



#np.savetxt("configurations/augmented_0.4_mpc",np.transpose([model_error_plot_est_1,model_error_plot_est_2,model_error_plot_est_3]))

if(cost_flag=='non_disturbed'):
    np.savetxt('configurations/cost_function_mpc_module_ndisturbed',cost_list)
    print("cost_values of the non disturbed model mpc")
    print(cost_list)

    ax1 = plt.subplot2grid(shape=(4,4), loc=(0,0),rowspan=2 ,colspan=2)
    plt.ylabel('s')
    plt.xlabel('N')
    plt.plot(range(N_main+1),s_plot,'r',label='true')
    plt.plot(range(N_main+1),s_plot_est,'g',label='est')
    plt.plot(range(N_main),x_reference[:,0],'b',label='reference')
    plt.tight_layout()
    plt.legend()

    ax2 = plt.subplot2grid(shape=(4,4), loc=(0,2),rowspan=2 , colspan=2)
    plt.ylabel('n')
    plt.xlabel('N')
    plt.plot(range(N_main+1),n_plot,'r',label='true')
    plt.plot(range(N_main+1),n_plot_est,'g',label='est')
    plt.plot(range(N_main),x_reference[:,1],'b',label='reference')
    plt.tight_layout()
    plt.legend()

    ax3 = plt.subplot2grid(shape=(4,4), loc=(2,0),rowspan=2 , colspan=2)
    plt.ylabel(r'$\alpha$')
    plt.xlabel('N')
    plt.plot(range(N_main+1),alpha_plot,'r',label='true')
    plt.plot(range(N_main+1),alpha_plot_est,'g',label='est')
    plt.plot(range(N_main),x_reference[:,2],'b',label='reference')
    plt.tight_layout()
    plt.legend()

    ax4 = plt.subplot2grid(shape=(4,4), loc=(2,2),rowspan=2 , colspan=2)
    plt.ylabel('v')
    plt.xlabel('N')
    plt.plot(range(N_main+1),velocity_plot,'r',label='true')
    plt.plot(range(N_main+1),velocity_plot_est,'g',label='est')
    plt.plot(range(N_main),x_reference[:,3],'b',label='reference') 
    plt.tight_layout()
    plt.legend()
    #plt.savefig('plots_save/true_and_est_plots.png')
    #plt.savefig('plots_save/true_and_est_plots estimated disturbances.png')
    
    plt.savefig('plots_save/non_disturbed/mpc_ekf_results/true_and_est_plots.pgf')
    plt.savefig('plots_save/non_disturbed/mpc_ekf_results/true_and_est_plots.pdf',bbox_inches='tight')
    plt.show()
    

    print(model_error_plot_true_1)
    """ d_cov=0.2*np.diag([1,1,1,1])


    for i in range(N_main):
        x_reference[[i],:]=np.random.multivariate_normal([0.1,0.1,0.1,0.1],d_cov,1)
    """
    x_track,y_track,psi_track,v_track=SpattoOrig(s_plot_est,n_plot_est,alpha_plot_est,velocity_plot_est,s0,xref,yref,psiref)
    x_track_true,y_track_true,psi_track_true,v_track_true=SpattoOrig(s_plot,n_plot,alpha_plot,velocity_plot,s0,xref,yref,psiref)
    plt.figure()
    plt.xlabel('x')
    plt.ylabel('y')
    plt.plot(x_track,y_track,label='Estimated output of the  car')
    plt.plot(x_track_true,y_track_true,label='Plant output of the  car')
    #plt.plot(x_track_1,y_track_1)
    plt.plot(xref,yref,label='track')
    plt.tight_layout()
    plt.legend()
    plt.savefig('plots_save/non_disturbed/track/track_performance.pdf',bbox_inches='tight')
    plt.savefig('plots_save/non_disturbed/track/track_performance.pgf')
    
    plt.show()

    plt.figure()
    ax1 = plt.subplot2grid(shape=(4,4), loc=(0,0), rowspan=2,colspan=2)
    # ax1.set(xlim=(0,N), ylim=(-0.5, 0.5))
    
    plt.ylabel(r'$d_s$')
    plt.xlabel('N')
    plt.plot(range(N_main),model_error_plot_true_1,'r',label='true_model_error_s')
    plt.plot(range(N_main),model_error_plot_est_1,'g',label='estimated_model_error_s')
    plt.tight_layout()
    plt.legend()


    ax2 = plt.subplot2grid((4,4), (0,2), rowspan=2,colspan=2)
    plt.ylabel(r'$d_n$')
    plt.xlabel('N')
    plt.plot(range(N_main),model_error_plot_true_2,'r',label='true')
    plt.plot(range(N_main),model_error_plot_est_2,'g',label='est')
    plt.tight_layout()
    plt.legend()

    ax3 = plt.subplot2grid((4,4), (2,1), rowspan=2,colspan=2)
    plt.ylabel(r'$d_\alpha$')
    plt.xlabel('N')
    plt.plot(range(N_main),model_error_plot_true_3,'r',label='true')
    plt.plot(range(N_main),model_error_plot_est_3,'g',label='est')
    plt.tight_layout()
    plt.legend()
    plt.show()

    plt.figure()
    ax1 = plt.subplot2grid(shape=(4,4), loc=(0,0), rowspan=2,colspan=2) 
    plt.xlabel('s reference')
    plt.ylabel(r'$d_s$')
    plt.plot(s0[0:N_main],model_error_plot_est_1,'g',label='est',)
    plt.plot(s0[0:N_main],model_error_plot_true_1,'r',label='true')
    plt.tight_layout()
    plt.legend()

    ax2 = plt.subplot2grid(shape=(4,4), loc=(0,2), rowspan=2,colspan=2) 
    plt.xlabel('s reference')
    plt.ylabel(r'$d_n$')
    plt.plot(s0[0:N_main],model_error_plot_est_2,'g',label='est')
    plt.plot(s0[0:N_main],model_error_plot_true_2,'r',label='true')
    plt.tight_layout()
    plt.legend()

    ax3 = plt.subplot2grid(shape=(4,4), loc=(2,1), rowspan=2,colspan=2)
    plt.xlabel('s reference')
    plt.ylabel(r'$d_\alpha$')
    plt.plot(s0[0:N_main],model_error_plot_est_3,'g',label='est')
    plt.plot(s0[0:N_main],model_error_plot_true_3,'r',label='true')
    plt.tight_layout()
    plt.legend()
    plt.savefig('plots_save/non_disturbed/disturbances_sref/disturbance vs s-rederence.pgf')
    plt.savefig('plots_save/non_disturbed/disturbances_sref/disturbance vs s-rederence.pdf',bbox_inches='tight')
    
    plt.show()

    
elif(cost_flag=='disturbed'):
    if (observer_model==observer_dyn):    
        np.savetxt('configurations/cost_function_mpc_module_disturbed',cost_list)
        #np.savetxt('configurations/cost_function_mpc_module_disturbed_offline',cost_list)
        print("cost_values of the disturbed model mpc")
        print(cost_list)

        ax1 = plt.subplot2grid(shape=(4,4), loc=(0,0), rowspan=2,colspan=2)
        plt.ylabel('s')
        plt.xlabel('N')
        plt.plot(range(N_main+1),s_plot,'r',label='true')
        plt.plot(range(N_main+1),s_plot_est,'g',label='est')
        plt.plot(range(N_main),x_reference[:,0],'b',label='reference')
        plt.tight_layout()
        plt.legend()

        ax2 = plt.subplot2grid(shape=(4,4), loc=(0,2), rowspan=2,colspan=2)
        plt.ylabel('n')
        plt.xlabel('N')
        plt.plot(range(N_main+1),n_plot,'r',label='true')
        plt.plot(range(N_main+1),n_plot_est,'g',label='est')
        plt.plot(range(N_main),x_reference[:,1],'b',label='reference')
        plt.tight_layout()
        plt.legend()

        ax3 = plt.subplot2grid(shape=(4,4), loc=(2,0), rowspan=2,colspan=2)
        plt.ylabel(r'$\alpha$')
        plt.xlabel('N')
        plt.plot(range(N_main+1),alpha_plot,'r',label='true')
        plt.plot(range(N_main+1),alpha_plot_est,'g',label='est')
        plt.plot(range(N_main),x_reference[:,2],'b',label='reference')
        plt.tight_layout()
        plt.legend()

        ax4 = plt.subplot2grid(shape=(4,4), loc=(2,2), rowspan=2,colspan=2)
        plt.ylabel('v')
        plt.xlabel('N')
        plt.plot(range(N_main+1),velocity_plot,'r',label='true')
        plt.plot(range(N_main+1),velocity_plot_est,'g',label='est')
        plt.plot(range(N_main),x_reference[:,3],'b',label='reference') 
        plt.tight_layout()
        plt.legend()

        plt.savefig('plots_save/disturbed/mpc_ekf_results/true_and_est_plots.pgf')
        plt.savefig('plots_save/disturbed/mpc_ekf_results/true_and_est_plots.pdf')
        #plt.savefig('plots_save/true_and_est_plots estimated disturbances.png')
        
        plt.show()
        #x_track,y_track,psi_track,v_track=SpattoOrig(x_reference[:,0],x_reference[:,1],x_reference[:,2],x_reference[:,3],s0,xref,yref,psiref)
        
        
        
        x_track,y_track,psi_track,v_track=SpattoOrig(s_plot_est,n_plot_est,alpha_plot_est,velocity_plot_est,s0,xref,yref,psiref)
        x_track_true,y_track_true,psi_track_true,v_track_true=SpattoOrig(s_plot,n_plot,alpha_plot,velocity_plot,s0,xref,yref,psiref)
        plt.figure()
        plt.xlabel('x')
        plt.ylabel('y')
        plt.plot(x_track,y_track,label='Estimated output of the  car')
        plt.plot(x_track_true,y_track_true,label='Plant output of the  car')
        #plt.plot(x_track_1,y_track_1)
        plt.plot(xref,yref,label='track')
        plt.legend()
        plt.savefig('plots_save/disturbed/track/track_performance.pgf')
        plt.savefig('plots_save/disturbed/track/track_performance.pdf')
        #plt.savefig('plots_save/track_performance estimated disturbances.png')
        #plt.plot(xref,yref)
        
        plt.show()

        plt.figure()
        ax1 = plt.subplot2grid(shape=(4,4), loc=(0,0), rowspan=2,colspan=2)
        # ax1.set(xlim=(0,N), ylim=(-0.5, 0.5))
        plt.ylabel(r'$d_\s')
        plt.xlabel('N')
        plt.plot(range(N_main),model_error_plot_true_1,'r',label='true')
        plt.plot(range(N_main),model_error_plot_est_1,'g',label='est')
        plt.tight_layout()
        plt.legend()


        ax2 = plt.subplot2grid((4,4), (0,2), rowspan=2,colspan=2)
        plt.ylabel(r'$d_n$')
        plt.xlabel('N')
        plt.plot(range(N_main),model_error_plot_true_2,'r',label='true')
        plt.plot(range(N_main),model_error_plot_est_2,'g',label='est')
        plt.tight_layout()
        plt.legend()

        ax3 = plt.subplot2grid((4,4), (2,1), rowspan=2,colspan=2)
        plt.ylabel(r'$d_\alpha$')
        plt.xlabel('N')
        plt.plot(range(N_main),model_error_plot_true_3,'r',label='true')
        plt.plot(range(N_main),model_error_plot_est_3,'g',label='est')
        plt.tight_layout()
        plt.legend()
        plt.show()

        plt.figure()
        ax1 = plt.subplot2grid(shape=(4,4), loc=(0,0), rowspan=2,colspan=2) 
        plt.xlabel('s reference')
        plt.ylabel(r'$d_s$')
        plt.plot(s0[0:N_main],model_error_plot_est_1,'g',label='est',)
        plt.plot(s0[0:N_main],model_error_plot_true_1,'r',label='true')
        plt.tight_layout()
        plt.legend()

        ax2 = plt.subplot2grid(shape=(4,4), loc=(0,2), rowspan=2,colspan=2) 
        plt.xlabel('s reference')
        plt.ylabel(r'$d_n$')
        plt.plot(s0[0:N_main],model_error_plot_est_2,'g',label='est')
        plt.plot(s0[0:N_main],model_error_plot_true_2,'r',label='true')
        plt.tight_layout()
        plt.legend()

        ax3 = plt.subplot2grid(shape=(4,4), loc=(2,1), rowspan=2,colspan=2)
        plt.xlabel('s reference')
        plt.ylabel(r'$d_\alpha$')
        plt.plot(s0[0:N_main],model_error_plot_est_3,'g',label='est')
        plt.plot(s0[0:N_main],model_error_plot_true_3,'r',label='true')
        plt.tight_layout()
        plt.legend()
        
        plt.savefig('plots_save/disturbed/disturbances_sref/disturbances-sref.pdf')
        plt.savefig('plots_save/disturbed/disturbances_sref/disturbances-sref.pgf')
        plt.show()
        print("going not in offline mode")

    elif(observer_model==observer_dyn_offline):
        np.savetxt('configurations/cost_function_mpc_module_disturbed_offline',cost_list)
        print("cost_values of the disturbed_offline model mpc")
        print(cost_list)

        ax1 = plt.subplot2grid(shape=(4,4), loc=(0,0), rowspan=2,colspan=2)
        plt.ylabel('s')
        plt.xlabel('N')
        plt.plot(range(N_main+1),s_plot,'r',label='true')
        plt.plot(range(N_main+1),s_plot_est,'g',label='est')
        plt.plot(range(N_main),x_reference[:,0],'b',label='reference')
        plt.tight_layout()
        plt.legend()

        ax2 = plt.subplot2grid(shape=(4,4), loc=(0,2), rowspan=2,colspan=2)
        plt.ylabel('n')
        plt.xlabel('N')
        plt.plot(range(N_main+1),n_plot,'r',label='true')
        plt.plot(range(N_main+1),n_plot_est,'g',label='est')
        plt.plot(range(N_main),x_reference[:,1],'b',label='reference')
        plt.tight_layout()
        plt.legend()

        ax3 = plt.subplot2grid(shape=(4,4), loc=(2,0), rowspan=2,colspan=2)
        plt.ylabel(r'$\alpha$')
        plt.xlabel('N')
        plt.plot(range(N_main+1),alpha_plot,'r',label='true')
        plt.plot(range(N_main+1),alpha_plot_est,'g',label='est')
        plt.plot(range(N_main),x_reference[:,2],'b',label='reference')
        plt.tight_layout()
        plt.legend()

        ax4 = plt.subplot2grid(shape=(4,4), loc=(2,2), rowspan=2,colspan=2)
        plt.ylabel('v')
        plt.xlabel('N')
        plt.plot(range(N_main+1),velocity_plot,'r',label='true')
        plt.plot(range(N_main+1),velocity_plot_est,'g',label='est')
        plt.plot(range(N_main),x_reference[:,3],'b',label='reference') 
        plt.tight_layout()
        plt.legend()

        #plt.savefig('plots_save/disturbed/offline/mpc_ekf_results/true_and_est_plots.pgf')
        plt.savefig('plots_save/disturbed/offline/mpc_ekf_results/true_and_est_plots.pdf')
        #plt.savefig('plots_save/true_and_est_plots estimated disturbances.png')
        
        plt.show()
        #x_track,y_track,psi_track,v_track=SpattoOrig(x_reference[:,0],x_reference[:,1],x_reference[:,2],x_reference[:,3],s0,xref,yref,psiref)
        x_track,y_track,psi_track,v_track=SpattoOrig(s_plot_est,n_plot_est,alpha_plot_est,velocity_plot_est,s0,xref,yref,psiref)
        x_track_true,y_track_true,psi_track_true,v_track_true=SpattoOrig(s_plot,n_plot,alpha_plot,velocity_plot,s0,xref,yref,psiref)
        plt.figure()
        plt.xlabel('x')
        plt.ylabel('y')
        plt.plot(x_track,y_track,label='Estimated output of the  car')
        plt.plot(x_track_true,y_track_true,label='Plant output of the  car')
        #plt.plot(x_track_1,y_track_1)
        plt.plot(xref,yref,label='track')
        plt.savefig('plots_save/disturbed/offline/track/track_performance.pgf')
        plt.savefig('plots_save/disturbed/offline/track/track_performance.pdf')
        #plt.savefig('plots_save/track_performance estimated disturbances.png')
        #plt.plot(xref,yref)
        plt.show()

        plt.figure()
        ax1 = plt.subplot2grid(shape=(4,4), loc=(0,0), rowspan=2,colspan=2)
        # ax1.set(xlim=(0,N), ylim=(-0.5, 0.5))
        plt.ylabel(r'$d_\s')
        plt.xlabel('N')
        plt.plot(range(N_main),model_error_plot_true_1,'r',label='true')
        plt.plot(range(N_main),model_error_plot_est_1,'g',label='est')
        plt.tight_layout()
        plt.legend()


        ax2 = plt.subplot2grid((4,4), (0,2), rowspan=2,colspan=2)
        plt.ylabel(r'$d_n$')
        plt.xlabel('N')
        plt.plot(range(N_main),model_error_plot_true_2,'r',label='true')
        plt.plot(range(N_main),model_error_plot_est_2,'g',label='est')
        plt.tight_layout()
        plt.legend()

        ax3 = plt.subplot2grid((4,4), (2,1), rowspan=2,colspan=2)
        plt.ylabel(r'$d_\alpha$')
        plt.xlabel('N')
        plt.plot(range(N_main),model_error_plot_true_3,'r',label='true')
        plt.plot(range(N_main),model_error_plot_est_3,'g',label='est')
        plt.tight_layout()
        plt.legend()
        plt.show()

        plt.figure()
        ax1 = plt.subplot2grid(shape=(4,4), loc=(0,0), rowspan=2,colspan=2) 
        plt.xlabel('s reference')
        plt.ylabel('s')
        plt.plot(s0[0:N_main],model_error_plot_est_1,'g',label='s estimated_disturbance',)
        plt.plot(s0[0:N_main],model_error_plot_true_1,'r',label='s_true_disturbance')
        plt.tight_layout()
        plt.legend()

        ax2 = plt.subplot2grid(shape=(4,4), loc=(0,2), rowspan=2,colspan=2) 
        plt.xlabel('s reference')
        plt.ylabel('n')
        plt.plot(s0[0:N_main],model_error_plot_est_2,'g',label='est')
        plt.plot(s0[0:N_main],model_error_plot_true_2,'r',label='true')
        plt.tight_layout()
        plt.legend()

        ax3 = plt.subplot2grid(shape=(4,4), loc=(2,1), rowspan=2,colspan=2)
        plt.xlabel('s reference')
        plt.ylabel(r'$d_\alpha$')
        plt.plot(s0[0:N_main],model_error_plot_est_3,'g',label='est')
        plt.plot(s0[0:N_main],model_error_plot_true_3,'r',label='true')
        plt.tight_layout()
        plt.legend()
        
        plt.savefig('plots_save/disturbed/offline/disturbances_sref/disturbances-sref.pdf')
        plt.savefig('plots_save/disturbed/offline/disturbances_sref/disturbances-sref.pgf')
        plt.show()    
        print("going in offline mode")

""" 
np.savetxt("configurations/est-true_mpc",np.transpose([s_plot,n_plot,alpha_plot,velocity_plot,\
    s_plot_est,n_plot_est,alpha_plot_est,velocity_plot_est,\
]))
np.savetxt("configurations/model-error-est-true",np.transpose([model_error_plot_est_1,model_error_plot_est_2,model_error_plot_est_3,\
    model_error_plot_true_1,model_error_plot_true_2,model_error_plot_true_3]))
np.savetxt("configurations/u_data",np.transpose([u_D_plot,u_delta_plot])) 

np.savetxt("configurations/y_measured_vk",np.transpose([s_plot_meas_y,n_plot_meas_y,alpha_plot_meas_y]))



 """
