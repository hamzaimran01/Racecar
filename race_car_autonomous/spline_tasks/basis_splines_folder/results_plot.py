import sys
sys.path.append("../THESIS")
import numpy as np
from casadi import *
import matplotlib.pyplot as plt
#from spline_tasks.custom_splines_continous import N_MHE
from spline_tasks.plots import plotTrackProj
from spline_tasks.model.tracks.readDataFcn import getTrack

#N_MHE=593
N_MHE=1610
#N_MHE=1565
#X_measurement = np.loadtxt("final_data/0.8_training_data/states_final_alligned_0.8_training")
#X_estimated = np.loadtxt("final_data/0.8_training_data/estimated_data/simX_res.txt")


""" 
if want to change for different configuration change only 
->X_measurement,X_estimated
->savefig


 """


def plots_individual_states(N_MHE,X_measurement,X_estimated):

    #X_measurement = np.loadtxt("final_data/0.8_training_data/states_final_alligned_0.8_training")
    #X_estimated = np.loadtxt("final_data/0.8_training_data/estimated_data/simX_res.txt")


    #X_measurement = np.loadtxt("final_data/1.4_data/states_final_alligned_1.4")
    #X_estimated = np.loadtxt("final_data/1.4_data/estimated_data/simX_res.txt")
    #X_measurement = np.loadtxt("final_data/0.8_training_data/states_final_alligned_0.8_training")
    #X_estimated = np.loadtxt("final_data/0.8_training_data/estimated_data/simX_res.txt")

    result_s = X_estimated[0:N_MHE,0]-X_measurement[0:N_MHE,0]
    result_n = X_estimated[0:N_MHE,1]-X_measurement[0:N_MHE,1]
    result_alpha = X_estimated[0:N_MHE,2]-X_measurement[0:N_MHE,2]
    result_v = X_estimated[0:N_MHE,3]-X_measurement[0:N_MHE,3]
    s_plot=np.linspace(4.695507388859999909e-02,8.6502137881266751407e+00,N_MHE)

    fig,axs=plt.subplots(2,2)

    axs[0,0].set_xlabel(r'$s_{ref}[m]$ track')
    axs[0,0].set_ylabel(r'$s[m]$')
    #axs[0,0].set(xlabel=r'$s_{ref}[m]$ track', ylabel=r'$s[m]$')
    axs[0,0].bbox
    """ axs[0,0].vlines(x=1.38e+00,ymin=0,ymax=8.9,color='g',linestyles='dotted',label='reset 1',linewidth=2)
    axs[0,0].vlines(x=2.88e+00,ymin=0,ymax=8.9,color='r',linestyles='dotted',label='reset 2',linewidth=2)
    axs[0,0].vlines(x=3.64e+00,ymin=0,ymax=8.9,color='k',linestyles='dotted',label='reset 3',linewidth=2)
    axs[0,0].vlines(x=4.52e+00,ymin=0,ymax=8.9,color='m',linestyles='dotted',label='reset 4',linewidth=2)
    axs[0,0].vlines(x=6.8534e+00,ymin=0,ymax=8.9,color='r',linestyles='dotted',label='reset 5',linewidth=2) """
    
    #axs[0,0].vlines(x=6.866119772930000309e+00,ymin=0,ymax=8,color='r',linestyles='dotted',label='reset 5')
    
    axs[0,0].plot(s_plot,X_measurement[0:N_MHE,0],'k',label=r'$s$ measured',linewidth=4)
    axs[0,0].plot(s_plot,X_estimated[0:N_MHE,0],'y',alpha=0.7,label=r'$s$ with model mismatch ',linewidth=4)
    s_difference = X_measurement[0:N_MHE,0] - X_estimated[0:N_MHE,0]
    axs[0,0].plot(s_plot,s_difference,'b',label=r'difference ')
    axs[0,0].legend(ncol=2)





    axs[0,1].set_xlabel(r'$s_{ref}[m]$ track')
    axs[0,1].set_ylabel(r'$n[m]$')
    axs[0,1].bbox
    """ axs[0,1].vlines(x=1.38e+00,ymin=-0.20,ymax=0.20,color='g',linestyles='dotted',label='reset 1',linewidth=2)
    axs[0,1].vlines(x=2.88e+00,ymin=-0.20,ymax=0.20,color='r',linestyles='dotted',label='reset 2',linewidth=2)
    axs[0,1].vlines(x=3.64e+00,ymin=-0.20,ymax=0.20,color='k',linestyles='dotted',label='reset 3',linewidth=2)
    axs[0,1].vlines(x=4.52e+00,ymin=-0.20,ymax=0.20,color='m',linestyles='dotted',label='reset 4',linewidth=2)
    axs[0,1].vlines(x=6.8534e+00,ymin=-0.20,ymax=0.20,color='r',linestyles='dotted',label='reset 5',linewidth=2) """
    axs[0,1].plot(s_plot,X_measurement[0:N_MHE,1],'k',label=r'$n$ measurement',linewidth=4)
    axs[0,1].plot(s_plot,X_estimated[0:N_MHE,1],'y',alpha=0.7,label=r'$n$ with model mismatch',linewidth=4)
    
    n_difference = X_measurement[0:N_MHE,1] - X_estimated[0:N_MHE,1]
    axs[0,1].plot(s_plot,n_difference,'b',label=r'difference')
    axs[0,1].legend(ncol=2)


    axs[1,0].set_xlabel(r'$s_{ref}[m]$ track')
    axs[1,0].set_ylabel(r'$\alpha[\degree]$')
    axs[1,0].bbox
    """ axs[1,0].vlines(x=1.38e+00,ymin=-1,ymax=7,color='g',linestyles='dotted',label='reset 1',linewidth=2)
    axs[1,0].vlines(x=2.88e+00,ymin=-1,ymax=7,color='r',linestyles='dotted',label='reset 2',linewidth=2)
    axs[1,0].vlines(x=3.64e+00,ymin=-1,ymax=7,color='k',linestyles='dotted',label='reset 3',linewidth=2)
    axs[1,0].vlines(x=4.52e+00,ymin=-1,ymax=7,color='m',linestyles='dotted',label='reset 4',linewidth=2)
    axs[1,0].vlines(x=6.8534e+00,ymin=-1,ymax=1.5,color='r',linestyles='dotted',label='reset 5',linewidth=2) """

    axs[1,0].plot(s_plot,X_measurement[0:N_MHE,2],'k',label=r'$\alpha$ measurement',linewidth=4)
    axs[1,0].plot(s_plot,X_estimated[0:N_MHE,2],'y',label=r'$\alpha$ with model mismatch',linewidth=4)
    
    
    alpha_difference = X_measurement[0:N_MHE,2] - X_estimated[0:N_MHE,2]
    #np.savetxt("final_data/0.8_training_data/estimated_data/alpha_difference.txt",alpha_difference)
    #np.savetxt("final_data/0.8_test_data/estimated_data/alpha_difference.txt",alpha_difference)
    axs[1,0].plot(s_plot,alpha_difference,'b',label=r'difference')
    axs[1,0].legend(ncol=2)

    #plt.plot(s_plot,result_alpha,label='alpha_difference')
    
    
    
    
    
    
    axs[1,0].legend(ncol=2)


    axs[1,1].set_xlabel(r'$s_{ref}[m]$ track')
    axs[1,1].set_ylabel(r'$v$$[\frac{m}{s}]$')
    axs[1,1].bbox
    """ axs[1,1].vlines(x=1.38e+00,ymin=-0.2,ymax=0.7,color='g',linestyles='dotted',label='reset 1',linewidth=2)
    axs[1,1].vlines(x=2.88e+00,ymin=-0.2,ymax=0.7,color='r',linestyles='dotted',label='reset 2',linewidth=2)
    axs[1,1].vlines(x=3.64e+00,ymin=-0.2,ymax=0.7,color='k',linestyles='dotted',label='reset 3',linewidth=2)
    axs[1,1].vlines(x=4.52e+00,ymin=-0.2,ymax=0.7,color='m',linestyles='dotted',label='reset 4',linewidth=2)
    axs[1,1].vlines(x=6.8534e+00,ymin=0,ymax=0.5,color='r',linestyles='dotted',label='reset 5',linewidth=2) """
    
    axs[1,1].plot(s_plot,X_measurement[0:N_MHE,3],'k',label=r'$v$ measurement',linewidth=4)
    axs[1,1].plot(s_plot,X_estimated[0:N_MHE,3],'y',label=r'$v$ with model mismatch',linewidth=4)
    
    v_difference = X_measurement[0:N_MHE,3] - X_estimated[0:N_MHE,3]
    #np.savetxt("final_data/0.8_training_data/estimated_data/v_difference.txt",v_difference)
    #np.savetxt("final_data/0.8_test_data/estimated_data/v_difference.txt",v_difference)
    axs[1,1].plot(s_plot,v_difference,'b',label=r'difference')
    
    
    
    axs[1,1].legend(ncol=2)

    #plt.tight_layout()
    #plt.savefig('final_data/0.8_training_data/estimated_data/plots_training/0.8_trainings_n_alpha_measurement_estimated.pdf')
    

    figure = plt.gcf()  # get current figure
    figure.set_size_inches(24, 18) # set figure's size manually to your full screen (32x18)
    

    #plt.savefig('final_data/0.8_training_data/estimated_data/plots_training/0.8_trainings_n_alpha_measurement_estimated.pdf',bbox_inches='tight')
    plt.show() 

    track = "LMS_Track.txt"
    [s0,xref,yref,psiref,kapparef] = getTrack(track)
    kapparef_s=interpolant('kapparef_s','bspline',[s0],kapparef)

    plt.show()


""" 

    fig,axs=plt.subplots(2,2)

    axs[0,0].set_xlabel(r'$s_{ref}[m]$ track')
    axs[0,0].set_ylabel(r'$s[m]$')
    #axs[0,0].set(xlabel=r'$s_{ref}[m]$ track', ylabel=r'$s[m]$')
    axs[0,0].bbox
    s_difference = X_measurement[0:N_MHE,0] - X_estimated[0:N_MHE,0]
    #np.savetxt("final_data/0.8_training_data/estimated_data/s_difference.txt",s_difference)
    #np.savetxt("final_data/0.8_test_data/estimated_data/s_difference.txt",s_difference)
    axs[0,0].plot(s_plot,s_difference,'g',label=r'difference: $s$ measurement and $s$ with model mismatch ')
    axs[0,0].legend()


   

    axs[0,1].set_xlabel(r'$s_{ref}[m]$ track')
    axs[0,1].set_ylabel(r'$n$[m]')
    #axs[0,0].set(xlabel=r'$s_{ref}[m]$ track', ylabel=r'$s[m]$')
    axs[0,1].bbox
    n_difference = X_measurement[0:N_MHE,1] - X_estimated[0:N_MHE,1]
    #np.savetxt("final_data/0.8_training_data/estimated_data/n_difference.txt",n_difference)
    np.savetxt("final_data/0.8_test_data/estimated_data/n_difference.txt",n_difference)
    axs[0,1].plot(s_plot,n_difference,'g',label=r'difference: $n$ measurement and $n$ with model mismatch ')
    axs[0,1].legend()




    axs[1,0].set_xlabel(r'$s_{ref}[m]$ track')
    axs[1,0].set_ylabel(r'$\alpha$[m]')
    #axs[0,0].set(xlabel=r'$s_{ref}[m]$ track', ylabel=r'$s[m]$')
    axs[1,0].bbox
    alpha_difference = X_measurement[0:N_MHE,2] - X_estimated[0:N_MHE,2]
    #np.savetxt("final_data/0.8_training_data/estimated_data/alpha_difference.txt",alpha_difference)
    #np.savetxt("final_data/0.8_test_data/estimated_data/alpha_difference.txt",alpha_difference)
    axs[1,0].plot(s_plot,alpha_difference,'g',label=r'difference: $\alpha$ measurement and $\alpha$ with model mismatch ')
    axs[1,0].legend()




    axs[1,1].set_xlabel(r'$s_{ref}[m]$ track')
    axs[1,1].set_ylabel(r'$v$[m]')
    #axs[0,0].set(xlabel=r'$s_{ref}[m]$ track', ylabel=r'$s[m]$')
    axs[1,1].bbox
    v_difference = X_measurement[0:N_MHE,3] - X_estimated[0:N_MHE,3]
    #np.savetxt("final_data/0.8_training_data/estimated_data/v_difference.txt",v_difference)
    #np.savetxt("final_data/0.8_test_data/estimated_data/v_difference.txt",v_difference)
    axs[1,1].plot(s_plot,v_difference,'g',label=r'difference: $v$ measurement and $v$ with model mismatch ')
    axs[1,1].legend()

    #plt.tight_layout()
    
    
    #fig.set_size_inches((11, 9), forward=False)
    #figure = plt.gcf()  # get current figure
    #figure.set_size_inches(32, 18) # set figure's size manually to your full screen (32x18)
    
    figure = plt.gcf()  # get current figure
    figure.set_size_inches(24, 18) # set figure's size manually to your full screen (32x18)
    fig.savefig('final_data/0.8_test_data/estimated_data/plots_training/1.1_test_difference_s_n_alpha_v.pdf',bbox_inches='tight')

 """

    #plt.show()

    
#N_MHE=N_MHE=757
N_MHE=756
plt.rcParams.update({'font.size': 16})
x_measured=X_measurement=np.loadtxt("final_data/1.1_test_data/states_final_alligned_test_1.1")
simX_res=np.loadtxt('final_data/1.1_test_data/estimated_data/simX_res.txt')

#plots_individual_states(N_MHE,X_measurement,simX_res)

#plt.show()

print(x_measured.shape)
segments_1_simX_res=np.zeros((118,6))
residual_1=np.zeros((118,6))
segments_2_simX_res=np.zeros((250-118,6))
residual_2=np.zeros((250-118,6))
segments_3_simX_res=np.zeros((318-250,6))
residual_3=np.zeros((318-250,6))
segments_4_simX_res=np.zeros((394-318,6))
residual_4=np.zeros((394-318,6))
segments_5_simX_res=np.zeros((600-394,6))
residual_5=np.zeros((600-394,6))

segments_6_simX_res=np.zeros((756-600,6))
residual_6=np.zeros((756-600,6))
#segments_5_simX_res=np.zeros((756-600,6))
""" 
0.8_test_data 
1st segment=515


 """

for i in range(0,118):     

    if(simX_res[i,1]<=0.18 and simX_res[i,1]>=-0.18):
        segments_1_simX_res = simX_res[i,:]
        #segments_1_simX_res = simX_res[i,:]
        residual_1[i-0,0:6] = x_measured[i,0:6] - simX_res[i,:]
    else: 
        break

#[ 2.474916   -0.16148472 -1.30795378  0.28295276  0.46048051 -0.36909903]
for i in range(118,250):  
    if(simX_res[i,1]<=0.18 and simX_res[i,1]>=-0.18):
        segments_2_simX_res = simX_res[i,:]
        residual_2[i-118,0:6]=x_measured[i,0:6] - simX_res[i,:]
        #print("residual vlaue is")
       
    else: 
        break

    
for i in range(250,318):  

    if(simX_res[i,1]<=0.18 and simX_res[i,1]>=-0.18):
        segments_3_simX_res = simX_res[i,:]
        residual_3[i-250,0:6]=x_measured[i,0:6] - simX_res[i,:]
    else: 
        break


for i in range(318,394):  
    
    if(simX_res[i,1]<=0.18 and simX_res[i,1]>=-0.18):
        segments_4_simX_res = simX_res[i,:]
        #segments_4_simX_res = simX_res[i,:]
        residual_4[i-318,0:6]=x_measured[i,0:6] - simX_res[i,:]
    else: 
        break


for i in range(394,600):  
    
    if(simX_res[i,1]<=0.18 and simX_res[i,1]>=-0.18):
        segments_5_simX_res = simX_res[i,:]
        segments_5_simX_res = simX_res[i,:]
        residual_5[i-394,0:6]=x_measured[i,0:6] - simX_res[i,:]
        #print(simX_res[i,0:6])
        #residual_5[i-1129,0:6] = x_measured[i,0:6] - simX_res[i,:]
        
    else: 
        break

for i in range(600,756):  
    
    if(simX_res[i,1]<=0.18 and simX_res[i,1]>=-0.18):
        segments_6_simX_res = simX_res[i,:]
        segments_6_simX_res = simX_res[i,:]
        residual_6[i-600,0:6]=x_measured[i,0:6] - simX_res[i,:]
        #print(simX_res[i,0:6])
        #residual_5[i-1129,0:6] = x_measured[i,0:6] - simX_res[i,:]
        
    else: 
        break


""" 
for i in range(1129,1592):  
    
    if(simX_res[i,1]<=0.18 and simX_res[i,1]>=-0.18):
        segments_5_simX_res = simX_res[i,:]
        segments_5_simX_res = simX_res[i,:]
        #print(simX_res[i,0:6])
        #residual_5[i-1129,0:6] = x_measured[i,0:6] - simX_res[i,:]
        
    else: 
        break """

""" 
for i in range(1592,1680):  
    
    if(simX_res[i,1]<=0.18 and simX_res[i,1]>=-0.18):
        segments_6_simX_res = simX_res[i,:]
        segments_6_simX_res = simX_res[i,:]
        #print(simX_res[i,0:6])
        #residual_5[i-1129,0:6] = x_measured[i,0:6] - simX_res[i,:]
        
    else: 
        break

 """







#x_measured[i,0:6] - simX_res[i,:]
plt.plot(range(600,756), x_measured[600:756,3],label='measured_s')
plt.plot(range(600,756),simX_res[600:756,3],label='simulation_s')
plt.plot(range(600,756),x_measured[600:756,3]-simX_res[600:756,3],label='difference')
plt.legend()
plt.show()

s_diff=x_measured[600:756,3]-simX_res[600:756,3]
sum=0
s_mean=np.mean(s_diff,axis=0)
print("mean value is")
print(s_mean)


print("variance is")
cov_s = np.var(s_diff)
print(cov_s)

#s_diff=x_measured[1592:1680,0]-simX_res[1592:1680,0]
#print(x_measured[0:539,0])

