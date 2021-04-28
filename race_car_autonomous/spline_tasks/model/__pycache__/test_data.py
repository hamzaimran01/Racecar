import sys
sys.path.append("../THESIS")
import numpy as np
from casadi import *
#from spline_tasks.basis_splines_folder.spline_casadi_parameters import *
#from spline_tasks.optimization_initilizations import *
import matplotlib.pyplot as plt
from spline_tasks.casadi_integrator import casadi_integrator_fun as ci_f
from spline_tasks.readDataFcn import getTrack
#from plots import plotTrackProj
#from spline_tasks.basis_splines_folder.basis_splines import seg
from time import perf_counter 
#from spline_tasks.optimization_initilizations import *
from spline_tasks.model.tracks.readDataFcn import getTrack
from matplotlib import rc
from spline_tasks.plots import plotTrackProj
#from spline_tasks.basis_splines_folder.results_plot import plots_individual_states 
from time import perf_counter

#N_MHE=1560 #1.0
#N_MHE=1148#0.9
#N_MHE=1600#0.8
#N_MHE=1590
#N_MHE=1560#1.0
N_MHE=779#1.1
#N_MHE=1620

#X_measurement = np.loadtxt("final_data/0.8_training_data/states_final_alligned_0.8_training")
#X_measurement = np.loadtxt("final_data/1.0_data/states_final_alligned_1.0")
#X_estimated = np.loadtxt("final_data/1.0_data/estimated_data/simX_res.txt")
#con=np.loadtxt("final_data/0.8_training_data/control_final_alligned_0.8_training")
""" 
X_measurement = np.loadtxt("final_data/0.8_training_data/states_final_alligned_0.8_training")
#X_estimated = np.loadtxt("final_data/1.0_data/estimated_data/simX_res.txt")
con=np.loadtxt("final_data/0.8_training_data/control_final_alligned_0.8_training")

coefficeints_us=np.loadtxt('final_data/0.8_training_data/estimated_data/coeffcients_estimated.txt')
coefficients_s= np.zeros((4,258))


for i in range(200):
  #print("hi")  
  coefficients_s[:,i]  = coefficeints_us[4*i:4*i+4]
  print(coefficeints_us[4*i:4*i+4])  





basis_splines_table = np.loadtxt("spline_tasks/basis_splines")
t_segments_subpoints = np.loadtxt("spline_tasks/t_segments_subpoints")

nr_points=256
nr_sub_points=50

t=MX.sym('t',1)
t_1=MX.sym('t_1',1)
t_2=MX.sym('t_2',1)

f_1=(t-t_2)*(t-t_2)/2
f_2=0.5*(-2* (t-t_2)*(t-t_2) + 2*((t-t_2)) + 1)
f_3=0.5*(((t-t_2)*(t-t_2)) - 2*(t-t_2) + 1)
f_4=0
bspline_basis_1=Function('bspline_basis',[t,t_2],[f_1])
bspline_basis_2=Function('bspline_basis',[t,t_2],[f_2])
bspline_basis_3=Function('bspline_basis',[t,t_2],[f_3])
bspline_basis_4=Function('bspline_basis',[t,t_2],[f_4])
P_3=MX.sym('P_3',4,nr_points+2)
"""
#forming casadi function
""" 
f_1=(t-t_2)*(t-t_2)/2
f_2=0.5*(-2* (t-t_2)*(t-t_2) + 2*((t-t_2)) + 1)
f_3=0.5*(((t-t_2)*(t-t_2)) - 2*(t-t_2) + 1)
f_4=0
bspline_basis_1=Function('bspline_basis',[t,t_2],[f_1])
bspline_basis_2=Function('bspline_basis',[t,t_2],[f_2])
bspline_basis_3=Function('bspline_basis',[t,t_2],[f_3])
bspline_basis_4=Function('bspline_basis',[t,t_2],[f_4])


#points and subrange
#nr_points=64
nr_points=256
nr_sub_points=50

bspline_pieces= []
t_segments_subpoints=np.empty((0,50),int)   
x_segments_subpoints=np.empty((0,50),int)    

P_3=MX.sym('P',4,nr_points+2)


for i in range(nr_points):
    
    t_range=np.linspace(2 +i, 3 +i, 4)
    
    t_sub_range=np.linspace(2 +i,3 +i, nr_sub_points)
    t_sub_range=np.array([t_sub_range])
    #print(t_sub_range)
    #print(t_sub_range)
    x_1=[]
    x_2=[]
    x_3=[]
    
    t_segments_subpoints=np.append(t_segments_subpoints,t_sub_range,axis=0)
    
    
    
    for j in t_sub_range:
        #print(j)
        t_1=1*bspline_basis_1(j,2+i)
        t_1=t_1.full()
        t_2=1*bspline_basis_2(j,2+i)
        t_2=t_2.full()
        t_3=1*bspline_basis_3(j,2+i)
        t_3=t_3.full()
        
        x_1=np.append(x_1,t_1)    
        x_2=np.append(x_2,t_2)
        x_3=np.append(x_3,t_3)
        x_1=np.array([x_1])
        x_2=np.array([x_2])
        x_3=np.array([x_3])
        
            
    #appending the basis functions   
    for i in range(1):
        x_segments_subpoints=np.append(x_segments_subpoints,x_1,axis=0)
        x_segments_subpoints=np.append(x_segments_subpoints,x_2,axis=0)
        x_segments_subpoints=np.append(x_segments_subpoints,x_3,axis=0)

plt.xlabel('t')
plt.ylabel('$N_{i,j}(t)$')
for i in range(nr_points):
    plt.plot(t_segments_subpoints[i,:],x_segments_subpoints[3*i,:])
    plt.plot(t_segments_subpoints[i,:],x_segments_subpoints[3*i+1,:])
    plt.plot(t_segments_subpoints[i,:],x_segments_subpoints[3*i+2,:])

plt.legend    
#plt.savefig('complete_range_bspline.pdf')    
plt.show()

s=8.71

s_points=np.linspace(0,s,nr_points-1)
s_points=np.append(s_points,0)
s_segments=np.empty([nr_points-1,2])

t_points=np.linspace(2,nr_points,nr_points-1)
t_points=np.append(t_points,0)
t_segments=np.empty((nr_points-1,2))


for i in range(nr_points-1):
    s_segments[i,:]=s_points[0+i:2+i]
    t_segments[i,:]=t_points[0+i:2+i]
    





#for knowing the that which t_segment equals which s_segment
t_s_segments=np.concatenate((t_segments,s_segments),axis=1)

#for getting the s_val and finding it's corresponding t value such that we could form the bspline segments when needed
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx],idx

 
#uses the find_nearest function to tell us the corresponding t value
def bspline_fun_choose(t_s_segments,s_val):
   
    for i in range(nr_points):
        
        if(  t_s_segments[i,2] <= s_val <= t_s_segments[i,3] ):
            #subpoint_range=(t_s_segments[i,2]-t_s_segments[i,0])
            subrange_t=np.linspace(t_s_segments[i,0],t_s_segments[i,1],nr_sub_points)
            subrange_s=np.linspace(t_s_segments[i,2],t_s_segments[i,3],nr_sub_points)
            value_s,index=find_nearest(subrange_s,s_val)
            value_t=subrange_t[index]
            
            
            bspline_segment=t_s_segments[i,:]
            key=bspline_segment[0]
            return bspline_segment ,index,value_s,value_t,key



def bspline_choose(bspline_fun_choose, s_val, t_s_segments,t_segments_subpoints):
    for i in range(nr_points):
        
        bspline_segment, index, value, value_t, key = bspline_fun_choose(t_s_segments, s_val)
        
        #print(bspline_segment)

        t_range=np.linspace(2 +i, 3 +i, 4)
        
        t_sub_range=np.linspace(2 +i,3 +i, nr_sub_points)
        t_sub_range=np.array([t_sub_range])

        x_1=[]
        x_2=[]
        x_3=[]
        
        t_segments_subpoints=np.append(t_segments_subpoints,t_sub_range,axis=0)
        
        
        
        for j in t_sub_range:
            #print(j)
            t_1=1*bspline_basis_1(j,2+i)
            t_1=t_1.full()
            t_2=1*bspline_basis_2(j,2+i)
            t_2=t_2.full()
            t_3=1*bspline_basis_3(j,2+i)
            t_3=t_3.full()
            
            x_1=np.append(x_1,t_1)    
            x_2=np.append(x_2,t_2)
            x_3=np.append(x_3,t_3)
            x_1=np.array([x_1])
            x_2=np.array([x_2])
            x_3=np.array([x_3])
            
        
        if(2+i==key):
            seg_1=P_3[:,i] *bspline_basis_1(value_t,2+i)  + P_3[:,i+1] *bspline_basis_2(value_t,2+i) + P_3[:,i+2] *bspline_basis_3(value_t,2+i)
            #seg_1=P[i,:] *bspline_basis_1(value_t,2+i)  + P[i+1,:] *bspline_basis_2(value_t,2+i) + P[i+2] *bspline_basis_3(value_t,2+i)
            
            return seg_1, bspline_basis_1, bspline_basis_2, bspline_basis_3    
               
seg=[]
for i in range(N_MHE):
    
    z, bspline_basis_1, bspline_basis_2, bspline_basis_3=bspline_choose(bspline_fun_choose, X_measurement[i,0], basis_splines_table,t_segments_subpoints)  
    seg = horzcat(seg,z)#it tells us the total number of measurements vs the dimension of the control points
    print(z)
seg=seg.T
#print(seg.shape)
#print(seg)




#task is to put the coefficients of the 0.8 generated, in the segments of the 1.0 m/s



Tf = 1.0  # prediction horizon
N = 50  # number of discretization steps
T=0.8


n_steps=1 #change it for later for 3 steps
DT=Tf/N
Nsim=int(N*T/Tf)
print(T)

n_states=6
n_controls=2

m = 0.043
C1=0.5
C2=15.5
Cm1=0.28
Cm2=0.05
Cr0=0.011
Cr2=0.006
track = "LMS_Track.txt"
[s0,xref,yref,psiref,kapparef] = getTrack(track)
kapparef_s=interpolant('kapparef_s','bspline',[s0],kapparef)


X=MX.sym('X',n_states,(N_MHE+1))
st=X[:,0]
#initiliazing the last values as the starting points
#G=vertcat(st-P_meas[6*(N_MHE):6*(N_MHE)+6])

D_values=con[:,0]
der_values=con[:,1]
#forming coefficients over here
f_coefficients=Function('f_coefficients',[P_3],[seg])

model_m_t = f_coefficients(coefficients_s)
print(model_m_t.shape)
print(model_m_t)
#for i in range(500):
  #  print(model_m_t[i,:])




st = X_measurement[1,0:6]
X_estimated=np.zeros((N_MHE,6))
for k in range(0,N_MHE):  

    for i in range(n_steps): 
        model_m = model_m_t[i,:]
        
        
        model_m = np.append(model_m,[0,0])
        model_m=model_m.T
        
        #print(model_m.shape)
        #print(model_m)
        


        s = MX.sym('s')
        n = MX.sym('n')
        alpha = MX.sym('alpha')
        v   = MX.sym('v')
        D = MX.sym('D')
        delta = MX.sym('delta')
        x = vertcat(s,n,alpha,v,D,delta)

        Fxd = (Cm1 - Cm2 * v) * D - Cr2 * v * v - Cr0 * tanh(5 * v)
        sdota = ( v* np.cos( alpha+ C1* delta)) / (1- kapparef_s (s) * n)

        xdot = vertcat( sdota, \
        v * sin( alpha + C1 * delta),\
        v * C2 * delta - kapparef_s(s) * sdota,\
        Fxd/m* cos( C1* delta),\
        D_values[k],\
        der_values[k]
        )

        f_rk4 = Function('f_rk4',[x],[xdot])


        X0 = MX.sym('X0',6,1)
        U0 = MX.sym('U0',2,1)
        DT = DT/n_steps
            
        k1 = f_rk4(X0)
        k2 = f_rk4(X0 + DT/2 * k1)
        k3 = f_rk4(X0 + DT/2 * k2)
        k4 = f_rk4(X0 + DT * k3)
        X_out = X0 + DT/6*(k1 +2*k2 +2*k3 +k4)		

        integrator_fun=Function('integrator_fun',[X0],[X_out])
        
        #st_next = X[:,k+1]
        st=integrator_fun(st)
        st = integrator_fun(st) 
        st[0]= st[0]*1
        st[1]= st[1]
        st[2]= st[2]
        st[3]= st[3]
        #st[4]= st[4]*0.0001
        #st[5]= st[5]*0.5
        print(st)
        X_estimated[i:i+1,:]= st.T

        
        #st_next_euler = integrator_fun(st) + model_m 

plotTrackProj(N_MHE,X_estimated, X_measurement, s0, filename="LMS_Track.txt", T_opt='None')
plt.show() """


#X_measurement = np.loadtxt("final_data/1.1_test_data/states_final_alligned_test_1.1")
#X_measurement = np.loadtxt("final_data/1.1_test_data/states_final_alligned_test_1.1")

#X_estimated = np.loadtxt("final_data/1.0_data/estimated_data/simX_res.txt")
#con=np.loadtxt("final_data/1.1_test_data/control_final_alligned_test_1.1")


#X_measurement = np.loadtxt("final_data/1.1_data/states_final_alligned_1.2")
#X_estimated = np.loadtxt("final_data/1.2_data/estimated_data/simX_res.txt")
#con=np.loadtxt("final_data/1.2_data/control_final_alligned_1.2")


#X_measurement = np.loadtxt("final_data/1.4_data/states_final_alligned_1.4")
#X_estimated = np.loadtxt("final_data/1.0_data/estimated_data/simX_res.txt")
#con=np.loadtxt("final_data/1.4_data/control_final_alligned_1.4")


X_measurement = np.loadtxt("final_data/1.1_test_data/states_final_alligned_test_1.1")
con=np.loadtxt("final_data/1.1_test_data/control_final_alligned_test_1.1")

#X_measurement = np.loadtxt("final_data/0.8_training_data/states_final_alligned_0.8_training")
#con=np.loadtxt("final_data/0.8_training_data/control_final_alligned_0.8_training")


#con=np.loadtxt("final_data/0.8_training_data/control_final_alligned_0.8_training")
coefficeints_us=np.loadtxt('final_data/0.8_training_data/estimated_data/coeffcients_estimated.txt')
coefficients_s= np.zeros((4,258))


for i in range(258):
  #print("hi")  
  coefficients_s[:,i]  = coefficeints_us[4*i:4*i+4]
  #print(coefficeints_us[4*i:4*i+4])  





basis_splines_table = np.loadtxt("spline_tasks/basis_splines")
t_segments_subpoints = np.loadtxt("spline_tasks/t_segments_subpoints")

nr_points=256
nr_sub_points=50

t=MX.sym('t',1)
t_1=MX.sym('t_1',1)
t_2=MX.sym('t_2',1)

f_1=(t-t_2)*(t-t_2)/2
f_2=0.5*(-2* (t-t_2)*(t-t_2) + 2*((t-t_2)) + 1)
f_3=0.5*(((t-t_2)*(t-t_2)) - 2*(t-t_2) + 1)
f_4=0
bspline_basis_1=Function('bspline_basis',[t,t_2],[f_1])
bspline_basis_2=Function('bspline_basis',[t,t_2],[f_2])
bspline_basis_3=Function('bspline_basis',[t,t_2],[f_3])
bspline_basis_4=Function('bspline_basis',[t,t_2],[f_4])
P_3=MX.sym('P_3',4,nr_points+2)
"""
#forming casadi function
""" 
f_1=(t-t_2)*(t-t_2)/2
f_2=0.5*(-2* (t-t_2)*(t-t_2) + 2*((t-t_2)) + 1)
f_3=0.5*(((t-t_2)*(t-t_2)) - 2*(t-t_2) + 1)
f_4=0
bspline_basis_1=Function('bspline_basis',[t,t_2],[f_1])
bspline_basis_2=Function('bspline_basis',[t,t_2],[f_2])
bspline_basis_3=Function('bspline_basis',[t,t_2],[f_3])
bspline_basis_4=Function('bspline_basis',[t,t_2],[f_4])


#points and subrange
#nr_points=64
nr_points=256
nr_sub_points=50

bspline_pieces= []
t_segments_subpoints=np.empty((0,50),int)   
x_segments_subpoints=np.empty((0,50),int)    

P_3=MX.sym('P',4,nr_points+2)


for i in range(nr_points):
    
    t_range=np.linspace(2 +i, 3 +i, 4)
    
    t_sub_range=np.linspace(2 +i,3 +i, nr_sub_points)
    t_sub_range=np.array([t_sub_range])
    #print(t_sub_range)
    #print(t_sub_range)
    x_1=[]
    x_2=[]
    x_3=[]
    
    t_segments_subpoints=np.append(t_segments_subpoints,t_sub_range,axis=0)
    
    
    
    for j in t_sub_range:
        #print(j)
        t_1=1*bspline_basis_1(j,2+i)
        t_1=t_1.full()
        t_2=1*bspline_basis_2(j,2+i)
        t_2=t_2.full()
        t_3=1*bspline_basis_3(j,2+i)
        t_3=t_3.full()
        
        x_1=np.append(x_1,t_1)    
        x_2=np.append(x_2,t_2)
        x_3=np.append(x_3,t_3)
        x_1=np.array([x_1])
        x_2=np.array([x_2])
        x_3=np.array([x_3])
        
            
    #appending the basis functions   
    for i in range(1):
        x_segments_subpoints=np.append(x_segments_subpoints,x_1,axis=0)
        x_segments_subpoints=np.append(x_segments_subpoints,x_2,axis=0)
        x_segments_subpoints=np.append(x_segments_subpoints,x_3,axis=0)

#plt.xlabel('t')
#plt.ylabel('$N_{i,j}(t)$')
#for i in range(nr_points):
    #plt.plot(t_segments_subpoints[i,:],x_segments_subpoints[3*i,:])
    #plt.plot(t_segments_subpoints[i,:],x_segments_subpoints[3*i+1,:])
    #plt.plot(t_segments_subpoints[i,:],x_segments_subpoints[3*i+2,:])

#plt.legend    
#plt.savefig('complete_range_bspline.pdf')    
#plt.show()

s=8.71

s_points=np.linspace(0,s,nr_points-1)
s_points=np.append(s_points,0)
s_segments=np.empty([nr_points-1,2])

t_points=np.linspace(2,nr_points,nr_points-1)
t_points=np.append(t_points,0)
t_segments=np.empty((nr_points-1,2))


for i in range(nr_points-1):
    s_segments[i,:]=s_points[0+i:2+i]
    t_segments[i,:]=t_points[0+i:2+i]
    





#for knowing the that which t_segment equals which s_segment
t_s_segments=np.concatenate((t_segments,s_segments),axis=1)

#for getting the s_val and finding it's corresponding t value such that we could form the bspline segments when needed
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx],idx

 
#uses the find_nearest function to tell us the corresponding t value
def bspline_fun_choose(t_s_segments,s_val):
   
    for i in range(nr_points):
        
        if(  t_s_segments[i,2] <= s_val <= t_s_segments[i,3] ):
            #subpoint_range=(t_s_segments[i,2]-t_s_segments[i,0])
            subrange_t=np.linspace(t_s_segments[i,0],t_s_segments[i,1],nr_sub_points)
            subrange_s=np.linspace(t_s_segments[i,2],t_s_segments[i,3],nr_sub_points)
            value_s,index=find_nearest(subrange_s,s_val)
            value_t=subrange_t[index]
            
            
            bspline_segment=t_s_segments[i,:]
            key=bspline_segment[0]
            return bspline_segment ,index,value_s,value_t,key



def bspline_choose(bspline_fun_choose, s_val, t_s_segments,t_segments_subpoints):
    for i in range(nr_points):
        
        bspline_segment, index, value, value_t, key = bspline_fun_choose(t_s_segments, s_val)
        
        #print(bspline_segment)

        t_range=np.linspace(2 +i, 3 +i, 4)
        
        t_sub_range=np.linspace(2 +i,3 +i, nr_sub_points)
        t_sub_range=np.array([t_sub_range])

        x_1=[]
        x_2=[]
        x_3=[]
        
        t_segments_subpoints=np.append(t_segments_subpoints,t_sub_range,axis=0)
        
        
        
        for j in t_sub_range:
            #print(j)
            t_1=1*bspline_basis_1(j,2+i)
            t_1=t_1.full()
            t_2=1*bspline_basis_2(j,2+i)
            t_2=t_2.full()
            t_3=1*bspline_basis_3(j,2+i)
            t_3=t_3.full()
            
            x_1=np.append(x_1,t_1)    
            x_2=np.append(x_2,t_2)
            x_3=np.append(x_3,t_3)
            x_1=np.array([x_1])
            x_2=np.array([x_2])
            x_3=np.array([x_3])
            
        
        if(2+i==key):
            seg_1=P_3[:,i] *bspline_basis_1(value_t,2+i)  + P_3[:,i+1] *bspline_basis_2(value_t,2+i) + P_3[:,i+2] *bspline_basis_3(value_t,2+i)
            #seg_1=P[i,:] *bspline_basis_1(value_t,2+i)  + P[i+1,:] *bspline_basis_2(value_t,2+i) + P[i+2] *bspline_basis_3(value_t,2+i)
            
            return seg_1, bspline_basis_1, bspline_basis_2, bspline_basis_3    
               
seg=[]

time1=perf_counter()
print("for dataset 2 starting time for getting segments")
print(time1)
for i in range(N_MHE):
    
    
    z, bspline_basis_1, bspline_basis_2, bspline_basis_3=bspline_choose(bspline_fun_choose, X_measurement[i,0], basis_splines_table,t_segments_subpoints)  
    seg = horzcat(seg,z)#it tells us the total number of measurements vs the dimension of the control points
    #print(z)
seg=seg.T


time2=perf_counter()
print("for dataset 2 end time for getting segments is")
print(time2)
print(time2-time1)
    

#print(seg.shape)
#print(seg)


f_coefficients=Function('f_coefficients',[P_3],[seg])



time3=perf_counter()
print("for dataset 2 start time for mm is")
print(time3)


model_m_t = f_coefficients(coefficients_s)

time4=perf_counter()
print("for dataset 2 end time for m is")
print(time4)
print(time4-time3)




#print(model_m_t)

track="LMS_Track.txt"
Tf = 1.0  # prediction horizon
N = 50  # number of discretization steps
T =3 # maximum simulation time[s]
DT=Tf/N
n_steps=1
DT=DT/1
Nsim=int(N*T/Tf)

m = 0.043
C1=0.5
C2=15.5
Cm1=0.28
Cm2=0.05
Cr0=0.011
Cr2=0.006

track = "LMS_Track.txt"
[s0,xref,yref,psiref,kapparef] = getTrack(track)
kapparef_s=interpolant('kapparef_s','bspline',[s0],kapparef)


#dynamic system

s = MX.sym('s')
n = MX.sym('n')
alpha = MX.sym('alpha')
v   = MX.sym('v')
D = MX.sym('D')
delta = MX.sym('delta')
x = vertcat(s,n,alpha,v,D,delta)

# controls

derD = MX.sym('derD')
derDelta = MX.sym('derDelta')
u = vertcat(derD,derDelta)

Fxd = (Cm1 - Cm2 * v) * D - Cr2 * v * v - Cr0 * tanh(5 * v)
sdota = ( v* np.cos( alpha+ C1* delta)) / (1- kapparef_s (s) * n)


xdot = vertcat( sdota,\
v * sin( alpha + C1 * delta),\
v * C2 * delta - kapparef_s(s) * sdota,\
Fxd/m* cos( C1* delta),\
derD,\
derDelta
)


#applying rk4 function scheme here
f_rk4 = Function('f_rk4',[x,u],[xdot])

X0 = MX.sym('X0',6,1)
U0 = MX.sym('U0',2,1)
    
k1 = f_rk4(X0, U0)
k2 = f_rk4(X0 + DT/2 * k1, U0)
k3 = f_rk4(X0 + DT/2 * k2, U0)
k4 = f_rk4(X0 + DT * k3, U0)
X_out = X0 + DT/6*(k1 +2*k2 +2*k3 +k4)		

integrator_fun=Function('integrator_fun',[X0,U0],[X_out])

simX=np.zeros((N_MHE,6))

x0=X_measurement[0,0:6]

simX[0,:]=x0
for i in range(N_MHE):
    model_m = model_m_t[i,:]    
    model_m = np.append(model_m,[0,0])
    model_m=model_m.T
        
    x0=integrator_fun(x0,con[i,0:2]) + model_m						
    x0 = np.ravel(x0)
    #print(x0)
                
    simX[i+1:,0:6] = x0
    #print(x0)

#plotTrackProj(N_MHE,simX, X_measurement, s0, filename="LMS_Track.txt", T_opt='None')

plt.show()
        


#print(model_m_t.shape)




def Function_test(simX,x0,idx):

    #simX = np.zeros((N_MHE,6))
    #x0 = X_measurement[0,0:6]
    simX[0,0:6]=x0
    def find_nearest(measurements,value):
        s_meas=measurements[:,0]#taking the s _measurements
        idx=(np.abs(s_meas-value)).argmin()
        return measurements[idx],idx #returning measurements at least difference

    out_of_bound_c=0
    for i in range(N_MHE-1):
        
        if(idx+i==N_MHE):
            break

        
        if (x0[1] < 0.16 or x0[1] > -0.16):    

            
                    
            model_m = model_m_t[idx+i,:]    
            model_m = np.append(model_m,[0,0])
            model_m=model_m.T
                
            x0=integrator_fun(x0,con[i+idx,0:2]) + model_m						
            x0 = np.ravel(x0)
            #print(x0)
                        
            simX[i+1:,0:6] = x0

        if (x0[1] > 0.16 or x0[1] < -0.16):
            #print("break point")
            #print(x0)
            meas_len=X_measurement.shape[0]
            x0_s=x0[0]
            closest,idx = find_nearest(X_measurement,x0_s)
            near_point = closest[0:6]
            x0 = near_point
            out_of_bound_c= out_of_bound_c + 1
            
        #simX[i+1:,:]=x0
        
            break

    return simX,x0,idx

idx=0
simX_0 = np.zeros((N_MHE,6))
x0 = X_measurement[0,0:6]
simX_0,x0,idx = Function_test(simX_0,x0,idx)
#print(simX)
#plotTrackProj(N_MHE,simX_0, X_measurement, s0, filename="LMS_Track.txt", T_opt='None')
#plt.show()

#print(idx)

#for i in range(N_MHE):
    #print(simX_0[i,:])
#print(simX_0)

simX_1 = np.zeros((N_MHE,6))
#print("xo value is")
#print(x0)
simX_1,x0,idx=Function_test(simX_1,x0,idx)
#plotTrackProj(N_MHE,simX_1, X_measurement, s0, filename="LMS_Track.txt", T_opt='None')
#plt.show()
#print(idx)


simX_2 = np.zeros((N_MHE,6))
#print("xo value is")
#print(x0)
simX_2,x0,idx=Function_test(simX_2,x0,idx)
#plotTrackProj(N_MHE,simX_2, X_measurement, s0, filename="LMS_Track.txt", T_opt='None')
#plt.show()
#print(idx)


simX_3 = np.zeros((N_MHE,6))
#print("xo value is")
#print(x0)
simX_3,x0,idx=Function_test(simX_3,x0,idx)
#plotTrackProj(N_MHE,simX_3, X_measurement, s0, filename="LMS_Track.txt", T_opt='None')
#plt.show()
#print(idx)


simX_4 = np.zeros((N_MHE,6))
#print("xo value is")
#print(x0)
simX_4,x0,idx=Function_test(simX_4,x0,idx)
#plotTrackProj(N_MHE,simX_0,simX_1,simX_2,simX_3,simX_4, X_measurement, s0, filename="LMS_Track.txt", T_opt='None')
#plt.show()
#print(idx)

simX_5 = np.zeros((N_MHE,6))
#print("xo value is")
#print(x0)
simX_5,x0,idx=Function_test(simX_5,x0,idx)
#plotTrackProj(N_MHE,simX_5, X_measurement, s0, filename="LMS_Track.txt", T_opt='None')
#plt.show()

#print(idx)


simX_6 = np.zeros((N_MHE,6))
#print("xo value is")
#print(x0)
simX_6,x0,idx=Function_test(simX_6,x0,idx)



#simX_6 = np.zeros((N_MHE,6))
#print("xo value is")
#print(x0)
#simX_6,x0,idx=Function_test(simX_6,x0,idx)



plotTrackProj(N_MHE, simX_0, simX_1, simX_2, simX_3, simX_4, simX_5,simX_6, X_measurement, s0, filename="LMS_Track.txt", T_opt='None')
plt.show()


simX_res=np.concatenate((simX_0,simX_1,simX_2,simX_3,simX_4,simX_5),axis=0)
#print(simX_res.shape)
#simX_res=np.loadtxt('final_data/1.1_test_data/estimated_data/simX_res.txt')
#np.savetxt('final_data/1.1_test_data/estimated_data/simX_res.txt',simX_res)

#plots_individual_states(N_MHE,X_measurement,simX_res)
#plt.show()