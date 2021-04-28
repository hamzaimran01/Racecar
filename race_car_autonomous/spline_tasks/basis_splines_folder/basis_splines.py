import sys
sys.path.append("../THESIS")
import numpy as np
from casadi import *
from spline_tasks.basis_splines_folder.spline_casadi_parameters import *
import matplotlib.pyplot as plt
from time import perf_counter 



np.set_printoptions(precision=5,
                       threshold=10000,
                       linewidth=150)
#############################################################################################
""" this part is for forming the bspline basis functions  """
n_controls = 2
n_states = 6
basis_splines_table = np.loadtxt("spline_tasks/basis_splines")
t_segments_subpoints = np.loadtxt("spline_tasks/t_segments_subpoints")

#forming casadi variables

""" t=MX.sym('t',1)
t_1=MX.sym('t_1',1)
t_2=MX.sym('t_2',1)

#forming casadi function

f_1=(t-t_2)*(t-t_2)/2
f_2=0.5*(-2* (t-t_2)*(t-t_2) + 2*((t-t_2)) + 1)
f_3=0.5*(((t-t_2)*(t-t_2)) - 2*(t-t_2) + 1)
f_4=0
bspline_basis_1=Function('bspline_basis',[t,t_2],[f_1])
bspline_basis_2=Function('bspline_basis',[t,t_2],[f_2])
bspline_basis_3=Function('bspline_basis',[t,t_2],[f_3])
bspline_basis_4=Function('bspline_basis',[t,t_2],[f_4])
 """

#points and subrange
#nr_points=64
nr_points=256
nr_sub_points=50

bspline_pieces= []
t_segments_subpoints=np.empty((0,50),int)   
x_segments_subpoints=np.empty((0,50),int)    

P_3=MX.sym('P',3,nr_points+2)


for i in range(nr_points):
    
    #t_range=np.linspace(2 +i, 3 +i, 4)
    #print(t_range)
    
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

for i in range(nr_points):
    plt.plot(t_segments_subpoints[i,:],x_segments_subpoints[3*i,:])
    plt.plot(t_segments_subpoints[i,:],x_segments_subpoints[3*i+1,:])
    plt.plot(t_segments_subpoints[i,:],x_segments_subpoints[3*i+2,:])
#plt.show()
#plt.show()
#defined for 1/4th of the track
s=8.71
#s=8.6

s_points=np.linspace(0,s,nr_points-1)
s_points=np.append(s_points,0)
s_segments=np.empty([nr_points-1,2])

t_points=np.linspace(2,nr_points,nr_points-1)
t_points=np.append(t_points,0)
t_segments=np.empty((nr_points-1,2))


for i in range(nr_points-1):
    s_segments[i,:]=s_points[0+i:2+i]
    
    t_segments[i,:]=t_points[0+i:2+i]
    


#for knowing the that which t_segment basis equals which s_segment basis range
t_s_segments=np.concatenate((t_segments,s_segments),axis=1)


np.savetxt("spline_tasks/basis_splines_folder/saved/basis_splines",t_s_segments)
#np.savetxt("spline_tasks/basis_splines_folder/saved/t_segments_subpoints",t_segments_subpoints)


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



#loops through the basis function again, with the table of t and s segments such that it can accept a s_value and than form
#its corresponding t_value basis function and the control points(coefficients relevent to it)


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


X_measurement=np.loadtxt("final_data/states_final_alligned")
con=np.loadtxt("final_data/control_final_alligned") 


#sref,xref,yref,psiref,kapparef=getTrack("LMS_Track.txt")




#N_MHE being the time horizon until the open loop goes out of bound
N_MHE=980

X_meas=X_measurement[0:N_MHE,0:n_states]

seg=[]


""" accumulated the relevent bsplines and concatenate them w.r.t to y_meas[i,0] that is the s_value of the experiment.
-> seg= [d_s,] *[basis_1 + basis_2 +basis*3] + [d_n,] *[basis_1 + basis_2 +basis*3] + [d_alpha,] *[basis_1 + basis_2 +basis*3]
->3*1 array for seg appended in each loop

 """
t1_start = perf_counter()
print(t1_start) 
for i in range(N_MHE):
   
   
    z, bspline_basis_1, bspline_basis_2, bspline_basis_3=bspline_choose(bspline_fun_choose, X_meas[i,0], basis_splines_table,t_segments_subpoints)  
    seg = horzcat(seg,z)

    
   
seg=seg.T    
t1_end = perf_counter()
print(t1_end)
#print(seg.shape)


