#from spline_utils import
import sys
sys.path.append("../THESIS")
import numpy as np
from casadi import *
import matplotlib.pyplot as plt
from spline_tasks.casadi_integrator import casadi_integrator_fun as ci_f
from spline_tasks.readDataFcn import getTrack
from plots import plots



# i defines the the i'th index, it equals the total number of control points,now we can take the first control point only 

np.set_printoptions(precision=5,
                       threshold=10000,
                       linewidth=150)
#############################################################################################
""" this part is for forming the bspline basis functions  """


#forming casadi variables

t=MX.sym('t',1)
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


#points and subrange
nr_points=64
nr_sub_points=50

bspline_pieces= []
t_segments_subpoints=np.empty((0,50),int)   
x_segments_subpoints=np.empty((0,50),int)    
P=SX.sym('P',nr_points+2)
P_3=SX.sym('P',3,nr_points+2)


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






#defined for 1/4th of the track
s=2.175

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
               
            
           



##################################################################################################

""" 
creating an optimization problem based on the bspline segments the s_value falls.the s_value is the experimental s_value
"""


#loading the measurements,controls not used till now 
y_measurement=np.loadtxt("final_data/states_final_alligned")
con=np.loadtxt("final_data/control_final_alligned") 
y_measurement=y_measurement[46:150,:]
con=con[46:150,:]

sref,xref,yref,psiref,kapparef=getTrack("LMS_Track.txt")


DT=0
#casasdi open loop function without the disturbances such that we can have for a time horizon all the values
simX,measurements,s0,xref,yref,psiref=ci_f(y_measurement,con,DT)

#N_MHE being the time horizon until the open loop goes out of bound
N_MHE=40
con=con[0:N_MHE,:]
y_meas=y_measurement[0:N_MHE,:]

seg=[]
count=0

""" accumulated the relevent bsplines and concatenate them w.r.t to y_meas[i,0] that is the s_value of the experiment.
-> seg= [d_s,] *[basis_1 + basis_2 +basis*3] + [d_n,] *[basis_1 + basis_2 +basis*3] + [d_alpha,] *[basis_1 + basis_2 +basis*3]
->3*1 array for seg appended in each loop

 """
for i in range(N_MHE):
    
    z, bspline_basis_1, bspline_basis_2, bspline_basis_3=bspline_choose(bspline_fun_choose, y_meas[i,0], t_s_segments,t_segments_subpoints)  
    seg = horzcat(seg,z)
    print(seg.shape)
    
    

seg=seg.T



""" 
forming decision variables
->P_meas consists of 6 measurement values and 6 simulated values for a unit of time expanded through the horizon
"""


P_meas=SX.sym('P_meas',12*(N_MHE-1)+12)#(measurements + simulation)+ state for starting and 
#state for measurement which is the same

P_meas=P_meas.T


obj=0
V=DM([[0.000000000000001, 0, 0, 0, 0, 0], [0, 0.0001, 0 , 0, 0, 0],[0, 0, 0.00000001, 0, 0, 0], [0, 0, 0, 0.0000001, 0, 0], [0, 0, 0, 0, 0.0001, 0], [0, 0, 0, 0, 0, 0.000005]])


#this includes d_s,d_n,d_alpha for a single semgent control 

"""
Formation of the correspoding optimization problem


  """
for k in range(N_MHE):
    #0-6 - measurements , 6-12 -simulation results stored
    y_meas = P_meas[12*k  : 12*k+6 ]
    x_sim = P_meas[12*k +6 :12*k+ 12]
    dis=seg[k,:]
    #print(dis)
    dis=dis.T
    dis=vertcat(dis, 0, 0, 0)
    
    
    obj = obj+mtimes(mtimes((y_meas.T-x_sim.T -dis).T,V),(y_meas.T-x_sim.T -dis))
    




   
#66*3

#opt_variables are the coefficients defined i the start. which make up the segments  
   
OPT_variables = vertcat(reshape((P_3),3*(nr_points+2),1))
nlp = {'f':obj,'x':OPT_variables,'p':P_meas}

solver=nlpsol('solver','ipopt',nlp)







y_measurement=y_measurement[0:N_MHE,0:6]
y_m=y_measurement.flatten('C') #all in one line
simX_m=simX.flatten('C')


p=np.zeros((12*(N_MHE-1)+12))

args={}


""" ->filling the P_measurement with the measurement and the open loop simulated values
    ->so 0:6 are the measurement, 6:12 are the simulated values """

for k in range(N_MHE):
    
    p[12*k  : 12*k+6]=y_m[6*k:6*k+6]
    p[12*k +6  : 12*k+12 ]=simX_m[6*k :6*k+6]
    
    
#saves the measurements and the simulation values
args['p']= p

"""

-> initilization of state decision variables=total seubpoints on the grid*[d_s,d_n.d_alpha] = 0 
"""

X0_1=np.zeros(((nr_points+2,3)))


X0_1[:,0]=0
X0_1[:,1]=0
X0_1[:,2]=0


#initilization value of the optimization problem
args['x0']=reshape((X0_1.T),3*(nr_points+2),1)

#solving the optimization proble,
sol=solver(x0=args['x0'],p=args['p'])

#
""" -> try it for one segment and the second segment,
-> to check it's authenticity put it in the optimization problem back again with the bspline and plot it on the graph """



P_cal=np.zeros((3,66))

""" 
forming a function depending on the ouput of the optimization problem i.e the control points(Coefficients) for 
the numerical value i.e of dimension [3,1] at it's corresponding bspline segment
->"basis_choose_cp" function
 """

for i in range(66):
    P_cal[ :, i:i+1] = sol['x'][3*i:3*i+3]


basis_choose_cp=Function('basis_choose', [P_3], [seg])
val=basis_choose_cp(P_cal)



simX_res=y_measurement[0:1,:]


""" to find if our open loop simulation + disturbance is close to the measurement or not
-> basically reverse engineering here and savinf it in simX_res 
"""
for k in range(N_MHE):
    val_choose=np.zeros((1,6))
    val_choose[:,0:3]=val[k,:]
  
    
    new_s=simX[k:k+1,:].T + val_choose.T
    new_s=new_s.T
    simX_res = np.append(simX_res,new_s,axis=0)
    
""" plotting the three different graphs """   
#plots(simX_res,0,s0,xref,yref,psiref,keyword="si")
#plt.show()
#plots(simX,0,s0,xref,yref,psiref,keyword="open loop")
#plt.show()   
#plots(y_measurement,0,s0,xref,yref,psiref,keyword="measurements")
#plt.show()
#plots(simX,0,s0,xref,yref,psiref)


#print(simX.shape)
#print(simX)
#print(simX_res)
print(y_measurement)

