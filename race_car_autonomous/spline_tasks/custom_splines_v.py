#from spline_utils import
import sys
sys.path.append("../THESIS")
import numpy as np
from casadi import *
#from spline_tasks.basis_splines_folder.spline_casadi_parameters import *
#from spline_tasks.optimization_initilizations import *
import matplotlib.pyplot as plt
from spline_tasks.casadi_integrator import casadi_integrator_fun as ci_f
from spline_tasks.readDataFcn import getTrack
from plots import plotTrackProj
#from spline_tasks.basis_splines_folder.basis_splines import seg
from time import perf_counter 
#from spline_tasks.optimization_initilizations import *
from spline_tasks.model.tracks.readDataFcn import getTrack
from matplotlib import rc
import time



basis_splines_table  =  np.loadtxt("spline_tasks/basis_splines")
t_segments_subpoints  =  np.loadtxt("spline_tasks/t_segments_subpoints")




Tf  =  1.0  # prediction horizon
N  =  50  # number of discretization steps
T = 0.8


n_steps = 1 #change it for later for 3 steps
DT = Tf/N
Nsim = int(N*T/Tf)
print(T)

n_states = 6
n_controls = 2

m  =  0.043
C1 = 0.5
C2 = 15.5
Cm1 = 0.28
Cm2 = 0.05
Cr0 = 0.011
Cr2 = 0.006
track  =  "LMS_Track.txt"
[s0,xref,yref,psiref,kapparef]  =  getTrack(track)
kapparef_s = interpolant('kapparef_s','bspline',[s0],kapparef)





# i defines the the i'th index, it equals the total number of control points,now we can take the first control point only 

np.set_printoptions(precision = 5,
                       threshold = 10000,
                       linewidth = 150)
#############################################################################################
""" this part is for forming the bspline basis functions  """


#forming casadi variables
#points and subrange
#nr_points = 64
nr_points = 256
nr_sub_points = 50

t = MX.sym('t',1)
t_1 = MX.sym('t_1',1)
t_2 = MX.sym('t_2',1)

f_1 = (t-t_2)*(t-t_2)/2
f_2 = 0.5*(-2* (t-t_2)*(t-t_2) + 2*((t-t_2)) + 1)
f_3 = 0.5*(((t-t_2)*(t-t_2)) - 2*(t-t_2) + 1)
f_4 = 0
bspline_basis_1 = Function('bspline_basis',[t,t_2],[f_1])
bspline_basis_2 = Function('bspline_basis',[t,t_2],[f_2])
bspline_basis_3 = Function('bspline_basis',[t,t_2],[f_3])
bspline_basis_4 = Function('bspline_basis',[t,t_2],[f_4])
P_3 = MX.sym('P',4,nr_points+2)
"""
#forming casadi function
""" 
f_1 = (t-t_2)*(t-t_2)/2
f_2 = 0.5*(-2* (t-t_2)*(t-t_2) + 2*((t-t_2)) + 1)
f_3 = 0.5*(((t-t_2)*(t-t_2)) - 2*(t-t_2) + 1)
f_4 = 0
bspline_basis_1 = Function('bspline_basis',[t,t_2],[f_1])
bspline_basis_2 = Function('bspline_basis',[t,t_2],[f_2])
bspline_basis_3 = Function('bspline_basis',[t,t_2],[f_3])
bspline_basis_4 = Function('bspline_basis',[t,t_2],[f_4])


#points and subrange
#nr_points = 64
nr_points = 256
nr_sub_points = 50

bspline_pieces =  []
t_segments_subpoints = np.empty((0,50),int)   
x_segments_subpoints = np.empty((0,50),int)    
P = MX.sym('P',nr_points+2)
P_3 = MX.sym('P',4,nr_points+2)
figure = plt.gcf()
figure.set_size_inches(24, 18)

for i in range(nr_points):
    
    t_range = np.linspace(2 +i, 3 +i, 4)
    
    t_sub_range = np.linspace(2 +i,3 +i, nr_sub_points)
    t_sub_range = np.array([t_sub_range])
    #print(t_sub_range)
    #print(t_sub_range)
    x_1 = []
    x_2 = []
    x_3 = []
    
    t_segments_subpoints = np.append(t_segments_subpoints,t_sub_range,axis = 0)
    
    
    
    for j in t_sub_range:
        #print(j)
        t_1 = 1*bspline_basis_1(j,2+i)
        t_1 = t_1.full()
        t_2 = 1*bspline_basis_2(j,2+i)
        t_2 = t_2.full()
        t_3 = 1*bspline_basis_3(j,2+i)
        t_3 = t_3.full()
        
        x_1 = np.append(x_1,t_1)    
        x_2 = np.append(x_2,t_2)
        x_3 = np.append(x_3,t_3)
        x_1 = np.array([x_1])
        x_2 = np.array([x_2])
        x_3 = np.array([x_3])
        
            
    #appending the basis functions   
    for i in range(1):
        x_segments_subpoints = np.append(x_segments_subpoints,x_1,axis = 0)
        x_segments_subpoints = np.append(x_segments_subpoints,x_2,axis = 0)
        x_segments_subpoints = np.append(x_segments_subpoints,x_3,axis = 0)

#plt.xlabel('t')
#plt.ylabel('$N_{i,j}(t)$')
#for i in range(nr_points):
    #plt.plot(t_segments_subpoints[i,:],x_segments_subpoints[3*i,:])
    #plt.plot(t_segments_subpoints[i,:],x_segments_subpoints[3*i+1,:])
    #plt.plot(t_segments_subpoints[i,:],x_segments_subpoints[3*i+2,:])

#
# plt.legend()

#plt.savefig('complete_range_bspline.pdf',bbox_inches = 'tight',dpi = 100)    
#plt.show()
#plt.show()
#defined for 1/4th of the track
s = 8.71
#s = 8.6

s_points = np.linspace(0,s,nr_points-1)
s_points = np.append(s_points,0)
s_segments = np.empty([nr_points-1,2])

t_points = np.linspace(2,nr_points,nr_points-1)
t_points = np.append(t_points,0)
t_segments = np.empty((nr_points-1,2))


for i in range(nr_points-1):
    s_segments[i,:] = s_points[0+i:2+i]
    t_segments[i,:] = t_points[0+i:2+i]
    





#for knowing the that which t_segment equals which s_segment
t_s_segments = np.concatenate((t_segments,s_segments),axis = 1)

#for getting the s_val and finding it's corresponding t value such that we could form the bspline segments when needed
def find_nearest(array, value):
    array  =  np.asarray(array)
    idx  =  (np.abs(array - value)).argmin()
    return array[idx],idx

 
#uses the find_nearest function to tell us the corresponding t value
def bspline_fun_choose(t_s_segments,s_val):
   
    for i in range(nr_points):
        
        if(  t_s_segments[i,2] <=  s_val <=  t_s_segments[i,3] ):
            #subpoint_range = (t_s_segments[i,2]-t_s_segments[i,0])
            subrange_t = np.linspace(t_s_segments[i,0],t_s_segments[i,1],nr_sub_points)
            subrange_s = np.linspace(t_s_segments[i,2],t_s_segments[i,3],nr_sub_points)
            value_s,index = find_nearest(subrange_s,s_val)
            value_t = subrange_t[index]
            
            
            bspline_segment = t_s_segments[i,:]
            key = bspline_segment[0]
            return bspline_segment ,index,value_s,value_t,key



#loops through the basis function again, with the table of t and s segments such that it can accept a s_value and than form
#its corresponding t_value basis function and the control points(coefficients relevent to it)


def bspline_choose(bspline_fun_choose, s_val, t_s_segments,t_segments_subpoints):
    for i in range(nr_points):
        
        bspline_segment, index, value, value_t, key  =  bspline_fun_choose(t_s_segments, s_val)
        
        #print(bspline_segment)

        t_range = np.linspace(2 +i, 3 +i, 4)
        
        t_sub_range = np.linspace(2 +i,3 +i, nr_sub_points)
        t_sub_range = np.array([t_sub_range])

        x_1 = []
        x_2 = []
        x_3 = []
        
        t_segments_subpoints = np.append(t_segments_subpoints,t_sub_range,axis = 0)
        
        
        
        for j in t_sub_range:
            #print(j)
            t_1 = 1*bspline_basis_1(j,2+i)
            t_1 = t_1.full()
            t_2 = 1*bspline_basis_2(j,2+i)
            t_2 = t_2.full()
            t_3 = 1*bspline_basis_3(j,2+i)
            t_3 = t_3.full()
            
            x_1 = np.append(x_1,t_1)    
            x_2 = np.append(x_2,t_2)
            x_3 = np.append(x_3,t_3)
            x_1 = np.array([x_1])
            x_2 = np.array([x_2])
            x_3 = np.array([x_3])
            
        
        if(2+i == key):
            seg_1 = P_3[:,i] *bspline_basis_1(value_t,2+i)  + P_3[:,i+1] *bspline_basis_2(value_t,2+i) + P_3[:,i+2] *bspline_basis_3(value_t,2+i)
            #seg_1 = P[i,:] *bspline_basis_1(value_t,2+i)  + P[i+1,:] *bspline_basis_2(value_t,2+i) + P[i+2] *bspline_basis_3(value_t,2+i)
            
            
            return seg_1, bspline_basis_1, bspline_basis_2, bspline_basis_3    
               
            
           



##################################################################################################

""" 
creating an optimization problem based on the bspline segments the s_value falls.the s_value is the experimental s_value
"""


#loading the measurements,controls not used till now 
#X_measurement = np.loadtxt("final_data/1.0_data/states_final_alligned_1.0")
#con = np.loadtxt("final_data/1.0_data/control_final_alligned_1.0") 

#X_measurement = np.loadtxt("final_data/0.8_training_data/states_final_alligned_1.0")
#con = np.loadtxt("final_data/0.8_training_data/control_final_alligned_1.0") 


#X_measurement = np.loadtxt("final_data/1.4_data/states_final_alligned_1.4")
#con = np.loadtxt("final_data/1.4_data/control_final_alligned_1.4") 

X_measurement = np.loadtxt("final_data/0.8_training_data/states_final_alligned_0.8_training")
con = np.loadtxt("final_data/0.8_training_data/control_final_alligned_0.8_training") 


#sref,xref,yref,psiref,kapparef = getTrack("LMS_Track.txt")




#N_MHE being the time horizon until the open loop goes out of bound
#N_MHE = 980
#N_MHE = 779 #for 0.8
#N_MHE = 1000
N_MHE = 1625
#N_MHE = 1770 #for 1.0
#N_MHE = 596 # for 1.4
con = con[0:N_MHE,0:n_controls]
X_meas = X_measurement[0:N_MHE,0:n_states]
print(X_measurement.shape)




seg = []
count = 0

""" accumulated the relevent bsplines and concatenate them w.r.t to y_meas[i,0] that is the s_value of the experiment.
-> seg =  [d_s,] *[basis_1 + basis_2 +basis*3] + [d_n,] *[basis_1 + basis_2 +basis*3] + [d_alpha,] *[basis_1 + basis_2 +basis*3]
->3*1 array for seg appended in each loop

 """


t1_start  =  perf_counter()
print(t1_start) 

for i in range(N_MHE):
    
    z, bspline_basis_1, bspline_basis_2, bspline_basis_3 = bspline_choose(bspline_fun_choose, X_meas[i,0], basis_splines_table,t_segments_subpoints)  
    seg  =  horzcat(seg,z)

t1_end  =  perf_counter()
print(t1_end)

print("the segmentation accumulation time is")
print(t1_start-t1_end)
t_1_2 = perf_counter()
print(t_1_2)
P_meas = MX.sym('P_meas',6*(N_MHE)+6)#(measurements )+state for starting  

#state for measurement which is the same

""" 
-> formation of multiple shooting problem
"""
#multiple shooting


D_values = con[:,0]
der_values = con[:,1]
G = []

#states optimization variables along the horizon
X = MX.sym('X',n_states,(N_MHE))
st = X[:,0]
#initiliazing the last values as the starting points
#G = vertcat(st-P_meas[6*(N_MHE):6*(N_MHE)+6])
G = vertcat(st-P_meas[0:6])
D_values = con[:,0]
der_values = con[:,1]

start2 = perf_counter()
print("start 2 value is")
print(start2)


for k in range(0,N_MHE-1):  

    for i in range(n_steps): 
        dis  =  seg[:,k]

        #dis  =  dis.T
        print(dis.shape)
        dis  =  vertcat(dis, 0, 0)
        print(dis.shape)
        st  =  X[:,k]

        s  =  MX.sym('s')
        n  =  MX.sym('n')
        alpha  =  MX.sym('alpha')
        v    =  MX.sym('v')
        D  =  MX.sym('D')
        delta  =  MX.sym('delta')
        x  =  vertcat(s,n,alpha,v,D,delta)

        Fxd  =  (Cm1 - Cm2 * v) * D - Cr2 * v * v - Cr0 * tanh(5 * v)
        sdota  =  ( v* np.cos( alpha+ C1* delta)) / (1- kapparef_s (s) * n)

        xdot  =  vertcat( sdota,\
        v * sin( alpha + C1 * delta),\
        v * C2 * delta - kapparef_s(s) * sdota,\
        Fxd/m* cos( C1* delta),\
        D_values[k],\
        der_values[k]
        )

        f_rk4  =  Function('f_rk4',[x],[xdot])


        X0  =  MX.sym('X0',6,1)
        U0  =  MX.sym('U0',2,1)
        DT  =  DT/n_steps
            
        k1  =  f_rk4(X0)
        k2  =  f_rk4(X0 + DT/2 * k1)
        k3  =  f_rk4(X0 + DT/2 * k2)
        k4  =  f_rk4(X0 + DT * k3)
        X_out  =  X0 + DT/6*(k1 +2*k2 +2*k3 +k4)		

        integrator_fun = Function('integrator_fun',[X0],[X_out])

        st_next  =  X[:,k+1]

        st_next_euler  =  integrator_fun(st) + dis 
            
    G = vertcat(G,st_next-st_next_euler)
        


""" 
forming decision variables
->P_meas consists of 6 measurement values and 6 simulated values for a unit of time expanded through the horizon
"""


obj = 0
#V = DM([[0.0001, 0, 0, 0, 0, 0], [0, 0.0000001, 0 , 0, 0, 0],[0, 0, 0.0001, 0, 0, 0], [0, 0, 0, 0.000001, 0, 0], [0, 0, 0, 0, 0.0001, 0], [0, 0, 0, 0, 0, 0.000005]])
V = DM([[0.0000001, 0, 0, 0, 0, 0], [0, 0.0001, 0 , 0, 0, 0],[0, 0, 0.00000001, 0, 0, 0], [0, 0, 0, 0.0000001, 0, 0], [0, 0, 0, 0, 0.0001, 0], [0, 0, 0, 0, 0, 0.000005]])

#this includes d_s,d_n,d_alpha for a single semgent control 

"""
Formation of the correspoding optimization problem

  """
print("seg shape is")
print(seg.shape)
for k in range(0,N_MHE-1 ):
    
    X_m = P_meas[6*k:6*k+6]
    st = X[:,k]
    dis = seg[:,k]
    
    #dis = dis.T
    dis = vertcat(dis, 0, 0, 0)
    
    obj  =  obj+mtimes(mtimes((X_m-st).T,V),(X_m-st))
   
#66*3

#opt_variables are the coefficients defined i the start. which make up the segments  

#consists of x variables and the bspline points
OPT_variables  =  vertcat(reshape(X,6*(N_MHE),1),reshape((P_3),4*(nr_points+2),1))
nlp  =  {'f':obj,'x':OPT_variables,'g':G, 'p':P_meas}


solver = nlpsol('solver','ipopt',nlp)

p = np.zeros( (6*(N_MHE) + 6) )

for k in range(N_MHE):
    
    p[6*k:6*k+6]  =  X_meas[k,:]
    
    


args = {}
args['p'] = p

lbg = np.zeros((6*(N_MHE)))
ubg = np.zeros((6*(N_MHE)))
#ubg[0:4*(N):1] = 0
#lbg[0:4*(N):1] = 0
args['lbg'] = lbg
args['ubg'] = ubg

ubx = np.zeros((6*(N_MHE)))
lbx = np.zeros((6*(N_MHE)))
#s
lbx[0:6*(N_MHE):6] = -inf
ubx[0:6*(N_MHE):6] = +inf
#n
lbx[1:6*(N_MHE):6] = -0.12
ubx[1:6*(N_MHE):6] = +0.12
#alpha
lbx[2:6*(N_MHE):6] = -inf
ubx[2:6*(N_MHE):6] = +inf
#v
lbx[3:6*(N_MHE):6] = -inf
ubx[3:6*(N_MHE):6] = +inf
#deltaD

lbx[4:6*(N_MHE):6] = -10
ubx[4:6*(N_MHE):6] = +10

#delta delta
lbx[5:6*(N_MHE):6] = -2
ubx[5:6*(N_MHE):6] = +2


lbx = np.append(lbx,-inf*np.ones((4*(nr_points+2))))
ubx = np.append(ubx,inf*np.ones((4*(nr_points+2))))


args['lbx'] = lbx
args['ubx'] = ubx

#X0_1  =  np.zeros(((N_MHE+1,6)))
#print(N_MHE)
#print(X0_1.shape)
X0_1  =  X_measurement[0:N_MHE,0:n_states]
DIS_1 = np.zeros(((nr_points+2,4)))

args['x0'] = np.append(reshape((X0_1.T),6*(N_MHE),1),reshape((DIS_1.T),4*(nr_points +2),1))

sol = solver(x0 = args['x0'],p = args['p'],lbg = args['lbg'],ubg = args['ubg'],lbx = args['lbx'],ubx = args['ubx'])



start3 = perf_counter()
print("start 3 value is")
print(start3)

print("time taken for the solver problem is")
print(start3-start2)


""" 
-> sol['x'] first 246 values are the 41 states 
->rest are the control points
 """





""" 
forming a function depending on the ouput of the optimization problem i.e the control points(Coefficients) for 
the numerical value i.e of dimension [3,1] at it's corresponding bspline segment
->"basis_choose_cp" function
 """

simX_res = X_measurement[0:1,0:6]
simX_a = X_measurement[0:1,0:6]



  
for k in range(N_MHE):
    
    simX  =  sol['x'][6*k:6*k+6]
    #print("simulation value is")
    
    #print(val_choose)
    #print(val_choose)


    #previously added val_choose over here!!!!!!

    new_s = simX 
    
    #print(new_s.shape)
    new_s = new_s.T
    simX = simX.T
    simX_a  =  np.append(simX_a,simX,axis = 0)
    simX_res  =  np.append(simX_res,new_s,axis = 0)

#print(simX_res)

print("states are")
states  =   sol['x'][0:6*(N_MHE+1)]
print(states)
coefficients_estimated  =  sol['x'][6*N_MHE:6*N_MHE+4*(nr_points +2)]
print("estimated coefficients are")
print(coefficients_estimated)
print(simX_res.shape)

np.savetxt('final_data/0.8_training_data/estimated_data/coeffcients_estimated.txt',coefficients_estimated)
plotTrackProj(N_MHE,simX_res ,0,0,0,0,0,simX_a, s0, filename = 'LMS_Track.txt', T_opt = 'None')

#np.savetxt('final_data/1.1_test_data/estimated_data/simX_res.txt',simX_res)
np.savetxt('final_data/0.8_training_data/estimated_data/simX_res.txt',simX_res)

plt.show()
#print(X_meas)













