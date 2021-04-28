import sys
sys.path.append("../THESIS")

import numpy as np
from casadi import *
from plots import plots
from spline_tasks.readDataFcn import getTrack
from spline_tasks.casadi_integrator import casadi_integrator_fun

""" 
->storing the measurements, controls(done) 
->mapping it onto the track(done)
"""
#track=np.loadtxt("spline_tasks/LMS_Track.txt")
y_meas=np.loadtxt("final_data/states_final_alligned")
con=np.loadtxt("final_data/control_final_alligned") 
y_meas=y_meas[0:150,:]
con=con[0:150,:]

sref,xref,yref,psiref,kapparef=getTrack("LMS_Track.txt")

""" 
x_measurements=0
plots(y_meas,x_measurements,sref,xref,yref,psiref)
"""

""" 
s recorded(y_measutement is from 0.4148 to 1.63)
 """

""" 
-> import the intergrator casadi(done)
->make the measurements as the  sample points for s 
-> create the optimization problem of a single bspline joints 
"""

DT=0

simX,measurements,s0,xref,yref,psiref=casadi_integrator_fun(y_meas,con,DT)
print(simX)


