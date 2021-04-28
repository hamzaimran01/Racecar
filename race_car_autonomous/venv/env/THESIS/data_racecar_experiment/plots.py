

import time
import sys,os
import numpy as np
import sys
sys.path.append("../THESIS")

from coord_transform.spattoOrig import *
from casadi import *
import matplotlib.pyplot as plt
from model.track_load import *

# TODO: make this plot function

x_simulation=np.loadtxt('x_simulation')
measurements=np.loadtxt('states_columns',dtype=float)
measurements=measurements[0:2000,0:4] # TODO(): just plot same amount of measurements

xtrack,ytrack,psiref,vref=SpattoOrig(x_simulation[:,0],x_simulation[:,1],x_simulation[:,2],x_simulation[:,3],s0,xref,yref,psiref)

xmeas,ymeas,psiref,vref=SpattoOrig(measurements[:,0],measurements[:,1],measurements[:,2],measurements[:,3],s0,xref,yref,psiref)
plt.plot(xmeas,ymeas,label='measurement')
plt.plot(xref,yref,label='track')
plt.plot(xtrack,ytrack,label='prediction')
plt.legend()
plt.show()





