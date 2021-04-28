

import sys
sys.path.insert(0,"../THESIS")
import numpy as np
from casadi import *
from model.track_load import *
from scipy.interpolate import BSpline
import matplotlib.pyplot as plt
from coord_transform.spattoOrig import SpattoOrig
import scipy.interpolate as si
sref_points=[]





import scipy.interpolate as intrp

#x = np.linspace(0., 1., 100)

sref=np.array(256)
for x in s0:
    sref=np.append(sref,x)
       
    if(x>8.71):
        break

print(sref.shape)
#this defines us the total range from 0 to 1
t=np.linspace(0 ,1, 2560)
#these are the discrete points of the track defined by 256 knots
t_knots=np.linspace(0,8.71,257)

# starting and end defined fot t_knots
print(t_knots.shape)
numpy_t_knots=np.concatenate(([0,0,0],t_knots,[256,256,256]))

#this defines that 100 subpoints points and 258 total points
x_basis = np.zeros((t.shape[0], len(t_knots)+2))


for i in range(259):
    
    x_basis[:,i]=intrp.BSpline(numpy_t_knots, (np.arange(len(t_knots)+2)==i).astype(float), 3, extrapolate=False)(t)
#print(x_basis.shape)

plt.plot(t,x_basis)
plt.title('In SciPy')
plt.show()