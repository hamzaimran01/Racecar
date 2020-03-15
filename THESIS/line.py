import numpy as np
from casadi import *
import matplotlib.pyplot as plt

xoffset=0.1
yoffset=-0.9
sref=np.linspace(0,1,50)    # s values
psiref=np.zeros(50)
xref=np.linspace(0,1,50) +xoffset                                     # xref values
yref=-np.linspace(0,1,50) +yoffset                                      # yref values
kapparef=np.zeros(50)

np.savetxt("line.txt",np.transpose([sref,xref,yref,psiref,kapparef]))
plt.plot(xref,yref)
plt.show()
