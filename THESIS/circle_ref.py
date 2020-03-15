import numpy as np
import matplotlib.pyplot as plt

R = 0.7                     # radius of circle
circumference = 2*np.pi*R   # circumference of circle
ncircle = 50                # data points per circle
rounds = 1                  # rounds
yoffset=-0.9                # offset of circle in y direction
xoffset=0.1                 # offset in x direction

sref=np.linspace(0,circumference*rounds,ncircle*rounds)    # s values
psiref=np.linspace(0,rounds*2*np.pi,ncircle*rounds)               # psi values
xref=R*np.sin(psiref) +xoffset                                     # xref values
yref=-R*np.cos(psiref)+yoffset                                      # yref values
kapparef=1/R*np.ones(round(ncircle*rounds))

np.savetxt("circle_07.txt",np.transpose([sref,xref,yref,psiref,kapparef]))
plt.plot(xref,yref)
plt.show()

