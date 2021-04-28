import numpy as np
from scipy import misc
from math import *
from scipy import optimize
from casadi import *

m = 0.043
C1=0.5
C2=15.5
Cm1=0.58
Cm2=0.15
Cr0=0.029
Cr2=0.006
kapparef=1.4285714285714286
s=0
n=0
alpha=0
v=1.0



from scipy.optimize import fsolve, root



from scipy.optimize import fsolve, root


def fsolve_function(arguments):
    derD = arguments[0]
    derDelta = arguments[1]
    
    """ out=[(v*np.cos(alpha+C1*derDelta))/(1-kapparef*n)]
    out.append(v*sin(alpha+C1*derDelta))
    out.append((v*C2*derDelta)-kapparef*((v*np.cos(alpha+C1*derDelta))/(1-kapparef*n))) 
    out.append(((Cm1-Cm2*v)*derD-Cr2*v*v-Cr0*tanh(5*v))/m*cos(C1*derDelta)) """
    out=[((Cm1-Cm2*v)*derD-Cr2*v*v-Cr0*tanh(5*v))/m*cos(C1*derDelta)]
    out.append((v*C2*derDelta)-kapparef*((v*np.cos(alpha+C1*derDelta))/(1-kapparef*n)))
    out.append(v*sin(alpha+C1*derDelta))
    out.append((v*np.cos(alpha+C1*derDelta))/(1-kapparef*n))
    return out

initialGuess = [1, 1]
result = root(fsolve_function, initialGuess, method='lm')
print(result.x)

roots=[0.08138923, 0.09206817]
print(fsolve_function(roots))

 #v=1,[0.08138923, 0.09206817] [0.9989406186101207, 0.0460178279673152, 3.547213989207876e-08, 4.591374406234403e-08]
 #v=1.5,[0.11971826 0.09206817]

#print(fsolve_function(roots))

