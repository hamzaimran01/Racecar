import sys
sys.path.append("../THESIS")
from casadi import *

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
P_3=MX.sym('P',3,nr_points+2)


""" from multiprocessing import Pool

from casadi import *

N = 250
rs = np.linspace(1,3,N)

x = SX.sym('x')
y = SX.sym('y')

v = vertcat(x,y)
f = (1-x)**2+(y-x**2)**2
g = x**2+y**2
nlp = {'x': v, 'f': f, 'g': g}

# Create IPOPT solver object
solver = nlpsol('solver', 'ipopt', nlp)

def optimize(r):
  res = solver(x0=[2.5,3.0],      # solution guess
               lbx=-inf,          # lower bound on x
               ubx=inf,           # upper bound on x
               lbg=-inf,          # lower bound on g
               ubg=r)         # upper bound on g
   
  
  return float(res["f"])

b = Pool(2)
fsol = b.map(optimize,rs)  """