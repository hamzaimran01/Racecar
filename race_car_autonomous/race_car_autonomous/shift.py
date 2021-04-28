from linear_optimization_Example import *

def shift(x0,u,f):
    st=x0
    con=u[1,:].T
    st=integrator(st,con)
    x0=st.full()
    u0 = np.append(u[2:np.size(u,1),:],u(np.size(u,1)))