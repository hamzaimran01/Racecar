import sys
from casadi import *


N_MHE=980

nr_points=256
nr_sub_points=50
sys.path.append("../THESIS")
P_meas=MX.sym('P_meas',6*(N_MHE)+6)#(measurements )+state for starting 




args={}
p=np.zeros( (6*(N_MHE) + 6) )
args['p']=p

lbg=np.zeros((6*(N_MHE+1)))
ubg=np.zeros((6*(N_MHE+1)))
#ubg[0:4*(N):1]=0
#lbg[0:4*(N):1]=0
args['lbg']=lbg
args['ubg']=ubg

ubx=np.zeros((6*(N_MHE+1)))
lbx=np.zeros((6*(N_MHE+1)))
#s
lbx[0:6*(N_MHE+1):6]=0
ubx[0:6*(N_MHE+1):6]=+inf
#n
lbx[1:6*(N_MHE+1):6]=-0.12
ubx[1:6*(N_MHE+1):6]=+0.12
#alpha
lbx[2:6*(N_MHE+1):6]=-inf
ubx[2:6*(N_MHE+1):6]=+inf
#v
lbx[3:6*(N_MHE+1):6]=-1
ubx[3:6*(N_MHE+1):6]=+1
#deltaD

lbx[4:6*(N_MHE+1):6]=-1
ubx[4:6*(N_MHE+1):6]=+1

#delta delta
lbx[5:6*(N_MHE+1):6]=-0.4
ubx[5:6*(N_MHE+1):6]=+0.4

#lbx[4*(N_MHE+1):4*(N_MHE+1)+3*nr_points+2:1]=-inf
#ubx[4*(N_MHE+1):4*(N_MHE+1)+3*nr_points+2:1]=+inf

lbx=np.append(lbx,-inf*np.ones((3*(nr_points+2))))
ubx=np.append(ubx,inf*np.ones((3*(nr_points+2))))


args['lbx']=lbx
args['ubx']=ubx


X0_1 = np.zeros(((N_MHE+1,6)))
#print(N_MHE)
#print(X0_1.shape)

DIS_1=np.zeros(((nr_points+2,3)))

args['x0']=np.append(reshape((X0_1.T),6*(N_MHE+1),1),reshape((DIS_1.T),3*(nr_points +2),1))