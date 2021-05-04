from casadi import *
from readDataFcn import getTrack
from obstacle_handler import *


def bycicle_model(track='LMS_Track.txt'):
    model_name = 'Spatialbycicle_model'

    # load track parameters
    [s0,xref,yref,psiref,kapparef]=getTrack(track)
    length=len(s0)
    pathlength=s0[-1]
    # copy loop to beginning and end
    s0=np.append(s0,[s0[length-1] + s0[1:length]])
    kapparef=np.append(kapparef,kapparef[1:length])
    s0 = np.append([-s0[length-2] + s0[length-81:length-2]],s0)
    kapparef = np.append(kapparef[length-80:length-1],kapparef)

    # compute spline interpolations
    kapparef_s=interpolant('kapparef_s','bspline',[s0],kapparef)

    ## Race car parameters
    m = 0.043
    C1=0.5
    C2=15.5
    Cm1=0.58
    Cm2=0.15
    Cr0=0.029
    Cr2=0.006
    mu_x = 0.8 # coefficient for maximum longitudinal force Fx = mu_x * F_gravity
    mu_y = 0.5 # coefficient for maximum lateral force Fy = mu_y * F_gravity

    ## CasADi Model
    # set up states & controls
    s   = MX.sym('s')
    n   = MX.sym('n')
    alpha = MX.sym('alpha')
    v   = MX.sym('v')
    D = MX.sym('D')
    delta=MX.sym('delta')
    x = vertcat(s,n,alpha,v,D,delta)

    # controls
    derD = MX.sym('derD')
    derDelta=MX.sym('derDelta')
    u = vertcat(derD,derDelta)

    # xdot
    sdot= MX.sym('sdot')
    ndot= MX.sym('ndot')
    alphadot= MX.sym('alphadot')
    vdot = MX.sym('vdot')
    Ddot= MX.sym('Ddot')
    deltadot = MX.sym('deltadot')
    xdot = vertcat(sdot,ndot,alphadot,vdot,Ddot,deltadot)

    # algebraic variables
    z = vertcat([])

    # boundary splines
    nknots=120                                  # number of knots
    scon_l=pathlength/nknots*2                  # distance between knots
    scon_min=-scon_l*(nknots-1)/4               # starting knot (hast to be 0 with current impl)
    scon_max=scon_l*(nknots-1)*3/4              # end knot
    scon=np.linspace(scon_min,scon_max,nknots)  # knots of boundary splines (equidistant with 10cm)
    nl=createSpline(scon,scon_l)
    nr=createSpline(scon,scon_l)

    # parameters
    pnl = MX.sym('pnl',2*(len(scon)-1))
    pnr = MX.sym('pnr',2*(len(scon)-1))
    p=vertcat(pnl,pnr)

    # dynamics
    Fxd=(Cm1-Cm2*v)*D-Cr2*v*v-Cr0*tanh(5*v)# I tested the following formulation for better modelling of breaking -> 15/(v+7)**2*(-D)*(D<0)
    sdota=(v*np.cos(alpha+C1*delta))/(1-kapparef_s(s)*n)

    
    params=type('status',(object,),{})()
    params.C1=C1
    params.C2=C2
    params.Cm1=Cm1
    params.Cm2=Cm2
    params.Cr0=Cr0
    params.Cr2=Cr2
    model=type('status',(object,),{})()
    model.f_impl_expr = xdot-f_expl
    model.f_expl_expr = f_expl
    model.x = x
    model.xdot = xdot
    model.u = u
    model.z = z
    model.p = p
    model.name = model_name
    model.params=params

    #constraint on lateral force
    a_lat = C2*v*v*delta+Fxd*sin(C1*delta)/m
    a_long = Fxd / m

    # define constraints
    constraint=type('status',(object,),{})()
    constraint.alat=Function('a_lat',[x,u],[a_lat])
    constraint.scon_min=scon_min
    constraint.scon_max=scon_max
    constraint.pathlength=pathlength
    constraint.nl=Function('nl',[s,p],[nl(s,pnl)])
    constraint.nr=Function('nr',[s,p],[nr(s,pnr)])
    constraint.expr=vertcat(a_long, a_lat ,n-nl(s,pnl),nr(s,pnr)-n)
    constraint.x=x
    constraint.u=u
    constraint.p=p
    constraint.nc=constraint.expr.shape[0]
    constraint.name="con_bycicle"

    # Model bounds
    model.n_min = -0.12         # width of the track [m]
    model.n_max = 0.12          # width of the track [m]

    # state bounds
    model.throttle_min = -1.0
    model.throttle_max = 1.0

    model.delta_min = -0.40     # minimum steering angle [rad]
    model.delta_max = 0.40      # maximum steering angle [rad]

    # input bounds
    model.ddelta_min = -2.0     # minimum change rate of stering angle [rad/s]
    model.ddelta_max = 2.0      # maximum change rate of steering angle [rad/s]
    model.dthrottle_min = -10   #-10.0  # minimum throttle change rate
    model.dthrottle_max = 10    #10.0  # maximum throttle change rate

    # nonlinear constraintNOTE constraints (subject to change)
    constraint.alat_min = -4    # maximum lateral force [m/s^2]
    constraint.alat_max = 4     # maximum lateral force [m/s^1]

    constraint.along_min = -4   # maximum lateral force [m/s^2]
    constraint.along_max = 4    # maximum lateral force [m/s^2]
    # define initial obstacle from 0.5 to 1 which reduces track width from 0.1 to 0 on the left half plane
    plength=round(p.shape[0]/4)
    model.p0=np.concatenate([model.n_max*np.ones(plength),model.n_max*np.ones(plength),model.n_min*np.ones(plength),model.n_min*np.ones(plength)])
    # Define initial conditions
    model.x0 = np.array([-2,0,0,0,0,0])

    return model,constraint


def createSpline(svec,l):
    n = len(svec)               # number of /knots
    slow = svec[:len(svec)-1]   # lower knots
    shigh = svec[1:len(svec)]   # higher knots
    s = MX.sym('s')             # progress
    nl=MX.sym('nl',n-1)         # distance of lower knot
    nh=MX.sym('nh',n-1)         # distance of higher knot
    p = vertcat(nl,nh)          # parameter vector
    false_vec = np.zeros(n -1)  # false vector (zeros when not in spline region)
    # spline definitions
    y = if_else(logic_and(s >= slow, s<shigh) , -2*(nh-nl)/l**3*(s-slow)**3+3*(nh-nl)/l**2*(s-slow)**2+nl, false_vec)
    y = dot(transpose(y) , transpose(np.ones(n-1)))
    fspline = Function('spline', [s,p], [y]).expand() # spline casadi function
    return fspline

