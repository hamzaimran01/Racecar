
import sys
sys.path.append("../THESIS")


from casadi import *
import numpy as np

def track_load(filename):
	track_file=np.loadtxt("spline_tasks/model/tracks/LMS_Track.txt")
	s0=track_file[:,0]
	xref=track_file[:,1]
	yref=track_file[:,2]
	psiref=track_file[:,3]
	kapparef=track_file[:,4]

	#length=len(s0)
	#pathlength=s0[-1]


	#length=len(s0)
	#pathlength=s0[-1]
	#copy loop to beginning and end
	#s0=np.append(s0,[s0[length-1] + s0[1:length]])
	#kapparef=np.append(kapparef,kapparef[1:length])
	#s0 = np.append([-s0[length-2] + s0[length-81:length-2]],s0)
	#kapparef = np.append(kapparef[length-80:length-1],kapparef)
	#kapparef_s=interpolant('kapparef_s','bspline',[s0],kapparef)

	#xref=np.append(xref,[xref[length-1] +xref[1:length]])
	#yref=np.append(yref,[yref[length-1]+yref[1:length]])
	#psiref=np.append(psiref,[psiref[length-1]+psiref[1:length]])
	return(s0,kapparef,xref,yref,psiref)


	



