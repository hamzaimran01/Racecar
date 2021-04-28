

import sys
sys.path.append("../THESIS")

import numpy as np
import matplotlib.pyplot as plt
from common_parameters.common import N_main
from model.track_load import *
cost_function_disturbed=np.loadtxt('configurations/cost_function_mpc_module_disturbed')
cost_function_ndisturbed=np.loadtxt('configurations/cost_function_mpc_module_ndisturbed')
cost_function_offline=np.loadtxt('configurations/cost_function_mpc_module_disturbed_offline')
#print(cost_function_offline.size)
length=len(cost_function_disturbed)
result=[]
count_p=0
count_n=0
print(length)
#print(cost_function_offline)
plt.xlabel('sref')
plt.ylabel('cost')
plt.plot(s0[0:N_main+1],cost_function_disturbed,label='disturbed model mpc')
#plt.plot(s0[0:N_main-1],cost_function_ndisturbed,label='not disturbed model mpc')
plt.plot(s0[0:N_main],cost_function_offline,label='offline mpc') 
plt.legend()   
plt.grid()
plt.savefig('plots_save/costfun-result_varying_disturbedvsoffline.pdf')
plt.show()