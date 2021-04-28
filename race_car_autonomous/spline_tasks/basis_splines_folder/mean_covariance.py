import sys
sys.path.append("../THESIS")
import numpy as np
from casadi import *
import matplotlib.pyplot as plt


#s_difference= np.loadtxt("final_data/0.8_training_data/estimated_data/s_difference.txt")
#n_difference= np.loadtxt("final_data/0.8_training_data/estimated_data/n_difference.txt")
#alpha_difference= np.loadtxt("final_data/0.8_training_data/estimated_data/alpha_difference.txt")
#v_difference= np.loadtxt("final_data/0.8_training_data/estimated_data/v_difference.txt")

s_difference= np.loadtxt("final_data/1.0_data/estimated_data/s_difference.txt")
n_difference= np.loadtxt("final_data/1.0_data/estimated_data/n_difference.txt")
alpha_difference= np.loadtxt("final_data/1.0_data/estimated_data/alpha_difference.txt")
v_difference= np.loadtxt("final_data/1.0_data/estimated_data/v_difference.txt")
N_MHE=1565

s_mean=np.mean(s_difference,axis=0)
#print(s_mean)
cov_s=np.var(s_difference)
#print(cov_s)

n_mean=np.mean(n_difference,axis=0)
#print(n_mean)
cov_n=np.var(n_difference)
#print(cov_n)


alpha_mean=np.mean(alpha_difference,axis=0)
#print(alpha_mean)
cov_alpha=np.var(alpha_difference)
#print(cov_alpha)


v_mean=np.mean(v_difference,axis=0)
#print(v_mean)
cov_v=np.var(v_difference)
#print(cov_v)

print(s_difference.shape)

l = []
for _ in range(1):
    l.append(np.random.uniform(0, 1, 100))
    
mean_s =s_difference  
standard_dev_s = 0.0013638 
N_MHE=len(s_difference)
plt.xlabel(r'$s_{ref}[m]$ track')
plt.ylabel(r'$ s $[m]')
plt.plot(mean_s,label=r'difference: $s$ measurement and $s$ estimate ')
plt.fill_between(range(N_MHE),mean_s-standard_dev_s,mean_s+standard_dev_s,color='purple',alpha=0.3,label='standard deviation')
plt.legend()
plt.show()


mean_n =n_difference  
standard_dev_n = 0.0007127
N_MHE=len(n_difference)

plt.xlabel(r'$s_{ref}[m]$ track')
plt.ylabel(r'$n$[m]')
plt.plot(mean_n,label=r'difference: $n$ measurement and $n$ estimate ')
plt.fill_between(range(N_MHE),mean_n-standard_dev_n,mean_n+standard_dev_n,color='purple',alpha=0.3,label='standard deviation')
plt.legend()
plt.show()


mean_alpha =alpha_difference  
standard_dev_alpha = 0.04
N_MHE=len(alpha_difference)
plt.xlabel(r'$s_{ref}[m]$ track')
plt.ylabel(r'$\alpha [\degree]$ ')
plt.plot(mean_alpha,label=r'difference: $\alpha$ measurement and $\alpha$ estimate ')
plt.fill_between(range(N_MHE),mean_alpha-standard_dev_alpha,mean_alpha+standard_dev_alpha,color='purple',alpha=0.3,label='standard deviation')
plt.legend()
plt.show()


mean_v =v_difference  
standard_dev_v = 0.00767
N_MHE=len(v_difference)
plt.xlabel(r'$s_{ref}[m]$ track')
plt.ylabel(r'$v$$[\frac{m}{s}]$ ')
plt.plot(mean_v,label=r'difference: $v$ measurement and $v$ estimate ')
plt.fill_between(range(N_MHE),mean_v-standard_dev_v,mean_v+standard_dev_v,color='purple',alpha=0.3,label='standard deviation')
plt.legend()
plt.show()