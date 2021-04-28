import sys
sys.path.append("../THESIS")

import matplotlib.pyplot as plt
import numpy as np
import os
from common_parameters.common import N_main
from model.track_load import *
#noise=np.random.normal(0,0.1,N_main)
#length_p=len(N_main)
#x=np.arange(0,2*np.pi,0.04)
#x=np.arange(0,15*np.pi,0.16)
x=s0[0:N_main]
x=x[0:N_main]
print(x.size)
#y=0.5*np.sin(x)

cwd=os.getcwd()
#print(cwd)
y_plot=[]
N_main=500
print(N_main)

for i in range(N_main+1):
    if(i<5):
        y_plot=np.append(y_plot,0.2)
        
    if(i>5 and i<=10):
        y_plot=np.append(y_plot,0.05)
    
    if(i>10 and i<=20):
        y_plot=np.append(y_plot,-0.15)
    
    if(i>20 and i <=40):
        y_plot=np.append(y_plot,-0.10)
        
    if(i>40 and i <=60):
        y_plot=np.append(y_plot,-0.17)
   
    if(i>60 and i <=80):
        y_plot=np.append(y_plot,0.08)
    if(i>80 and i <=100):
        y_plot=np.append(y_plot,0.25)
    
    if(i>100 and i <=120):
        y_plot=np.append(y_plot,0.20)
   
    if(i>120 and i <=150):
        y_plot=np.append(y_plot,-0.14)
    
    if(i>150 and i <=160):
        y_plot=np.append(y_plot,-0.18)
   
    if(i>160 and i <=170):
        y_plot=np.append(y_plot,0.27)
   
    if(i>170 and i <=195):
        y_plot=np.append(y_plot,-0.12)
    if(i>195 and i <=205):
        y_plot=np.append(y_plot,-0.18)
    if(i>205 and i <=220):
        y_plot=np.append(y_plot,0.25)
   
    if(i>220 and i<240):
        y_plot=np.append(y_plot,0.1)
   
    if(i>240 and i<260):
        y_plot=np.append(y_plot,0.2)

    if(i>260 and i<280):
        y_plot=np.append(y_plot,-0.15)

    if(i>280 and i<300):
        y_plot=np.append(y_plot,-0.05)

    if(i>300 and i<320):
        y_plot=np.append(y_plot,0.12)
       
    if(i>320 and i<340):
        y_plot=np.append(y_plot,0.12)
    
    if(i>340 and i<350):
        y_plot=np.append(y_plot,-0.10)
    

    if(i>350 and i<370):
        y_plot=np.append(y_plot,-0.15)
    

    if(i>370 and i<400):
        y_plot=np.append(y_plot,0.06)


    
    if(i>400 and i<420):
        y_plot=np.append(y_plot,0.10)
        
    

    if(i>420 and i<450):
        y_plot=np.append(y_plot,-0.05)
    

    if(i>450 ):
        y_plot=np.append(y_plot,-0.24)
    

    if(i>450 and i<500 ):
        y_plot=np.append(y_plot,-0.24)
    

    if(i>500):
        y_plot=np.append(y_plot,0.1)
     
""" for i in range(N_main):

    if(i<125):
        y_plot=np.append(y_plot,-0.25)

    if(i>=125 and i<250):
        y_plot=np.append(y_plot,0.25)
    

    if(i>=250 and i<375):
        y_plot=np.append(y_plot,-0.25)

    if(i>=375 and i<500):
        y_plot=np.append(y_plot,0.25)
    
 """    


    #y_plot=np.append(y_plot,0.25*np.sin(x[i]))
    #y_plot=np.append(y_plot,0.25*np.sin(x[i]))
    #y_plot=np.append(y_plot,0.25*np.sin(x[i]))
    #y_plot=np.append(y_plot,0.25*np.sin(x[i]))
    #y_plot=np.append(y_plot,0.25*np.sin(x[i]))
    

print("N-main value is")

plt.xlabel('s ref')
plt.ylabel('disturbance magnitude')
plt.step(s0[0:(N_main)],y_plot[0:N_main],label='disturbance')
#plt.step(s0[0:(N_main)],y_plot[0:N_main],label='disturbance')
plt.legend()
plt.grid()
plt.savefig('plots_save/disturbance_varying_sin.pdf')
plt.show()
np.savetxt('configurations/disturbance_varying_sin',y_plot)
print(y_plot)