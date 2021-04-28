
import numpy as np

####control separation Ddelta
f_read_control="data_racecar_experiment/test_controls.txt"
f_write_control_der="data_racecar_experiment/dD_ddelta_control.txt"
f_write_control="data_racecar_experiment/D_delta_control.txt"

#lines=read_file.readlines()
#print(lines)

delete_list_c=["D:","delta:","ddelta:","dD:",]
#delete_list_c= [" "]


with open(f_read_control) as fr, open(f_write_control_der,"w") as fw:
    for line in fr:
        if  (line.startswith("dD:") or line.startswith("ddelta:")   ):
            line=line.replace("dD: ","")
            line=line.replace("ddelta: ","")
            #line=line.replace("ddD: ","")
            #line=line.replace("dddelta: ","")
            #print(line)
            fw.write(line)
with open(f_read_control) as fr, open(f_write_control,"w") as fw:
    for line in fr:
            
            if  (line.startswith("D:") or line.startswith("delta:")   ):
                line=line.replace("D: ","")
                line=line.replace("delta: ","")
                #line=line.replace("ddD: ","")
                #line=line.replace("dddelta: ","")
                
                fw.write(line)

fw.close()        
#########state separation snalpha v
f_read_s="data_racecar_experiment/test_states.txt"
f_write_s="data_racecar_experiment/states.txt"
delete_list_s=["s:","n:","alpha:","v:"]

with open(f_read_s) as fr_s, open(f_write_s,"w") as fw_c:
    for line in fr_s:
        if(line.startswith("s:") or line.startswith("n:") or line.startswith("alpha:")or line.startswith("v:")):
            for word in delete_list_s:
                if word in line:
                    line=line.replace(word,"")
                    fw_c.write(line)


states_data=np.zeros((2065,4))
control_data_der_D=np.zeros((1033,1))
control_data_der_delta=np.zeros((1033,1))
control_data_D=np.zeros((1033,1))
control_data_delta=np.zeros((1033,1))
control_data=np.zeros((1033,2))

data_states=np.loadtxt("data_racecar_experiment/states.txt",dtype=float)
data_control_der=np.loadtxt("data_racecar_experiment/dD_ddelta_control.txt",dtype=float)
data_control=np.loadtxt("data_racecar_experiment/D_delta_control.txt",dtype=float)
#print(data_control)

#print(data_states.size)
#print(data_states[2])
#print(states_data.shape)

lst_controls=list(range(1033))
for i in range(1,1033):
    
    states_data[i,0]=data_states[8*i]#4,8,12
    states_data[i,1]=data_states[8*i+1]#5,9,13
    states_data[i,2]=data_states[8*i+2]#6,10,14
    states_data[i,3]=data_states[8*i+3]#7,11,15
   
np.savetxt("data_racecar_experiment/states_columns",states_data)    
    
control_data_D=np.array(data_control[0:2066:2])
control_data_delta=data_control[1:2066:2]
    
control_data_der_D=np.array(data_control_der[0:2066:2])
control_data_der_delta=np.array(data_control_der[1:2066:2])

control_data_der=np.concatenate((control_data_D.reshape(-1,1),control_data_delta.reshape(-1,1)),axis=1)
control_data_der=np.concatenate((control_data_D.reshape(-1,1),control_data_delta.reshape(-1,1)),axis=1)


np.savetxt("data_racecar_experiment/control_columns_der",control_data_der)
np.savetxt("data_racecar_experiment/control_columns",control_data)    

print(control_data_D.shape)
print(control_data_delta.shape)