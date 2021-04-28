
import numpy as np
import re

####control separation D,delta
#f_read_control="rawdata/controls_1.0.txt"
f_read_control="rawdata/hamza_controls_lap3.txt"
#f_read_control="rawdata/hamza_controls_lap3.txt"
f_write_control_der="garbage/dD_ddelta_control.txt"
f_write_control="garbage/D_delta_control.txt"
f_write_t_control="garbage/time_combined_control.txt"






###1.extracting seconds and nanoseconds from the control file to a  single column###########

time_ns = []
linenum = 0
time_string = "nsecs:"# for seconds detection as string
time_string_sec = "secs:"#for nano seconds detection as a string
pattern = re.compile(time_string, re.IGNORECASE)  # Compile a case-insensitive regex
pattern_s = re.compile(time_string_sec, re.IGNORECASE)
time=[]
with open(f_read_control,'r') as fr,open(f_write_t_control,'w') as fwt:    
	for line in fr:
		if pattern.search(line) != None or pattern_s.search(line) != None:    # If a match is found copy the seconds and nanoseconds to time_combined_control file
			line.rstrip()			
			line=line.replace(time_string,"")
			line=line.replace(time_string_sec,"")
			fwt.write(line)
			time.append(line)

fwt.close()


#seconds and nano seconds intilizations

time = np.loadtxt("garbage/time_combined_control.txt")

#half of it is seconds and the other half is nanoseconds
len_time = len(time)
len_time_s=int(0.5*len_time)
len_time_ns=int(0.5*len_time)

sec = np.zeros((len_time_s,1))
nsec = np.zeros((len_time_ns,1))

sec = time[0:len_time:2]
nsec = time[1:len_time:2]
np.savetxt("final_data/sec_controls",sec)
np.savetxt("final_data/nsec_controls",nsec)
		

with open(f_read_control) as fr, open(f_write_control_der,"w") as fw:
	for line in fr:
		if  (line.startswith("dD:") or line.startswith("ddelta:")   ):
			line = line.replace("dD: ","")
			line = line.replace("ddelta: ","")
			fw.write(line)

		if  (line.startswith("D:") or line.startswith("delta:")   ):
			line = line.replace("D: ","")
			line = line.replace("delta: ","")
			fw.write(line)

fw.close()

			
#print("has a column of all the controls and dercontrols")

#with open(f_read_control) as fr, open(f_write_control,"w") as fw:
#    for line in fr:
#            
#            if  (line.startswith("D:") or line.startswith("delta:")   ):
#                line = line.replace("D: ","")
#                line = line.replace("delta: ","")
#                #line=line.replace("ddD: ","")
#                #line=line.replace("dddelta: ","")
#                
#                fw.write(line)
        




#data of control and control states in one column loaded which was closed above
data_control = np.loadtxt("garbage/dD_ddelta_control.txt",dtype=float)


#separate columns of state controls and controls and divide them into D,delta,derD and derDelta
len_control = len(data_control)
len_control_D = int(0.25*len_control)
len_control_delta = int(0.25*len_control)
len_control_der_D = int(0.25*len_control)
len_control_der_delta = int(0.25*len_control)




control_data_D = np.zeros((len_control_D,1))
control_data_delta = np.zeros((len_control_delta,1))
control_data_der_D = np.zeros((len_control_der_D,1))
control_data_der_delta = np.zeros((len_control_der_delta,1))

control_data_D = np.array(data_control[0:len_control:4])
control_data_delta = data_control[1:len_control:4]
control_data_der_D = np.array(data_control[2:len_control:4])
control_data_der_delta = np.array(data_control[3:len_control:4])



#final form of the controls and state cotnrols 2 columns


data_control=np.concatenate((control_data_D.reshape(-1,1),control_data_delta.reshape(-1,1),control_data_der_D.reshape(-1,1),control_data_der_delta.reshape(-1,1),sec.reshape(-1,1),nsec.reshape(-1,1)),axis=1)



#np.savetxt("final_data/control_der_final",control_data_der)
   





#########state separation s,n,alpha, v
#f_read_s="rawdata/states_1.0.txt"
f_read_s = "rawdata/hamza_state_measurements_lap3.txt"
f_write_s = "garbage/states.txt"
f_write_time_states="garbage/states_time.txt"

delete_list_s = ["s:","n:","alpha:","v:"]



time_ns = []
linenum = 0
time_string = "nsecs:"# for seconds detection as string
time_string_sec = "secs:"#for nano seconds detection as a string
pattern = re.compile(time_string, re.IGNORECASE)  # Compile a case-insensitive regex
pattern_s = re.compile(time_string_sec, re.IGNORECASE)

with open(f_read_s,'r') as fr_s,open(f_write_time_states,'w') as fw_ts:    
	for line in fr_s:
		if pattern.search(line) != None or pattern_s.search(line) != None:    # If a match is found 
			line=line.replace(time_string,"")
			line=line.replace(time_string_sec,"")
			fw_ts.write(line)
fw_ts.close()

states_time=np.loadtxt("garbage/states_time.txt",dtype=float)
len_states_time=len(states_time)
print(len_states_time)
###################################check this

states_time_sec_l=int(0.5*len_states_time)
states_time_nsec_l=int(0.5*len_states_time)

states_time_sec=np.zeros((states_time_sec_l,1))
states_time_nsec=np.zeros((states_time_nsec_l,1))


with open(f_read_s) as fr_s, open(f_write_s,"w") as fw_c:
    for line in fr_s:
        if(line.startswith("s:") or line.startswith("n:") or line.startswith("alpha:")or line.startswith("v:")):
            for word in delete_list_s:
                if word in line:
                    line=line.replace(word,"")
                    fw_c.write(line)






#all states in one column
data_states_column = np.loadtxt("garbage/states.txt",dtype=float)
len_states = len(data_states_column)
#print(len_states)


#separate columns of states
states_time_sec=np.array(states_time[0:len_states:2])
states_time_nsec=np.array(states_time[1:len_states:2])    
states_data_s = np.array(data_states_column[0:len_states:4])
states_data_n = np.array(data_states_column[1:len_states:4])
states_data_alpha = np.array(data_states_column[2:len_states:4])
states_data_v = np.array(data_states_column[3:len_states:4])
print(states_data_s.shape)
print(states_data_n.shape)
print(states_data_alpha.shape)
print(states_data_v.shape)

states_data = np.concatenate((states_data_s.reshape(-1,1),states_data_n.reshape(-1,1),states_data_alpha.reshape(-1,1),states_data_v.reshape(-1,1),states_time_sec.reshape(-1,1),states_time_nsec.reshape(-1,1)),axis=1)


#states=np.savetxt("final_data/states_final",states_data)    



#states_data
#control_data_der

#######################################starting of the data comparison part

len_states=states_data.shape[0]
len_controls=data_control.shape[0]
states_final=np.zeros((len_states,8))

for i in range(len_controls):
	nsecs_control=data_control[i,5]
	for k in range(len_states):
		nsecs_states=states_data[k,5]
		if(nsecs_control==nsecs_states):
			#print("match found")
			#print(nsecs_states)
			states_final[i,0:4]=states_data[k,0:4]
			states_final[i,4:6]=data_control[i,0:2]
			states_final[i,6:8]=data_control[i,4:6]

		else:
			pass
			#print("match not found")
			#print(nsecs_states)

np.savetxt("final_data/states_final_alligned",states_final)
np.savetxt("final_data/control_final_alligned",data_control[:,2:6])


#for i in range(len_controls):
#	if(control_data_der[i,1]==states_final[i,5]):
#		print("good till now ..")


#	elif(control_data_der[i,1]!=states_final[i,5]):
#		print("fuck")
#		break





