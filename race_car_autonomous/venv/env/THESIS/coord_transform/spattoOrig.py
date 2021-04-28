
import numpy as np

def SpattoOrig(s_plot,n_plot,alpha_plot,velocity_plot,s0,xref,yref,psiref):
    track_length=s0[-1]
    si=s_plot%track_length
    idxmindist=findClosestS(si,s0)
    idxmindist2=findSecondClosestS(si,s0,idxmindist)
    #print(idxmindist)
    #print(idxmindist)  
    t=(si-s0[idxmindist])/(s0[idxmindist2]-s0[idxmindist])
    x0=(1-t)*xref[idxmindist]+t*xref[idxmindist2]  
    y0=(1-t)*yref[idxmindist]+t*yref[idxmindist2]
    psi0=(1-t)*psiref[idxmindist]+t*psiref[idxmindist2]
    x=x0-n_plot*np.sin(psi0)  
    y=y0+n_plot*np.cos(psi0)
    psi=psi0+alpha_plot  
    v=velocity_plot
    return x,y,psi,v
 
def findClosestS(si,sref):
     idxmindist=[]
     # create iter (because si can be array or float)
     siiter = (si,) if not isinstance(si.tolist(), list) else si
     for siel in siiter:
          di=abs(siel-sref)
          idxmindist.append(np.argmin(di))
     idxmindist = np.where(idxmindist==sref.size,1,idxmindist)
     idxmindist = np.where(idxmindist<1,sref.size-1,idxmindist)
     
     return idxmindist.astype(int)
  
def dist2D(x1,x2,y1,y2): 
    return np.sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2))  

def findClosestNeighbour(x,y,xref,yref,idxmindist):
    distBefore=dist2D(x,xref[idxmindist-1],y,yref[idxmindist-1])
    distAfter=dist2D(x,xref[idxmindist+1],y,yref[idxmindist+1])
    if(distBefore<distAfter):
        idxmindist2=idxmindist-1
    else:
        idxmindist2=idxmindist+1
    if(idxmindist2<0):
        idxmindist2=xref.size-1  
    elif(idxmindist==xref.size):
        idxmindist2=0 
    return idxmindist2                 



def findSecondClosestS(si,sref,idxmindist): 
     d1=abs(si-sref[idxmindist-1])            
     d2=abs(si-sref[(idxmindist+1)%sref.size])
     idxmindist2 = np.where(d1>d2,idxmindist+1,idxmindist-1)  
     idxmindist2 = np.where(idxmindist2==sref.size,0,idxmindist2)    
     idxmindist2 = np.where(idxmindist2<0,sref.size-1,idxmindist2)
     return idxmindist2  