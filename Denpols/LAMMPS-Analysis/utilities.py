#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 11:13:17 2019

@author: manos
"""
import numpy as np
from numba import njit,jit
from matplotlib import pyplot as plt 








def plotrg(inputfile) :
    
    data0=np.loadtxt(inputfile)




    b=data0
    a=np.mean(b[:,1],axis=0)
    c=np.std(b[:,1],axis=0)
    f1=plt.figure()
    ax1=f1.add_subplot(111)



    ax1.plot(b[:,0],b[:,1])
    ax1.set_xlabel("time")
    ax1.set_ylabel("$R_g$")
    ax1.set_title( "Mean= {}, \n STD={}".format(a,c))
    plt.show()
    f1.savefig("my_rg.png",dpi=360)
    
    return




#Calculate Center of Mass
#   @njit(nopython=True)
def calc_rg (out, mybox, npart,time):
    comx=np.sum(out[:,2]+(out[:,5]*mybox))
    comy=np.sum(out[:,3]+(out[:,6]*mybox))
    comz=np.sum(out[:,4]+(out[:,7]*mybox))
    com = np.array([comx,comy,comz])/(npart)
    #Calculate Rg squared
    rgsqx=np.sum(np.dot((out[:,2]+(out[:,5]*mybox)-com[0]),(out[:,2]+(out[:,5]*mybox)-com[0])))
    rgsqy=np.sum(np.dot((out[:,3]+(out[:,6]*mybox)-com[1]),(out[:,3]+(out[:,6]*mybox)-com[1])))
    rgsqz=np.sum(np.dot((out[:,4]+(out[:,7]*mybox)-com[2]),(out[:,4]+(out[:,7]*mybox)-com[2])))
    rgsq=(rgsqx+rgsqy+rgsqz)/npart
    rg=rgsq**0.5;print ("Rg is ",rg)
    #rgsave=np.array([time,rg])
    #np.append(rgsave,[time,rg])
    #rgsave=rgsave.reshape(1,2)
    
    return rg,rgsqx,rgsqy,rgsqz

@njit (error_model='numpy',fastmath=True)
def calc_denpol_rdf (unwrap,dr,bins,N,nmon,gen,Linv,rcut,mybox):
    hist1=np.zeros((bins,nmon))

    for i in range (N):
            for j in range (nmon):
             if i%(1+2**gen) == 0 and (1+2**gen)*j == i:
                 #print (i, j, 17*j)
                 continue
             dx=unwrap[i,0]- unwrap[(1+2**gen)*j,0]#comx
             dy=unwrap[i,1]- unwrap[(1+2**gen)*j,1]#comy
             dz=unwrap[i,2]- unwrap[(1+2**gen)*j,2]#comz
             #print "||" ,dx, "||", dy, "||" ,dz
             #print "==============================================================================="
             dx=dx-mybox*np.rint(dx*Linv)
             dy=dy-mybox*np.rint(dy*Linv)
             dz=dz-mybox*np.rint(dz*Linv)
             rabs=np.sqrt(dx*dx+dy*dy+dz*dz)
             if rabs < rcut:
                 idx=int(rabs/dr);#print (idx)
                 hist1[idx,j]=hist1[idx,j]+1
       

        
    return hist1
@njit
def calc_linear_rdf (unwrap,dr,bins,nmon,Linv,rcut,mybox):
    hist1=np.zeros((bins,nmon))

    for i in range (nmon):
            for j in range (nmon):
             if j == i:
                 #print (i, j, 17*j)
                 continue
             dx=unwrap[i,0]- unwrap[j,0]#comx
             dy=unwrap[i,1]- unwrap[j,1]#comy
             dz=unwrap[i,2]- unwrap[j,2]#comz
             #print "||" ,dx, "||", dy, "||" ,dz
             #print "==============================================================================="
             dx=dx-mybox*np.rint(dx*Linv)
             dy=dy-mybox*np.rint(dy*Linv)
             dz=dz-mybox*np.rint(dz*Linv)
             rabs=np.sqrt(dx*dx+dy*dy+dz*dz)
             if rabs < rcut:
                 idx=int(rabs/dr);#print (idx)
                 hist1[idx,j]=hist1[idx,j]+1
       

        
    return hist1


@jit (error_model='numpy',fastmath=True)
def calc_molrdf(out,unwrap,bins,time,npart,mybox,particlespermonomer,runs,framestep,nmon,rho,updates):
    step=time*1000
    #laststep=step
    print (step)
    hist=np.zeros(bins)
    dr=mybox/bins ; #print ("dr is =",dr)
    rdf_r=dr*(np.arange(bins)+0.5)
    nideal=4*np.pi*rho*dr*rdf_r**2
    #updates=0
    #bins=512
    #print (bins, time ,npart,mybox,particlespermonomer,runs,framestep,nmon,rho)
    if (step==0) :
        hist=np.zeros(bins)
        dr=mybox/bins ;
        rdf_r=dr*(np.arange(bins)+0.5);
        nideal=4*np.pi*rho*dr*rdf_r**2
        #updates=0; Define updates as zero at the script and function will return updates value and keep it in the loop
        print ("nay")
    elif  (step==10**6) :
        hist=np.zeros(bins); 
        dr=mybox/bins; print (dr)
        rdf_r=dr*(np.arange(bins)+0.5); 
        nideal=4*np.pi*rho*dr*rdf_r**2; 
        updates=0;
        #print (updates)
        gr=np.zeros(bins)
    elif  (step==1.1*10**7) :
        hist=np.zeros(bins)
        dr=mybox/bins ; #print ("dr is =",dr)
        rdf_r=dr*(np.arange(bins)+0.5)
        nideal=4*np.pi*rho*dr*rdf_r**2
        #updates=0
    else:
        N=out.shape[0]
        Linv=1/mybox ;#print (Linv)
        rcut=mybox
        hist1=np.zeros((bins,nmon))
        for i in range (N):
            for j in range (nmon):
             if i% particlespermonomer == 0 and particlespermonomer*j == i:
                 #print (i, j, 17*j)
                 continue
             dx=unwrap[i,0]- unwrap[particlespermonomer*j,0]#comx
             dy=unwrap[i,1]- unwrap[particlespermonomer*j,1]#comy
             dz=unwrap[i,2]- unwrap[particlespermonomer*j,2]#comz
             #print "||" ,dx, "||", dy, "||" ,dz
             #print "==============================================================================="
             dx=dx-mybox*np.rint(dx*Linv)
             dy=dy-mybox*np.rint(dy*Linv)
             dz=dz-mybox*np.rint(dz*Linv)
             rabs=np.sqrt(dx*dx+dy*dy+dz*dz)
             #print (dr)
             if rabs < rcut:
                 idx=int(rabs/dr)
                 hist1[idx,j]=hist1[idx,j]+1
        
#        hist[:]+= np.sum(hist1[:,:],axis=1)/nmon
#        updates=updates +1
#        print (updates)
#        gr=(hist)/(nideal*updates) 
 
    a=np.arange(bins)
    save=np.zeros((a.size,2))
    save[:,0]=rdf_r
    save[:,1]=gr
    return save,step

def average_frames_rdf(bins,runs,framestep,npart,rdfavgfile_path):
    out = list()    
    #with open (rdfavgfile_path,'w') as grlastfile:
        #grlastfile.write(("#Averaged over Frames g(r) \n"))
    with open (rdfavgfile_path,'w') as grlastfile:
        grlastfile.write("#Step is:None\n")
    
    f = open(r"rdf.txt", 'r')
    next(f)
    frames=runs/framestep
    
    
    i = 0
    k=0
    avgrdf=np.zeros((int(bins),1))
    toavg=np.zeros((int(bins),int(frames)))
    trash = f.readline();#print (trash)
    trash = f.readline(); #print (trash)
    
    while k < (int(frames)):
        k=k+1
    
        
       # print ("k is ",k)
        for i in range(int(bins)):
            #print (i)  
            line = f.readline()[1:]
            #elems = str.split(line.strip()," ");
            elems = str.split(line.strip()," ")
            #print i,"||",k
            for el in elems:
                    out.extend([float(el)])
    
                    if len(elems) > 0 : elsize = len(elems)
        trash=f.readline();#print (trash)
    
    
    
        out = np.reshape(np.asarray(out), (int(len(out)/elsize),elsize))
        rdf_r=np.reshape(out[:,0],(int(bins),1))
        toavg[:,k-1]+=out[:,1]
        #print (toavg[:,k-1])
        
        #print (out.shape)
    
        out=list()
        #toavg[:,k]+=out
    for j in range (int(bins)):
        avgrdf[j,:]+=np.mean(toavg[j,:])
    rdf=np.zeros((int(bins),2))
    rdf[:,0]=np.reshape(rdf_r,(int(bins),))
    rdf[:,1]=np.reshape(avgrdf,(int(bins),))
        
    with open (rdfavgfile_path,'a') as grfile:
      #grfile.write("#Step is:{}\n".format(step))
      np.savetxt(grfile, rdf,delimiter=" ")
        
    print (4*np.pi*1*np.trapz(avgrdf[:,0]*rdf_r[:,0]**2,rdf_r[:,0]),"||",npart) 
    f.close()
   # g.close()
    return



def read_lammpstrj(f):
    out1=list()
    out=[]
    out=np.asarray(out)
    trash = f.readline()
    trash = f.readline()
    elems = str.split(trash.strip()," ") 
  #  print (str.isdigit(elems[0]))
    if str.isdigit(elems[0])==True:
        #print ("yay")
        time = float(elems[0])/1
        trash = f.readline() 
        trash = f.readline()
        elems = str.split(trash," ") 
        npart = int(elems[0]);
        trash = f.readline()
        trash = f.readline()
        trash = f.readline()
        trash = f.readline()
        elems = str.split(trash," ")
        mybox = float(elems[1])-float(elems[0])
        trash = f.readline()
    
        for i in range(npart):
            
            line = f.readline()
            elems = str.split(line.strip()," ")
            for el in elems:
                    out1.extend([float(el)])
    
                    if len(elems) > 0 : elsize = len(elems);
        out1 = np.reshape(np.asarray(out1), (int(len(out1)/elsize),elsize))
        #print(out1)
        out= np.copy(out1)
        out = out[out[:,0].argsort()]
    
        #print (type(out))
        #print (out.shape)
        nparticles = out.shape[0]
    #    
    #    with open('out.txt', 'ab') as outfile:
    #        np.savetxt(outfile, out, fmt='%4.1f')
    
    
        out1 = list()
        #f.close()
        return out,time, mybox, nparticles

    else:
        print ("nay")
        out1 = list()
        nparticles=0
        mybox=0
        time=1E-16
        #npart
        return out,time, mybox, nparticles

    


def q_grid(numberofindices,mybox):
    dx=2*np.pi/mybox
    q=list()
    for  i in range (numberofindices):
        for j in range (numberofindices):
            for k in range (numberofindices):
                q.extend([i*dx,j*dx,k*dx])
                 
    q=np.reshape(np.asarray(q),((numberofindices)**3,3))
    return q

def read_box(f):
    trash = f.readline()
    trash = f.readline()
    elems = str.split(trash.strip()," ") 
    time = float(elems[0])/1000.0
    trash = f.readline() 
    trash = f.readline()
    elems = str.split(trash," ") 
    npart = int(elems[0]);
    trash = f.readline()
    trash = f.readline()
    trash = f.readline()
    trash = f.readline()
    elems = str.split(trash," ")
    mybox = float(elems[1])-float(elems[0])
    trash = f.readline()
    return mybox, npart,time
