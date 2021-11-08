from __main__ import *
import numpy as np 
#import MDAnalysis as mda
from matplotlib import pyplot as plt 
from numba import njit
import sys
from timeit import default_timer as timer
from scipy.stats import gaussian_kde,iqr
import glob
import ctypes
import h5py
def read_traj(f,npart):
    out=[]
    for i in range (npart):
            line = f.readline()
            elems=str.split(line.strip()," ")
            #print(elems)
            for el in elems: 
                out.extend([float(el)])
                if len(elems)>0: elsize=len(elems);
    out=np.reshape(np.asarray(out),(int(len(out)/elsize),elsize))
    f.readline()
    f.readline()
    return out


#################### MSDs Lib ############
@njit(fastmath=True)
def coms_times(pos,pos_coms,Dpolsm,nsmall,frames):
    for i in range (frames):
        data_small=pos[i]
        for j in range (nsmall):
            #print("Calculating small chain",j)
            workdata=data_small[j*Dpolsm:(j+1)*Dpolsm]
            comx=np.sum(workdata[:,0])
            comy=np.sum(workdata[:,1])
            comz=np.sum(workdata[:,2])
            com = np.array([comx,comy,comz])/(workdata[:,0].size)
            #com=CoM(workdata,workdata[:,0].size)
            pos_coms[i][j]=com
    return pos_coms



#### Create pos array and store as H5FD file ######
directory='h5_sets'
import os
if not os.path.exists(directory):
    os.makedirs(directory)


framestep=float(sys.argv[1])
tot_t=float(sys.argv[2])
frames=int(tot_t/framestep)
nchains=float(sys.argv[3])
Dpol=float(sys.argv[4])
npart=int(nchains*Dpol)
bonds=npart
angles=npart
dens=0.85
box=(npart/dens)**(1/3)
file=sys.argv[5]
f=open(file,'r')
pos=np.zeros((int(frames),npart,3))
for i in range (npart+bonds+angles+2):
    f.readline()
print(f.readline())
count_frames=0
start=timer()
while 1: 
    try: 
        g="Frame:%d"%(count_frames+1)
        #sys.stdout.write("\r" + str(g))
        #sys.stdout.flush()
        data=read_traj(f,npart)
        pos[count_frames]=data[:,1:]
        count_frames+=1
    except Exception as e:
        print("Except:",e)
        break
hf = h5py.File(directory+'/pos_%s'%(file[2:-4]), 'w')
hf.create_dataset('%s'%(file[2:-4]), data=pos)
hf.close()


nch=int(nchains)
pos_coms=np.zeros((frames,nch,3))
pos_coms=coms_times(pos,pos_coms,Dpol,nch,frames)
hf = h5py.File(directory+'/pos_coms_%s'%(file[2:-4]), 'w')
hf.create_dataset('%s'%(file[2:-4]), data=pos_coms)
hf.close()





