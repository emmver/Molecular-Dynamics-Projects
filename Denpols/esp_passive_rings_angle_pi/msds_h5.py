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
#### MSD Libs #######

def _autocorr_fft(x):
    """Compute the autocorrelation of 1d signal using FFT."""
    N = x.shape[0]
    # 2*N because of zero-padding to compute non-cyclic correlation
    f = np.fft.fft(x, n=2*N)
    power_spectrum = f * f.conjugate()
    result = np.fft.ifft(power_spectrum)
    # the autocorrelation in the usual convention B
    result = (result[:N]).real
    return result

@njit 
def _msd_fft_compute_s1(r):
    """Compute the non-FFT part of the MSD in the FFT method."""
    N = r.shape[0]
    # Compute the squared distances
    D = np.square(r).sum(axis=1)
    # Apply the recursive algorithm
    s1 = np.zeros(N)
    s1[0] = 2*D.sum()
    for n in range(1, N):
        s1[n] = s1[n-1] - D[n-1] - D[N-n]
    return s1

def self_fft(r):
    """Compute the self MSD from a single particle trajectory using FFT.

    Based on the algorithm outlined in Section 4.2 of the following paper:
    https://doi.org/10.1051/sfn/201112010
    """
    N = r.shape[0]
    # Compute the non-FFT part
    s1 = _msd_fft_compute_s1(r)
    # Compute the FFT-part separately over each position component
    s2 = np.zeros_like(s1)
    for i in range(r.shape[1]):
        s2 += _autocorr_fft(r[:, i])
    return (s1 - 2*s2) / np.arange(N, 0, -1)
###################################################################


############### Calculating Quantities MSDs #########################

#### Create pos array for Fast MSD Calculation ######
directory='msds_results'
import os
if not os.path.exists(directory):
    os.makedirs(directory)

h5_dir='h5_sets'
store_dir=h5_dir
framestep=float(sys.argv[1])
tot_t=float(sys.argv[2])
frames=int(np.floor(tot_t/framestep))
nchains=float(sys.argv[3])
Dpol=float(sys.argv[4])
npart=int(nchains*Dpol)
bonds=npart
angles=npart
dens=0.85
box=(npart/dens)**(1/3)
file=sys.argv[5]



############### g1 calc ##################
hf = h5py.File(store_dir+'/pos_%s'%(file[2:-4]), 'r')
n1 = hf.get('%s'%(file[2:-4]))
pos=np.array(n1)
hf.close()
mean_pos=pos.mean(axis=1).reshape(-1, 1, 3)
pos-=mean_pos


time = np.zeros(frames)
count=0
for tt in range(frames):
    time[tt]=count
    count+=1
    
msds = np.zeros((npart, frames))
for i in range(pos.shape[1]):
    start=timer()
    msds[i] = self_fft(pos[:,i])
    end=timer()
    g="Particle %d Wall time %.4f seconds"%(i, end-start)
    sys.stdout.write("\r" + str(g))
    #sys.stdout.write('\n')
    sys.stdout.flush()
#       np.savetxt('./'+file[:-4]+'_g1_all.dat',(np.stack((msds),axis=-1)))
np.savetxt(directory+'/'+file[2:-4]+'_g1.dat',(np.stack((time[1:],msds.mean(axis=0)[1:]),axis=-1)))
hf = h5py.File(h5_dir+'/msds_%s'%(file[2:-4]), 'w')
hf.create_dataset('%s'%(file[2:-4]),data=msds)
hf.close()
######### g3 Calc ######################


nch=int(nchains)
pos_coms=np.zeros((frames,nch,3))

hf = h5py.File(store_dir+'/pos_coms_%s'%(file[2:-4]), 'r')
n1 = hf.get('%s'%(file[2:-4]))
pos_coms=np.array(n1)
hf.close()
mean_pos=pos_coms.mean(axis=1).reshape(-1, 1, 3)
pos_coms-=mean_pos
#pos=np.zeros((len(files),npart,3))
g3s = np.zeros((nch,frames))
for i in range(pos_coms.shape[1]):
    start=timer()
    g3s[i] = self_fft(pos_coms[:,i])
    end=timer()
    g="Particle %d Wall time %.4f seconds"%(i, end-start)
    #sys.stdout.write("\r" + str(g))
    #sys.stdout.write('\n')
    #sys.stdout.flush()
np.savetxt(directory+'/'+file[2:-4]+'_g3.dat',(np.stack((time[1:],g3s.mean(axis=0)[1:]),axis=-1)))
hf = h5py.File(h5_dir+'/coms_msds_%s'%(file[2:-4]), 'w')
hf.create_dataset('%s'%(file[2:-4]),data=g3s)
hf.close()