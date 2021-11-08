import numpy as np 
import os
import glob
from timeit import default_timer as timer
import sys
from numba import njit

@njit(fastmath=True)
def coms_times(pos,pos_coms,DPsm,nsmall,frames):
    for i in range (frames):
        data_small=pos[i]
        for j in range (nsmall):
            #print("Calculating small chain",j)
            workdata=data_small[j*DPsm:(j+1)*DPsm]
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

def read_h (files):
    #### Read Header File, getting number of particles, box size, timestep #######
    f=open(files[0],'r')
    trash = f.readline()
    trash = f.readline()
    elems = str.split(trash.strip()," ") 
    #  print (str.isdigit(elems[0]))
    if str.isdigit(elems[0])==True:
        #print ("yay")
        time1 = float(elems[0])/1
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
    #### Read Header File of second file, getting number of particles, box size, timestep #######
    f=open(files[1],'r')
    trash = f.readline()
    trash = f.readline()
    elems = str.split(trash.strip()," ") 
    #  print (str.isdigit(elems[0]))
    if str.isdigit(elems[0])==True:
        #print ("yay")
        time2 = float(elems[0])/1
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
    framestep=int(time2-time1)
    return framestep, npart, mybox

######### Loaded Libraries + Functions #####################33



#### Build pos array #######
wdir=sys.argv[1]
wfile="test_pos*"
filename=wdir+wfile
files=glob.glob(filename)
files.sort(key=os.path.getmtime)
framestep,npart,box=read_h(files)
data=np.genfromtxt(files[0],skip_header=9)
data=data[np.argsort(data[:,0])]
data=data[np.where(data[:,1]==2)]
pos=np.zeros((len(files),data[:,0].size,3))
#print("Comp:",100*data[:,0].size/npart)
#print('\n')
DPsm=96;DPlarge=250
data=np.loadtxt(wdir+'test_pos0.dat',skiprows=9)
data=sorted(data,key= lambda part: part[0] ) ## Sorting data with increasing particle number
data=np.asarray(data)#### Converts data into an array, because it is initially loaded and sorted as a list
data_small=data[np.where(data[:,1]==1)]#### Finds where particle type is 1 (small) ###########
data_large=data[np.where(data[:,1]==2)]#### Finds where particle type is 2 (large) ###########
comp_small=data_small[:,0].size/(data_small[:,0].size+data_large[:,0].size)
nsmall=int(npart*comp_small/DPsm)
nlarge=int((npart-npart*comp_small)/DPlarge)
print('Small:',nsmall,'Large:',nlarge)
print("Composition small:",100*comp_small)

count=0
start=timer()
for i in files:  
    start_n=timer()
    data=np.genfromtxt(i,skip_header=9)
    data=data[np.argsort(data[:,0])]
    data=data[np.where(data[:,1]==2)]
    temp=data[:,2:5]+data[:,5:8]*box
    pos[count]=temp
    count+=1
    end=timer()
    tpf=end-start_n
    files_left=len(files)-count
    t_left=tpf*files_left
    if (end-start)%2<5e-2 and (end-start)>1:
            g="%.2f percent in %.3f seconds ETA:%.3f seconds for folder %s"%(100*count/len(files),end-start,t_left,sys.argv[1])
            sys.stdout.write("\r" + str(g))
            sys.stdout.write('\n')
mean_pos=pos.mean(axis=1).reshape(-1, 1, 3)
pos-=mean_pos


sys.stdout.write('\n')
sys.stdout.flush()
#### Calculate MSDs #######
time = np.zeros(len(files))
count=0
for tt in range(len(files)):
    time[tt]=count*framestep*0.01
    count+=1

    
msds = np.zeros((int(nlarge*DPlarge), len(files)))
for i in range(pos.shape[1]):
    start=timer()
    msds[i] = self_fft(pos[:,i])
    end=timer()
    g="Particle %d Wall time %.4f seconds"%(i, end-start)
    sys.stdout.write("\r" + str(g))
    sys.stdout.write('\n')
np.savetxt(wdir+'large_g1.dat',(np.stack((time[1:],msds.mean(axis=0)[1:]),axis=-1)))



#### COMs MSD #####

pos_coms=np.zeros((len(files),nlarge,3))
pos_coms=coms_times(pos,pos_coms,DPlarge,nlarge,len(files))
#pos=np.zeros((len(files),npart,3))
g3s = np.zeros((nlarge, len(files)))
for i in range(pos_coms.shape[1]):
    start=timer()
    g3s[i] = self_fft(pos_coms[:,i])
    end=timer()
    g="Particle %d Wall time %.4f seconds"%(i, end-start)
    sys.stdout.write("\r" + str(g))
    sys.stdout.write('\n')
    #sys.stdout.flush()
np.savetxt(wdir+'large_msd_coms.dat',(np.stack((time[1:],g3s.mean(axis=0)[1:]),axis=-1)))
