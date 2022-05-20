import numpy as np 
from matplotlib import pyplot as plt 
##from numba import njit
from mpl_toolkits.mplot3d import Axes3D
#from utilities import read_lammpstrj
import sys
import os.path
from numba import njit
import os       
import glob

import matplotlib as mpl
#############  Plots Setup #####################
mpl.rcParams['figure.dpi'] = 400
plt.rcParams["font.family"] = "Ubuntu"
plt.style.use('C:/Users/Toumba/Documents/plotstyle.mplstyle')
plt.rcParams['axes.linewidth'] = 4
plt.rcParams['xtick.major.size'] = 8 
plt.rcParams['ytick.major.size'] = 8 
plt.rcParams['xtick.labelsize']=20
plt.rcParams['ytick.labelsize']=20
plt.rcParams['xtick.minor.visible']=True
plt.rcParams['ytick.minor.visible']=True
plt.rcParams['xtick.minor.size'] = 5 
plt.rcParams['ytick.minor.size'] = 5
plt.rcParams['xtick.minor.width'] = 1.5
plt.rcParams['ytick.minor.width'] = 1.5
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['ytick.major.width'] = 2
plt.rcParams['xtick.major.pad']='8'
plt.rcParams['ytick.major.pad']='8'
#########################################################################################################

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

@njit(fastmath=True)
def rg_tens (r):
    
    N = r.shape[0]
    r2=np.zeros((3,3))
    for i in range (N):
        for j in range (N):
            r2[0,0]+=(r[i,0]-r[j,0])*(r[i,0]-r[j,0])#diag
            r2[1,0]+=(r[i,1]-r[j,0])*(r[i,1]-r[j,0])#ndiag
            r2[2,0]+=(r[i,2]-r[j,0])*(r[i,2]-r[j,0])#ndiag
            r2[0,1]+=(r[i,0]-r[j,1])*(r[i,0]-r[j,1])#ndiag
            r2[0,2]+=(r[i,0]-r[j,2])*(r[i,0]-r[j,2])#ndiag
            r2[1,1]+=(r[i,1]-r[j,1])*(r[i,1]-r[j,1])#diag
            r2[2,1]+=(r[i,2]-r[j,1])*(r[i,2]-r[j,1])#ndiag
            r2[1,2]+=(r[i,1]-r[j,2])*(r[i,1]-r[j,2])#ndiag
            r2[2,2]+=(r[i,2]-r[j,2])*(r[i,2]-r[j,2])# diag
    r2=r2/(2*N**2)
    #rg2=r2[0,0]**2+r2[1,1]**2+r2[2,2]**2
    return r2


def how_many_frames(f):
    count=0
    while 1:
        try:
            out, time, mybox, nparticles=read_lammpstrj(f)
            count+=1
        except:
            break
    return count




_file=sys.argv[1]
f=open (_file , "r")


out, time1, mybox, nparticles=read_lammpstrj(f)
out, time2, mybox, nparticles=read_lammpstrj(f)
framestep=time2-time1
print("Framestep:",framestep)

dt=0.001
f.close()

f=open (_file , "r")
path=os.path.realpath(f.name)
rpath=os.path.dirname(path)
print("Path:",path)
print("Rpath:",rpath)
directory=rpath+"/"+"rg_tens"
import os
if not os.path.exists(directory):
    os.makedirs(directory)
print("Created Directory",directory)
frames=how_many_frames(f)

f.close()
f=open (_file , "r")

print("Frames:",frames)
nbb=int(sys.argv[2])
gen=int(sys.argv[3])
count=0
for i in range(frames):
    out, time, mybox, nparticles=read_lammpstrj(f)
    unwrap = out[:,2:5]+out[:,5:]*mybox;
    rgtens=rg_tens(unwrap)
    np.savetxt(directory+"/rg_tens_frame_%d.dat"%(count),rgtens)
    print("Frame:",count)
    count+=1


tens_files=glob.glob(directory+"/rg_tens_*")
directory=rpath+"/"+"tens_res"
import os
if not os.path.exists(directory):
    os.makedirs(directory)
print("Created Directory",directory)
#print(len(tens_files))
rg2=[]
asph=[]
time=[]
kappa2=[]
tc=0
acyl=[]
for i in tens_files:
    data=np.genfromtxt(i)
    #print(data.shape)
    temp_rg2=data[0,0]+data[1,1]+data[2,2]
    temp_asph=(3/2)*data[2,2]-temp_rg2/2
    temp_k=(3/2)*((data[0,0]**2+data[1,1]**2+data[2,2]**2)/((data[0,0]+data[1,1]+data[2,2])**2))-(1/2)
    temp_c=data[1,1]-data[0,0]
    rg2.append(temp_rg2)
    asph.append(temp_asph)
    kappa2.append(temp_k)
    acyl.append(temp_c)
    tc+=1
    time.append(tc*framestep*dt)

time=np.array(time)
rg2=np.array(rg2)
rg=rg2**0.5
asph=np.array(asph)
kappa2=np.array(kappa2)
acyl=np.array(acyl)
np.savetxt(directory+'/tens_results.dat',np.stack((time,rg2,rg,asph,kappa2,acyl),axis=-1))



plt.plot(time,rg2)
plt.ylabel('R${_g}{^2}[σ^2]$')
plt.xlabel('t/$τ_0$')
plt.xscale('log')
plt.savefig(directory+"/rg2_time.jpg",dpi=300,bbox_inches='tight')
plt.clf()
plt.plot(time,rg)
plt.ylabel('R${_g}[σ]$')
plt.xlabel('t/$τ_0$')
plt.xscale('log')
plt.savefig(directory+"/rg_time.jpg",dpi=300,bbox_inches='tight')
plt.clf()
plt.plot(time,asph)
plt.ylabel('b[σ$^2$]')
plt.xlabel('t/$τ_0$')
plt.xscale('log')
plt.savefig(directory+"/asph_time.jpg",dpi=300,bbox_inches='tight')
plt.clf()
plt.plot(time,kappa2)
plt.ylabel('κ$^2$')
plt.xlabel('t/$τ_0$')
plt.xscale('log')
plt.savefig(directory+"/anis_time.jpg",dpi=300,bbox_inches='tight')
plt.clf()
plt.plot(time,kappa2)
plt.ylabel('c$[σ^2]$')
plt.xlabel('t/$τ_0$')
plt.xscale('log')
plt.savefig(directory+"/acyl_time.jpg",dpi=300,bbox_inches='tight')
plt.clf()








