# -*- coding: utf-8 -*-
"""
Created on Wed May 08 15:43:48 2019

@author: Lyra
"""

import numpy as np 
from matplotlib import pyplot as plt 
##from numba import njit
from mpl_toolkits.mplot3d import Axes3D
from utilities import read_lammpstrj
import sys
import os.path
_file=sys.argv[1]
f=open (_file , "r")
path=os.path.realpath(f.name)
rpath=os.path.dirname(path)
runs=float (input('number of runs?'))
framestep=float(input('Framestep?'))
filename=input("Filename")
skipframes=int(input("Frames to Skip?"))
#runs=10**9
#framestep=10**6
frames= int (runs/framestep)
nbb=int(input("Backbone Atoms?"))
gen=int(input("Generation"))
timesRN=np.zeros((frames,2))
with open(r'timeSim', 'w') as outfile:
    outfile.write("TIME")
avgendtoend=0
endtoend=0
endtoend0=0

endtoends=[]; times = []
for i in range (skipframes):
    out, time,mybox,nparticles=read_lammpstrj(f)

print ('After skipframes',time)
    
for i in range (skipframes,frames):
   # print (i)
   

    out, time, mybox, nparticles=read_lammpstrj(f)
    #print (time*1000/framestep)
    print("In Loop...", time)
    if i > np.ceil((time/framestep)):
        print(i,time/framestep)
        break
    if gen==0:
        if i==0:
            time0 = time
        
        unwrap = out[:,2:5]+out[:,5:]*mybox;
        endtoend=unwrap[(nbb-1)]-unwrap[0];
        #avgendtoend=avgendtoend+endtoend
        #print (np.linalg.norm(endtoend))
        endtoends.extend([endtoend[0], endtoend[1], endtoend[2]])
        times.extend([time-time0])

   
   
   
    else:
   
        if i==skipframes:
            time0 = time
        
        unwrap = out[:,2:5]+out[:,5:]*mybox;
        endtoend=unwrap[(nbb-1)*(1+2**(gen+1))]-unwrap[0];
        #avgendtoend=avgendtoend+endtoend
        #print (np.linalg.norm(endtoend))
        endtoends.extend([endtoend[0], endtoend[1], endtoend[2]])
        times.extend([time-time0])
        
f.close

eends = np.reshape(np.asarray(endtoends),(int(len(endtoends)/3),3))
_M =eends.shape[0]

acorr = np.zeros((_M,2), dtype = float)
acorr[:,0] = np.asarray(times)

for i in range(_M): ##tau
    for j in range(_M-i):  ##t0
        acorr[i,1] += np.dot(eends[i+j,:],eends[j,:])/np.dot(eends[j,:],eends[j,:])

for i in range(_M):
    acorr[i,1] /= float(_M -i )


fig = plt.figure()
ax = fig.add_subplot(111)


ax.plot(acorr[:,0], acorr[:,1],"ro--",linewidth=2)
ax.set_xscale("log")

plt.show()
#with open(r'correl', 'w') as outfile:
#    outfile.write("%.4f \n" % avgendtoend)

#with open(r'correl', 'a') as outfile:
outfiles=os.path.join(rpath,filename+".txt")
np.savetxt(outfiles, acorr, fmt='%.5e',delimiter=' ')    
