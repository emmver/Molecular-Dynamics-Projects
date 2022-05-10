#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 13 12:50:10 2019

@author: manos
"""

import numpy as np
#import string
from utilities import read_lammpstrj,plotrg,calc_rg
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 
import sys , os 


    

_file=sys.argv[1]
f=open (_file , "r")
#runs=float(input("Number of runs?"))
step=float(input("Framestep?"))
jj=0
for i in range (int(1E10)):
    #print (i)
    jj+=1

    out, time, mybox, nparticles=read_lammpstrj(f)
    if i > (time*1000/step):
           break
f.close()
_file=sys.argv[1]
f=open (_file , "r")
frames=jj
print ("Frames available",frames)
module_path = os.path.dirname(os.path.realpath(_file));print(module_path)

rgfile=input("Name of Rg file?")
rg_path = os.path.join(module_path, rgfile);print (rg_path)
rg=np.zeros((frames,7))
with open (rg_path,'w') as rgfile:
            rgfile.write("#Time,Rg/n")
for i in range (frames):
    print ("Frame is" ,i)
    jj=i

    out, time, mybox, nparticles=read_lammpstrj(f)
    if i > (time*1000/step):
           break
    rgsave,rgx,rgy,rgz=calc_rg (out, mybox, nparticles,time)
    rg[i,0]=time
    rg[i,1]=rgsave
    rg[i,2]=rgx
    rg[i,3]=rgy
    rg[i,4]=rgz
    rg[i,5]=np.sqrt(rgz/rgx)
    rg[i,6]=np.sqrt(rgy/rgx)
rg=rg[:frames-1,:frames-1]
print("Break at" ,jj)
#rg = np.reshape(np.asarray(rg),(int(len(rg)),2))
with open (rg_path,'w') as rgfile:
            np.savetxt(rgfile,rg,delimiter=" ",header="time Rg lx ly lz l3/l1 l2/l1")
plotrg(rg_path)
