import numpy as np # Library for multi-dimensional arrays and math operations on them 
import pandas as pd # Library for data manipulation and analysis
import os # Interfacing python with the operating system
import matplotlib.pyplot as plt #Plotting library
import seaborn as sns #Plotting library with some statistical tools 

from matplotlib.ticker import FuncFormatter
import glob
import matplotlib.pylab as pl
from ovito.io import import_file, export_file
from numba import njit
from tqdm import tqdm
from hoomd_utils import * ## custom functions
import sys
###### Data import ##########
global to_store, nk, kgrid,knorm, nkbins,Sskarr,Ssktrans,Sskbb,dr,Lmax,nbins,c_r,cr_bb,Nth,btheta,jj,totrep,mygen,Lbackbone,mybox,nparticles
os.chdir('/mnt/d/University/Simulations/Dendronized_polymer/Analysis/HOOMD_DATA/'+sys.argv[1])
print(os.getcwd())
filename = sys.argv[2]
to_store = './results_' +filename[2:-4]
if os.path.isdir(to_store)==False:
    os.mkdir(to_store)
pipeline=import_file(filename)
data = pipeline.source.compute()
################################
### My plotting style is inputted here #####
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 150
plt.rcParams["font.family"] = "Ubuntu"
plt.style.use('~/plotstyle.mplstyle')
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
#############################################
def avg_persistence():
    kk=0
    lp_proj_f=lp_proj[0]
    for i in range (1,pipeline.source.num_frames-1):
        lp_proj_f[:]=[sum(x) for x in zip(lp_proj_f,lp_proj[i])]
        kk+=1
    lp_proj_f[:]=[x/kk for x in lp_proj_f]
    outfiles=os.path.join(tostore,"persistency_g%d%d.txt"%(gen,Lbackbone))
    g=open(outfiles,"w")
    for i in range (len(bonds)):
        g.write("%.5f\t%.5f\n" % (bonds[i],lp_proj_f[i]))
    g.close()


#############################################

print(os.getcwd())
totrep = 0
lp_proj=[]
bonds=[]
mygen = int(sys.argv[3])
Lbackbone = int(sys.argv[4])
mybox = data.cell[0,0]
for i in range (0,Lbackbone-1):
	bonds.extend([(i+1)/Lbackbone])
for frame_index in tqdm(range(0,pipeline.source.num_frames),desc='G%d_N%d'%(mygen,Lbackbone)):
    data = pipeline.source.compute(frame_index)
    pos=data.particles.positions[:]
    tmp=[]

    for j in range (Nbb-1):
       # print(j*(1+2**(gen+1)))
        bondvect=unwrap[(j+1)*(1+2**(gen+1))]-unwrap[j*(1+2**(gen+1))]
        bondlength=np.linalg.norm(bondvect)
        calc=np.dot(bondvect,endtoend)/bondlength**2
        tmp.extend([calc])
    lp_proj.append(tmp)

    totrep+=1

avg_persistence()