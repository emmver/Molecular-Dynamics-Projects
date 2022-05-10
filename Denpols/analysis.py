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

### For analysis of structure factor ####
_sk = ctypes.CDLL('/mnt/c/Users/Toumba/Documents/GitHub/emmver/Denpols/struc_fact/libsk.so')
_sk.calc_sk.argtypes = (ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_int, ctypes.c_int)

nk = int(np.log(10./0.001)/np.log(1.05)); print ("How many k",nk)
kgrid, knorm, nkbins= generate_kgrid_rand(nk, 0.001, 1.05)

Sskarr = np.zeros((nkbins),dtype = np.double)
Ssktrans = np.zeros((nkbins),dtype = np.double)

Sskbb = np.zeros((nkbins),dtype = np.double)

dr = 0.1;Lmax = 200.;nbins=int(Lmax/dr)
c_r = np.zeros((nbins), dtype = float)
cr_bb = np.zeros((nbins), dtype = float)

jj = -1
totrep = 0
mygen = int(sys.argv[3])
Lbackbone = int(sys.argv[4])
mybox = data.cell[0,0]
#########################################################################

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


print(os.getcwd())
rgsq=[]
lambda_1=[]
lambda_2=[]
lambda_3=[]
jj = -1
totrep = 0
for frame_index in tqdm(range(0,pipeline.source.num_frames),desc='G%d_N%d'%(mygen,Lbackbone)):
    #print("frame:",frame_index)
    data = pipeline.source.compute(frame_index)
    pos=data.particles.positions[:]
    r2=rg_tens(pos)
   # rgsq.append(r2[0,0]+r2[1,1]+r2[2,2])
    np.savetxt(to_store+'/rg_tens_frame_%d.dat'%frame_index,r2)

    analyze_sk(data,Sskarr,mybox,Sskbb,kgrid,nkbins,_sk,Lbackbone)
    totrep+=1
norm_sk(totrep,Sskarr,Sskbb,Ssktrans,knorm,nkbins,mygen,Lbackbone,pos[:][:,0].size,to_store)


#gyr_tens=np.vstack([np.array(rgsq),np.array(lambda_1),np.array(lambda_2),np.array(lambda_3)])
#np.savetxt(to_store+'/gyration_tensor.txt',gyr_tens)
#plot_shape(pipeline.source.num_frames,gyr_tens)



