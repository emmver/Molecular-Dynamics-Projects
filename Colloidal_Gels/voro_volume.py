import warnings
import os
import sys
import freud
import matplotlib.cm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib as mpl
from scipy.stats import gaussian_kde,iqr
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import normalize
from sklearn.svm import SVC
from umap import UMAP
from unitcell import  UnitCell as unitcell
import glob
warnings.filterwarnings("ignore")

def get_features(box, positions, structure):
    voro = freud.locality.Voronoi()
    voro.compute(system=(box, positions))
    nlist = voro.nlist.copy()
    nlist.filter(nlist.weights > 0.1)
    features = {}
    for l in [4, 6, 8, 10, 12]:
        ql = freud.order.Steinhardt(l=l, weighted=True)
        ql.compute(system=(box, positions), neighbors=nlist)
        features[f"q{l}"] = ql.particle_order

    return features



def write_to_lammps(filename,box,n,positions):
    f=open(filename,'w')
    f.write('ITEM: TIMESTEP\n')
    f.write('0\n')
    f.write('ITEM: NUMBER OF ATOMS\n')
    f.write('%d\n'%n)
    f.write('ITEM: BOX BOUNDS xx yy zz pp pp pp\n')
    f.write('%.2f %.2f\n'%(-box,box))
    f.write('%.2f %.2f\n'%(-box,box))
    f.write('%.2f %.2f\n'%(-box,box))
    f.write('ITEM: ATOMS id x y z\n')
    for i in range(positions[:,0].size):
        f.write('%d %.3f %.3f %.3f\n'%(i,positions[i,0],positions[i,1],positions[i,2]))

    f.close()

def fred_diac(data_interp_rg):
    rgs=data_interp_rg
    rgs=np.asarray(rgs)
    rg_iqr=iqr(rgs)#;print(rg_iqr)
    bw=2*rg_iqr/rgs.size**(1/3)#;print(bw)
    nbins=(rgs.max()-rgs.min())/bw#
    return int(nbins)


def read_lammps_frame(filename):
    print('Reading:',filename)
    N = int(np.genfromtxt(filename, skip_header=3, max_rows=1))
    box_data = np.genfromtxt(filename, skip_header=5, max_rows=3)
    box = freud.box.Box.from_box(box_data[:, 1] - box_data[:, 0])
    data = np.genfromtxt(filename, skip_header=9, invalid_raise=False)
    data = data[~np.isnan(data).all(axis=1)].reshape(-1, N, 9) ### To properly reshape need number of particles N and number of columns in data file...in this case 9 
    data[..., :3] -= box.L / 2 ### shift to match the freud coordinate system
    #positions=data[0]
    return data[0], box, N
    
############################# Plotting style #############################################################
mpl.rcParams['figure.dpi'] = 120
plt.rcParams["font.family"] = "Ubuntu"
#plt.style.use('C:/Users/Toumba/Documents/plotstyle.mplstyle')
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

np. set_printoptions(threshold=np. inf)

########################################################################################################
################### Calculate voro volumes ############################
os.chdir(sys.argv[1])
files=glob.glob('100percent_*')
print(files)
cycl=[0,100,300,400]
p_vol=(4/3)*np.pi * 0.5**3
for i in range(0,4):
    for idx,file in enumerate(files):
        data,box,N=read_lammps_frame(file)
        if i==0:
            print('zero')
        elif i==1: 
            indices=np.where(data[:,9]<8)
            data=data[indices]
        elif i==2:
            indices=np.where(data[:,9]>=8)
            data=data[indices]
        elif i==3: 
            indices=np.where(np.where(data[:,9]>10)
            data=data[indices]
            
        positions=data[:,1:4]
        print(positions.shape)
        voro = freud.locality.Voronoi()
        voro.compute(system=(box, positions))
        vol_frac=1/(voro.volumes/p_vol)

        
        print(vol_frac.shape)
        plt.hist(vol_frac,bins=fred_diac(vol_frac),density=True,alpha=0.7,label='Cycle: %d'%cycl[idx]);

    plt.legend(frameon=False)
    for lh in plt.legend().legendHandles:
                lh.set_alpha(1)
    plt.xlabel(r'$\phi$')
    plt.ylabel(r'$P(\phi)$')
    plt.savefig('./histo_vol_frac_%d.jpg'%i,dpi=300,bbox_inches='tight')
    #plt.show();
    plt.clf();