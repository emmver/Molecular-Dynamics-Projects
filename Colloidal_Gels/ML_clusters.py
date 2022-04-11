'''
 Attempt to get cluster size by first tagging each particle with a structural vector based on Machine Learning approaches in 
 combination with Steinhardt Vectors for recognizing crystalline structure. 

 Biblio: 
 1) https://aiche.onlinelibrary.wiley.com/doi/abs/10.1002/aic.16157
 2) https://arxiv.org/abs/1802.03426
 3) https://aip.scitation.org/doi/10.1063/1.2977970
 4) https://freud.readthedocs.io/en/latest/gettingstarted/examples/examples/Using%20Machine%20Learning%20for%20Structural%20Identification.html

'''


####################################### Libraries & Functions ########################################
import numpy as np
import sys
import warnings
from matplotlib import pyplot as plt
import matplotlib as mpl
import os 
import glob 
import freud
import pandas as pd
warnings.filterwarnings("ignore")
############################# Plotting style #############################################################
mpl.rcParams['figure.dpi'] = 400
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

def read_lammps_frame(filename):
    print('Reading:',filename)
    N = int(np.genfromtxt(filename, skip_header=3, max_rows=1))
    box_data = np.genfromtxt(filename, skip_header=5, max_rows=3)
    box = freud.box.Box.from_box(box_data[:, 1] - box_data[:, 0])
    data = np.genfromtxt(filename, skip_header=9, invalid_raise=False)
    data = data[~np.isnan(data).all(axis=1)].reshape(-1, N, 9) ### To properly reshape need number of particles N and number of columns in data file...in this case 9 
    data[..., :3] -= box.L / 2 ### shift to match the freud coordinate system
    positions=data[0][:,1:4]
    return positions, box, N


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

os.chdir(sys.argv[1])

files=glob.glob('100percent_*')
cycl=[0,100,300,400]
names=['sc','bcc','fcc','hcp']
structure_features = {}
positions,box,N=read_lammps_frame(files[0])
for idx,file in enumerate(files):
    data,box,N=read_lammps_frame(file)
    positions=data[:,1:4]
    for name in names:
        structure_features[name] = get_features(box, positions, name)

