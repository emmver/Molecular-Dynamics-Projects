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

lines=['-', '--', '-.', ':']*3
markers=[ 'o', 'v', '^', '<', '>', '1', '2', '3', '4', '8', 's', 'p', '*', 'h', 'H', '+', 'x', 'D', 'd', '|', '_', 'P', 'X', 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 'None', None, ' ', '']

#### Average and store  ########

DPs=[10,16,36,100,1000]
Ns=[1,2,4,6,11]
gen=3
gens=[3,4,5]
avgs_asph=np.zeros((3,6))
avgs_acyl=np.zeros((3,6))
avgs_anis=np.zeros((3,6))
in_dir="/tens_res"
filename="/tens_results.dat"

##time,rg2,rg,asph,kappa2,acyl
#### Creating Plotting array ######
for i in range(len(gens)):
	avgs_asph[i,0]=gens[i]
	avgs_acyl[i,0]=gens[i]
	avgs_anis[i,0]=gens[i]


	for j in range(1,len(Ns)+1):
		direc='G%d_%d'%(gens[i],Ns[j-1]) + in_dir
		data=np.genfromtxt(direc+filename)
		temp_asph=np.mean(data[:,3])
		avgs_asph[i,j]=temp_asph
		temp_anis=np.mean(data[:,4])
		avgs_anis[i,j]=temp_anis
		temp_acyl=np.mean(data[:,5])
		avgs_acyl[i,j]=temp_acyl

### Plotting ######
count=0
for i in range (1,avgs_anis.shape[1]):
	plt.plot(avgs_anis[:,0],avgs_anis[:,i],marker=markers[count],linestyle=lines[count],linewidth=4,label='DP=%d'%(DPs[i-1]))
	count+=1


plt.ylabel('κ',fontsize=35)
plt.xlabel('Generation',fontsize=35)
plt.xscale('log')
#plt.yscale('log')
plt.legend(loc=(1.0,0.5),frameon=False)
plt.savefig('./shape_plots/avg_anis.jpg',dpi=300,bbox_inches='tight')
np.savetxt('./shape_plots/avg_anis.dat',avgs_anis)
plt.clf()
count=0
for i in range (1,avgs_acyl.shape[1]):
	plt.plot(avgs_acyl[:,0],avgs_acyl[:,i],marker=markers[count],linestyle=lines[count],linewidth=4,label='DP=%d'%(DPs[i-1]))
	count+=1


plt.ylabel('c[σ$^2$]',fontsize=35)
plt.xlabel('Generation',fontsize=35)
plt.xscale('log')
#plt.yscale('log')
plt.legend(loc=(1.0,0.5),frameon=False)
plt.savefig('./shape_plots/avg_acyl.jpg',dpi=300,bbox_inches='tight')
np.savetxt('./shape_plots/avg_acyl.dat',avgs_acyl)
