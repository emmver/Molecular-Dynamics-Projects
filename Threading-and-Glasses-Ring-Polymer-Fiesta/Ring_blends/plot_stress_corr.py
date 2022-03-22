import numpy as np 
#import MDAnalysis as mda
from matplotlib import pyplot as plt 
from numba import njit
import sys
from timeit import default_timer as timer
from scipy.stats import gaussian_kde,iqr
import glob
import ctypes
import matplotlib as mpl
import os
from numba import njit
import tarfile
from numba_progress import ProgressBar
from statsmodels.graphics.tsaplots import plot_acf
import statsmodels.api as sm
#############  Plots Setup #####################
mpl.rcParams['figure.dpi'] = 400
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
#########################################################################################################

data=np.genfromtxt('log.stress')[50:,:]
fig,(ax1,ax2,ax3,ax4)=plt.subplots(4,1,figsize=(10,9))
ax1.plot(data[:,1],data[:,2])
ax1.set_ylabel(r'$\sigma_{xy}$')
ax2.plot(data[:,1],data[:,3],'r-')
ax2.set_ylabel(r'$\sigma_{xz}$')
ax3.plot(data[:,1],data[:,4],'g-')
ax3.set_ylabel(r'$\sigma_{yz}$')
ax4.set_ylabel(r'$SFA(t)$')
ax4.set_xlabel(r'$t[\tau]$')
t_step=0.01
#fig.tight_layout()
#mov_avg=np.convolve(data[:,2], np.ones(100)/100, mode='same')
#print(data.shape, mov_avg.shape)
#plt.plot(data[:,1],mov_avg,'--');
#plt.show()
#fig.savefig('stress_compon.jpg',dpi=300,bbox_inches='tight')
#plt.clf()
#ax4=plt.axes()
_M=data[:,0].size
print(_M)
avg_list=[5,10,100,1000,2000,5000,10000,20000,50000,100000]
colors = plt.cm.cividis(np.linspace(0,1,len(avg_list)))
colors_xy=plt.cm.Blues(np.linspace(0,1,len(avg_list)+1)[1:])
colors_xz=plt.cm.Reds(np.linspace(0,1,len(avg_list)+1)[1:])
colors_yz=plt.cm.Greens(np.linspace(0,1,len(avg_list)+1)[1:])
for idx, i in enumerate(avg_list):
    print(i)
    acf_xy,ci=sm.tsa.acf(np.convolve(data[:,2], np.ones(int(i))/i, mode='same'),alpha=0.05,nlags=_M,fft=True)
    tacf=data[:,1]
    acf_xz,ci=sm.tsa.acf(np.convolve(data[:,3], np.ones(int(i))/i, mode='same'),alpha=0.05,nlags=_M,fft=True)
    acf_yz,ci=sm.tsa.acf(np.convolve(data[:,4], np.ones(int(i))/i, mode='same'),alpha=0.05,nlags=_M,fft=True)
    acf=(acf_xy+acf_xz+acf_yz)/3
    ax4.plot(tacf,np.abs(acf),label=r'$t_{avg}=%.2f$'%(i*t_step),color=colors[idx],marker='o')
    ax1.plot(data[:,1],np.convolve(data[:,2], np.ones(int(i))/i, mode='same'),color=colors_xy[idx])
    ax2.plot(data[:,1],np.convolve(data[:,3], np.ones(int(i))/i, mode='same'),color=colors_xz[idx])
    ax3.plot(data[:,1],np.convolve(data[:,4], np.ones(int(i))/i, mode='same'),color=colors_yz[idx])

plt.ylim(1e-3,2)
plt.xlim(10,1e5)
ax4.set_yscale('log')
ax4.set_xscale('log')
ax4.legend(loc=(1.0,0.2),frameon=False)
fig.savefig('stress_compon.jpg',dpi=300,bbox_inches='tight')
#acf,ci=sm.tsa.acf(mov_avg[:],alpha=0.05,nlags=_M,fft=True)
#tacf=data[:,1]
#plt.plot(tacf,acf,'--')
#print(acf)
#plt.legend(loc=(1.0,0.2),frameon=False)
#plt.xscale('log')
#plt.yscale('log')
#plt.savefig('stress_corr.jpg',dpi=300,bbox_inches='tight')
#plt.show()


