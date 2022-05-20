import numpy as np 
from matplotlib import pyplot as plt 
import sys
import os
import matplotlib.pylab as pl
import glob
from matplotlib.ticker import FuncFormatter
import matplotlib.pylab as pl

### My plotting style is inputted here #####
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 120
plt.rcParams["font.family"] = "Ubuntu"
plt.style.use('C:/Users/Toumba/Documents/plotstyle.mplstyle')
#plt.style.use('~/plotstyle.mplstyle')
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


N=np.linspace(10,300,1000)
A=np.linspace(0.1,0.5,10)
colors = pl.cm.jet(np.linspace(0,1,A.size))

for i in A:
    T=i/np.log(N)
    plt.plot(N,T,label='A=%.2f'%i)
plt.legend(loc=(1.0,0.0),frameon=False)
plt.show()
