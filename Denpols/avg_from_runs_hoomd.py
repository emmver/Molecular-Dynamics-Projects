'''' Script is made to be used with a specific data tree and paths in mind. 
     For example, there is a plotting parameters file called in line 20
     Specific path structures given in lines 46, 53,55. 
     Careful alteration of such parameters is needed according to the user's needs.
'''

''' This script takes analyzed data from different runs of a simulation and averages them.
    Properly averaged data can then be used to do further analysis
'''

''' There are parameters such as:
 points: Number of points to be plotted. This can be turned into a time of simulation 
 if one knows the framestep (how many steps are needed to record a snapshot) and of course the timestep of the simulations
 gen: Generation of the dendron,  in this case needed for path navigation 
 N: Number of atoms on the backbone, once again used for path navigation 
 Ultimately, if one adjusts the script for different usage, there might not be a need for such parameters.
 '''


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
#############################################

os.chdir(sys.argv[1])

if os.path.isdir('./AVG')==False:
    os.mkdir('./AVG')

gen=int(sys.argv[2])
N=int(sys.argv[3])
file_list=glob.glob('./G%d_%d/RUN_*/results*/gyration*.txt'%(gen,N))
print('Starting Now: Gen=%d N=%d'%(gen,N))
print('------------------------------')

points=int(sys.argv[4])

print(file_list)
print('=============')
den_list=glob.glob('./G%d_%d/RUN_*/results*/*%d%d'%(gen,N,gen,N))
print(den_list)
print('=============')
bb_list=glob.glob('./G%d_%d/RUN_*/results*/*%d%d_bb'%(gen,N,gen,N))
print(bb_list)
print('==============')
gyr_tens=np.zeros((points,5))
data_den=np.zeros((186,2))
data_bb=np.zeros((186,2))
for i in range(len(file_list)):
    file=file_list[i]
    print('Working on file:',file)
    data=np.genfromtxt(file)
    if N!=1000:
        data=np.genfromtxt(file)
        data=data
   #num_frames=data[:,0].size
    gyr_tens+=data[:points]
    #w_path=os.path.dirname(os.path.abspath(file))
    #print(w_path)
    #os.chdir(w_path)
for i in range (len(den_list)):
    file=den_list[i]
    data=np.genfromtxt(file)
    data_den+=data
    file=bb_list[i]
    data=np.genfromtxt(file)
    data_bb+=data

gyr_tens=gyr_tens/(len(file_list))
data_den=data_den/(len(den_list))
data_bb=data_bb/(len(bb_list))


np.savetxt('./AVG/GYR_G%dN%d.avg'%(gen,N),gyr_tens)
np.savetxt('./AVG/Sq_full+G%dN%d.avg'%(gen,N),data_den)
np.savetxt('./AVG/Sq_bb_G%dN%d.avg'%(gen,N),data_bb)

#plot_shape(points,gyr_tens)

#os.chdir('../')
#print(os.getcwd())



#plt.loglog(data_den[:,0],data_den[:,1],color='#32a844',linestyle='-',label='Denpol')
#plt.loglog(data_bb[:,0],data_bb[:,1],color='#a87332',linestyle='--',label='Backbone')
#plt.legend(frameon=False)
#plt.savefig(sys.argv[2][:-1]+'_sq.jpg',dpi=300,bbox_inches='tight')



