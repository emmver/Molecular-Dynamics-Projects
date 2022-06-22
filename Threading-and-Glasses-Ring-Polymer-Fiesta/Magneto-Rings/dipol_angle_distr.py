'''
This script finds and reads in the .dat files (see variables filename and files). These files are produced during analysis.
These files contain a table with dimensions ("separation distance", "number of frames"). Basically each row contains the average value corresponding to 
some separation distance for each frame (columns).

Usage: python3 dipol_angle_distr.py /path/to/files DP

DP is the degree of polymerization. There are other variables in the script which are very important:

1) how_many: percentage of magnetic monomers 
2) nchains: number of chains 
3) lambdas: list with lambdas to look for and calculate 
4) histo_steps: histograms to plot for specific separation distances 

The Plot style section can be erased without much changing in the script. 
Otherwise the provided plostyle.mpl file should be used. 
'''

import numpy as np 
from matplotlib import pyplot as plt
import sys
import matplotlib as mpl
import glob 
import os
import seaborn as sns
import warnings
from scipy.stats import iqr

warnings.filterwarnings("ignore")
####################################################################

def fred_diac(data):
    '''
    This function calculates the number of bins to be used in a histogram.
    It follows the Friedman-Diaconis rule --> https://en.wikipedia.org/wiki/Freedman%E2%80%93Diaconis_rule#:~:text=For%20a%20set%20of%20empirical,of%20the%20theoretical%20probability%20distribution.
    iqr: interquartile range, which is a measure of statistical dispersion
    bw: bin width 
    '''
    data_iqr=iqr(data)#;print(rg_iqr)
    bw=2*data_iqr/data.size**(1/3)#;print(bw)
    nbins=(data.max()-data.min())/bw#
    return nbins

def error_bands(data,labeled,hue):
    '''
    Function to produce plots with the error bars shown as a filled area around the curve.
    x: independent variable
    y: dependent variable
    yerr: error of the dependent variable
    '''
    x=data[:,0]
    y_mean=data[:,1]
    y_std=data[:,2]
    error = 0.5*y_std
    lower = y_mean - error
    upper = y_mean + error
    ax.plot(x, y_mean,color=hue,label=labeled)
    ax.plot(x, lower,lw=0.5,linestyle=':', color=hue, alpha=0.5)
    ax.plot(x, upper,lw=0.5,linestyle=':', color=hue, alpha=0.5)
    ax.fill_between(x, lower, upper,color=hue, alpha=0.2)
    ax.legend(frameon=False)
    ax.set_xlabel('s')
    ax.set_ylabel(r'$\left<cos(\theta_{dip})\right>$')
    fig.savefig(directory+'/angle_correlation.png',dpi=300,bbox_inches='tight')

#############  Plotting  Style #####################
mpl.rcParams['figure.dpi'] = 400
plt.rcParams["font.family"] = "Ubuntu"
plt.style.use('~/plotstyle.mplstyle') ## Here input the path to the plotstyle.mplstyle file
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
###################################################

os.chdir(sys.argv[1]) #work in the sys.argv[1] directory 
directory='dipole_angle_distribution'
if not os.path.exists(directory):
    os.makedirs(directory) #create directory to store files, if it does not exist
dict_file=['active'] 
nchains=1
lambdas=np.array([1,4,9,16])
how_many=0.5
colors_l=mpl.cm.jet(np.linspace(0,1,lambdas.size)) #define colormap to use
hist_list=[]
DP=int(sys.argv[2])


'''
Following loops go over the activity list (dict_file) to analyze different parts of the polymer (only active, only passive or total).
Then the second loop goes through lambdas available. And third loops goes through the files found with files=glob.glob(filename)
'''
for count,d in enumerate(dict_file):
    ax=plt.axes()
    for count_l,i in enumerate(lambdas):
        print('lambda: '%i)
        tostore=np.zeros((int(DP/2)-1,1))
        filename=d+'dipol_dist_RUN_*_lambda_%d.dat'%i
        files=glob.glob(filename)
        print(files)
        for file in files: 
            data=np.genfromtxt(file)
            tostore=np.append(tostore,data[:,1:],axis=1)
            #print("lambda %d"%(i))
            #print(data.shape)
        tostore=tostore[:,1:]
        np.savetxt(directory+'/fordistr_dipol_angles_lambda_%d.dat'%i,tostore)  


dist_files=glob.glob(directory+'/fordistr_dipol_angles_lambda_*.dat')
print(dist_files)
histo_steps=[80,60,50,25,10,5,3,2,0]
colors = plt.cm.copper(np.linspace(1,0,len(histo_steps)))
'''
Following loops go over the files created above and plot histograms for the distribution of cosinuses between dipoles along the chain. 
'''

for count_l,i in enumerate(lambdas):
    means=np.zeros((int(DP/2)-1,3))
    filename=directory+'/fordistr_dipol_angles_lambda_%d.dat'%i
    data=np.genfromtxt(filename)
    plt.figure()
    print(i)
    #step=10
    for histo_idx,ii in enumerate(histo_steps):
        nbins=int(fred_diac(data[ii,:]))
        sns.histplot(data[ii,:],label='$s=%d$'%ii,bins=nbins,color=colors[histo_idx],element='poly')
    for ii in range(data.shape[0]):
        means[ii,0]=ii
        means[ii,1]=np.mean(data[ii,:]) # mean value from all configurations
        means[ii,2]=np.mean(data[ii,:]) # mean error from all configurations

    plt.xlabel(r'$cos(\theta_{dip})$')
    plt.ylabel('PDF')
    plt.legend(loc=(1.0,0.43),frameon=False,ncol=1,fontsize=11)
    plt.savefig(directory+"/dist_progression_lambda_%d.png"%i,dpi=300,bbox_inches='tight') #Store histograms
    np.savetxt(directory+'/avg_dipol_angle_lambda_%d.dat'%i,means) #Store mean and error

'''
This part plots the average cosinus as a function of separation distance s along with error bars.
Loads files created just above.
'''

colors=colors = plt.cm.copper(np.linspace(1,0,len(lambdas)))

fig, ax = plt.subplots(figsize=(9,5))
for count_l,i in enumerate(lambdas):
    filename=directory+'/avg_dipol_angle_lambda_%d.dat'%i
    data=np.genfromtxt(filename)
    label=r'$\lambda=%d$'%i
    error_bands(data,label,colors[count_l])


