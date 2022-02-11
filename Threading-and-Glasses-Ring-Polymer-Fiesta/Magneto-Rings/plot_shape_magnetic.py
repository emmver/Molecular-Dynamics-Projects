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
import glob 
import os 
import sys
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import KFold
from sklearn.neighbors import KernelDensity
import joblib
import scipy.stats as st


import seaborn as sns 
os.chdir(sys.argv[1])
#################### My Functions #######################################
bws=list()
def kde_scipy(x, x_grid,bw_idx,**kwargs):
    """Kernel Density Estimation with Scipy"""
    # Note that scipy weights its bandwidth by the covariance of the
    # input data.  To make the results comparable to the other methods,
    # we divide the bandwidth by the sample standard deviation here.
    filename='./distributions/bw.ls'
    if os.path.isfile(filename)==True:
        bwls=np.genfromtxt(filename)
        bandwidth=bwls[bw_idx]
    else:
        bandwidths = np.linspace(0.01, 0.5, 50)
        print("BW search...")
        grid = GridSearchCV(KernelDensity(kernel='gaussian'),
                            {'bandwidth': bandwidths},
                            cv=3,verbose=2,n_jobs=3)
        grid.fit(x[:, None]);
        bandwidth=grid.best_params_.get('bandwidth')
        print('Bandwidth for KDE:',bandwidth)
    kde = gaussian_kde(x, bw_method=bandwidth/ x.std(ddof=1) , **kwargs)
    bws.append(bandwidth)
    return kde.evaluate(x_grid)


def plot_corr_rg(files,xmin,xmax,store_name):
    idx=1
    if key=='tsa':
        plt.ylabel(r'$C^{R_g} [\tau]$')
        x=np.linspace(xmin,xmax,10000)
        start_from=0
        to_store='correlations'
    elif key=='rg':
        plt.ylabel(r'$R_g (t) [\sigma]$')
        x=np.linspace(xmin,xmax,10000)
        start_from=2
        to_store='timeseries'

    colors = plt.cm.Dark2(np.linspace(0,1,len(files)))
    data_interp=np.zeros(x.shape)
    for file in files: 
        data=np.genfromtxt(file)
        data[:,0]-=data[0,0]
        data_interp+=np.interp(x,data[:,0],data[:,1])
        plt.plot(data[:,0],data[:,1],color=colors[idx-1],marker='o',markersize=6,linestyle='none',label='RUN %d'%idx)
        #plt.plot(x,data_interp,color=colors[idx-1],lw=2)
        plt.xlabel(r'$\tau [\tau_0]$')
        plt.xscale('log')
        plt.tick_params(axis='x', which='minor')
        idx+=1
    plt.plot(x[start_from:],data_interp[start_from:]/(idx-1),color='r',lw=3)
    plt.xlim(500,5e6)
    plt.legend(frameon=False,loc=(1.0,0.0))
    plt.savefig(store_name+'.jpg',dpi=300,bbox_inches='tight')
    np.savetxt(store_name+'.avg',np.stack((x[:],data_interp[:]/(idx-1)),axis=-1))
    #plt.show()
    plt.clf()

def plot_eigen(files,xmin,xmax,rgsq,store_name):
    #files=glob.glob('*_eigen*.dat')
    start_from=2
    colors = plt.cm.Dark2(np.linspace(0,1,len(files)))
    x=np.linspace(xmin,xmax,10000)
    plot_avg=np.zeros((x.shape[0],3))
    for i in range(1,4):
        idx=1
        data_interp=np.zeros(x.shape)
        for file in files: 
            #print(file)
            data=np.genfromtxt(file)
            data[:,0]-=data[0,0]
            data_interp+=np.interp(x,data[:,0],data[:,i])
            plt.plot(data[:,0],data[:,i],color=colors[idx-1],marker='o',markersize=6,linestyle='none',label='RUN %d'%idx)
            #plt.plot(x,data_interp,color=colors[idx-1],lw=2)
            plt.ylabel(r'$\lambda_%d [\sigma^2]$'%i)
            plt.xlabel(r'$\tau [\tau_0]$')
            plt.xscale('log')
            plt.tick_params(axis='x', which='minor')
            idx+=1
        plot_avg[:,i-1]=data_interp/(idx-1)
        plt.plot(x[start_from:],data_interp[start_from:]/(idx-1),color='r',lw=3)
        plt.legend(frameon=False,loc=(1.0,0.0))
        plt.xlim(500,5e6)
        plt.savefig(store_name+'_eigen_%d.jpg'%i,dpi=300,bbox_inches='tight')
        np.savetxt(store_name+'_eigen_%d.avg'%i,np.stack((x[:],data_interp[:]/(idx-1)),axis=-1))
        #plt.show()
        plt.clf()
    colors = plt.cm.Set1(np.linspace(0,1,plot_avg.shape[1]))
    for i in range (plot_avg.shape[1]):
        plt.plot(x[start_from:],plot_avg[start_from:,i]/rgsq,color=colors[i],label=r'$\lambda_%d$'%(i+1))
    plt.xlabel(r'$\tau[\tau_0]$')
    plt.ylabel(r'$\dfrac{\lambda_i (t)}{\left<R_g^2\right>} $')
    plt.xscale('log')
    plt.legend(frameon=False)
    plt.savefig(store_name+'.jpg',dpi=300,bbox_inches='tight')
    plt.clf()




def rg_dist(files_rg,xmin,xmax,store_name):
    x=np.linspace(xmin,xmax,10000)
    start_from=1
    #files_rg=glob.glob('*rg*.dat')
    data_interp_rg=np.zeros(x.shape)
    idx=0
    for file in files_rg:
        data_rg=np.genfromtxt(file)
        data_rg[:,0]-=data_rg[0,0]
        data_interp_rg+=np.interp(x,data_rg[:,0],data_rg[:,1])
        rgs=data_rg[:,1]
        rgs=np.asarray(rgs)
        rg_iqr=iqr(rgs)#;print(rg_iqr)
        bw=2*rg_iqr/rgs.size**(1/3)#;print(bw)
        nbins=(rgs.max()-rgs.min())/bw#
        hist,bins=np.histogram(rgs,bins=int(nbins),density=True)
        idx+=1
        plt.plot(bins[:-1],hist,drawstyle='steps')
        plt.fill_between(bins[:-1],hist,step='pre',alpha=0.4)
    data_interp_rg=data_interp_rg/idx
    rgs=data_interp_rg
    rgs=np.asarray(rgs)
    rg_iqr=iqr(rgs)#;print(rg_iqr)
    bw=2*rg_iqr/rgs.size**(1/3)#;print(bw)
    nbins=(rgs.max()-rgs.min())/bw#
    hist,bins=np.histogram(data_interp_rg,bins=int(nbins),density=True)
    print("Bins:",int(nbins))
    plt.plot(bins[:-1],hist,drawstyle='steps')
    plt.fill_between(bins[:-1],hist,step='pre',alpha=0.4)
    plt.xlabel(r'$R_g [\sigma]$')
    plt.ylabel(r'$P(R_g) $')
    #plt.savefig('./distributions/histograms_rg.jpg',dpi=300,bbox_inches='tight')
    plt.clf()
    x_grid = np.linspace(rgs.min()-0.1*rgs.min(), rgs.max()+0.1*rgs.max(), 100)
    pdf=kde_scipy(data_interp_rg,x_grid,bw_idx)
    print('Area:',np.trapz(pdf,x=x_grid))
    mean=np.trapz(x_grid*pdf,x=x_grid); print('Mean Rg:',mean)
    std=np.sqrt(np.trapz((x_grid-mean)**2*pdf,x=x_grid)); print('Standard deviation:',std)
    plt.plot(x_grid,pdf)
    plt.fill_between(x_grid,pdf,alpha=0.4)
    gauss=(1/(std*np.sqrt(2*np.pi)))*np.exp(-0.5*(((x_grid-mean)/std)**2));print(gauss.max())
    plt.plot(x_grid,gauss,'k--')
    plt.ylabel(r'$P(R_g) $')
    plt.xlabel(r'$R_g [\sigma]$')
    plt.savefig(store_name+'.jpg',dpi=300,bbox_inches='tight')
    np.savetxt(store_name+'.avg',np.stack((x_grid,pdf),axis=-1))
    plt.clf()
    return mean**2

def eigen_dist(files_rg,rgsq,store_name):
    #files_rg=glob.glob('./timeseries/lambda_*.avg')
    idx=1
    colors = plt.cm.Set1(np.linspace(0,1,len(files_rg)))
    for file in files_rg:
        print(file)
        data_interp_rg=np.genfromtxt(file)[:,1]/rgsq
        #print(data_interp_rg)
        rgs=data_interp_rg
        rgs=np.asarray(rgs)
        rg_iqr=iqr(rgs)#;print(rg_iqr)
        bw=2*rg_iqr/rgs.size**(1/3)#;print(bw)
        nbins=(rgs.max()-rgs.min())/bw#
        hist,bins=np.histogram(data_interp_rg,bins=int(nbins),density=True)
        print("Bins:",int(nbins))
        x_grid = np.linspace(rgs.min()-0.1*rgs.min(), rgs.max()+0.1*rgs.max(), 100)
        pdf=kde_scipy(data_interp_rg,x_grid,idx)
        print('Area:',np.trapz(pdf,x=x_grid))
        mean=np.trapz(x_grid*pdf,x=x_grid); print('Mean Rg:',mean)
        std=np.sqrt(np.trapz((x_grid-mean)**2*pdf,x=x_grid)); print('Standard deviation:',std)
        plt.plot(x_grid,pdf,label=r'$lambda_%d$'%idx,color=colors[idx-1])
        plt.fill_between(x_grid,pdf,alpha=0.4,color=colors[idx-1])
        np.savetxt(store_name+'eigen_%d.avg'%idx,np.stack((x_grid,pdf),axis=-1))
        idx+=1
    plt.legend(frameon=False)
    plt.ylabel(r'$P(\dfrac{\lambda_i (t)}{\left<R_g^2\right>}) $')
    plt.xlabel(r'$\dfrac{\lambda_i (t)}{\left<R_g^2\right>} $')
    plt.savefig(store_name+'.jpg',dpi=300,bbox_inches='tight')
    plt.clf()
    idx+=1

# def prob_2D(x,y,xlabel,ylabel,save_fig):
#     # Define the borders
#     scale_f=5
#     deltaX = (max(x) - min(x))/10
#     deltaY = (max(y) - min(y))/10
#     xmin = min(x) - deltaX
#     xmax = max(x) + deltaX
#     ymin = min(y) - deltaY
#     ymax = max(y) + deltaY
#     print(xmin, xmax, ymin, ymax)
#     # Create meshgrid
#     xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
#     print(xx.shape,yy.shape)
#     positions = np.vstack([xx.ravel(), yy.ravel()])
#     values = np.vstack([x, y])
#     kernel = gaussian_kde(values)
#     f = np.reshape(kernel(positions).T, xx.shape)
#     fig = plt.figure(figsize=(8,8))
#     ax = fig.gca()
#     ax.set_xlim(xmin, xmax)
#     ax.set_ylim(ymin*scale_f, ymax*scale_f)
#     cfset = ax.contourf(xx, yy, f, cmap='coolwarm')
#     ax.imshow(np.rot90(f), cmap='coolwarm',extent=[0, 8, 0, 8])
#     cset = ax.contour(xx, yy, f, colors='k')
#     print('Cset:',cset)
#     ax.clabel(cset, inline=1, fontsize=10)
#     ax.set_xlabel(xlabel)
#     ax.set_ylabel(ylabel)
#     plt.savefig(save_fig)#dpi=300,bbox_inches='tight')
#     plt.show()
#     plt.clf()
#     #plt.title('2D Gaussian Kernel density estimation')


def shape_avg(files,rgfile,store_name,store_name_2):
    lambda_1=np.genfromtxt(files[0])
    time=lambda_1[:,0]
    lambda_1=lambda_1[:,1]
    lambda_2=np.genfromtxt(files[1])[:,1]
    lambda_3=np.genfromtxt(files[2])[:,1]
    ### Invariants ######
    i_1=lambda_1+lambda_2+lambda_3
    i_2=lambda_1*lambda_2+lambda_2*lambda_3+lambda_3*lambda_1
    i_3=lambda_1*lambda_2*lambda_3
    ##### Shape Params #####
    anis=1-3 * i_2/i_1**2 ## Shape anisotropy 
    plt.plot(time,anis)
    plt.xscale('log')
    plt.ylabel(r'$\delta$')
    plt.xlabel(r'$\tau [\tau_0]$')
    plt.savefig(store_name'_anisotropy.jpg',dpi=300,bbox_inches='tight')
    plt.clf()
    np.savetxt(store_name+'anisotropy.avg',np.stack((time,anis),axis=-1))
    prolat=((3*lambda_1-i_1)*(3*lambda_2-i_1)*(3*lambda_3-i_1))/i_1**3 ## Shape prolateness 
    plt.plot(time,prolat)
    plt.xscale('log')
    plt.ylabel(r'S*')
    plt.xlabel(r'$\tau [\tau_0]$')
    plt.savefig(store_name+'prolateness.jpg',dpi=300,bbox_inches='tight')
    plt.clf()
    np.savetxt(store_name+'prolateness.avg',np.stack((time,prolat),axis=-1))
    rgs=np.genfromtxt(rgfile[0])[:,1]
    # print(rgs.shape);print(prolat.shape)
    # sns.kdeplot(x=rgs, y=prolat, cmap="Reds", shade=False, bw_adjust=np.array(bws).mean(),thresh=0)
    # #plt.xlabel()
    # #plt.ylabel(r)
    plt.hist2d(rgs,prolat,bins=50,density=True,cmap='jet')
    plt.colorbar()
    plt.savefig(store_name_2+'_prolat_rg.jpg',dpi=300,bbox_inches='tight')
    plt.show()

    #prob_2D(rgs,prolat,r'$R_g [\sigma]$',r'S*','./2d_distributions/prolat_rg.jpg')




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

dir_names=['timeseries','distributions','2d_distributions','correlations']
for nam in dir_names:
    if os.path.isdir(nam)==False:
        os.mkdir(nam)

 

xmin=0
xmax=9.9e6
lambdas=np.arange(1,5)
for i in lambdas:
    files=glob.glob('active_rg_RUN_*_lambda_%d.dat'%i)
    store_name='./timeseries/rg_lambda_%d'%i
    plot_corr_rg(files,'tsa',xmin,xmax,store_name)
    store_name='./timeseries/corr_lambda_%d'%i
    plot_corr_rg(files,'rg',xmin,xmax,store_name)
    print('Done with gyration timeseries and correlation')
    files=glob.glob('active_rg_RUN_*_lambda_%d.dat'%i)
    store_name='./distributions/rg_lambda_%d'%i
    rgsq=rg_dist(files,xmin,xmax,store_name)
    print('Done with gyration distribution')
    files=glob.glob('active_eigen_RUN_*_lambda_%d.dat'%i)
    store_name='./timeseries/lambda_%d'%i
    plot_eigen(files,xmin,xmax,rgsq,store_name)
    print('Done with eigen values timeseries')
    files=glob.glob('./timeseries/lambda_%d*.avg'%i)
    store_name='./distributions/'
    eigen_dist(rgsq)
    print('Done with eigen values distribution')
    files=glob.glob('./timeseries/lambda_%d*.avg'%i)
    rgfile=glob.glob('./timeseries/rg_lambda_%d.avg'%i)
    store_name='./distributions/lambda_%d'%i
    store_name_2='./2d_distributions/lambda_%d'%i
    print(files)
    print(rgfile)
    shape_avg(files,rgfile,store_name,store_name_2)
    # files_eigen=glob.glob('*eigen*.dat')
    # rgsq=mean**2
    # print(data.shape)
    # for i in range(1,4):
    #     for file in files_eigen:
    #         data_eigen=np.genfromtxt(file)
    #         data_eigen[:,0]-=data_rg[0,0]
    #         data_interp_eigen+=np.interp(x,data_eigen[:,0],data_eigen[:,i])
    bws=np.array(bws)
    np.savetxt('./distributions/bw.ls',bws)