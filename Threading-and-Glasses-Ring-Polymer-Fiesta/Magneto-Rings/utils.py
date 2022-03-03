import numpy as np 
#import MDAnalysis as mda
from matplotlib import pyplot as plt 
from numba import njit
import sys
from timeit import default_timer as timer
from scipy.stats import gaussian_kde,iqr
import ctypes
import matplotlib as mpl
from numba import njit
import glob 
import os 
import sys
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import KFold
from sklearn.neighbors import KernelDensity
import scipy.stats as st

############### Analysis Functions ############################
def read_traj(f,npart):
    out=[]
    for i in range (npart):
            line = f.readline()
            elems=str.split(line.strip()," ")
            #print(elems)
            for el in elems: 
                #print(el)
                out.extend([float(el)])
                if len(elems)>0: elsize=len(elems);
    out=np.reshape(np.asarray(out),(int(len(out)/elsize),elsize))
    f.readline()
    f.readline()
    return out

def find_xmax(files):
    xmax=1e50
    for i in files: 
        data=np.genfromtxt(i)
        if data[-1,0]<xmax:
            xmax=data[-1,0]
    return xmax

######### Calculate Gyration Tensor #######################3
@njit(fastmath=True)
def rg_tens (r):
    
    N = r.shape[0]
    r2=np.zeros((3,3))
    for i in range (N):
        for j in range (N):
            r2[0,0]+=(r[i,0]-r[j,0])*(r[i,0]-r[j,0])#diag
            r2[1,0]+=(r[i,1]-r[j,0])*(r[i,1]-r[j,0])#ndiag
            r2[2,0]+=(r[i,2]-r[j,0])*(r[i,2]-r[j,0])#ndiag
            r2[0,1]+=(r[i,0]-r[j,1])*(r[i,0]-r[j,1])#ndiag
            r2[0,2]+=(r[i,0]-r[j,2])*(r[i,0]-r[j,2])#ndiag
            r2[1,1]+=(r[i,1]-r[j,1])*(r[i,1]-r[j,1])#diag
            r2[2,1]+=(r[i,2]-r[j,1])*(r[i,2]-r[j,1])#ndiag
            r2[1,2]+=(r[i,1]-r[j,2])*(r[i,1]-r[j,2])#ndiag
            r2[2,2]+=(r[i,2]-r[j,2])*(r[i,2]-r[j,2])# diag
    r2=r2/(2*N**2)
    #rg2=r2[0,0]**2+r2[1,1]**2+r2[2,2]**2
    return r2

#################### MSDs Lib ############
@njit(fastmath=True)
def coms_times(pos,pos_coms,Dpolsm,nsmall,frames):
    for i in range (frames):
        data_small=pos[i]
        for j in range (nsmall):
            workdata=data_small[j*Dpolsm:(j+1)*Dpolsm]
            #print(workdata.shape)
            comx=np.sum(workdata[:,0])
            comy=np.sum(workdata[:,1])
            comz=np.sum(workdata[:,2])
            com = np.array([comx,comy,comz])/(workdata[:,0].size)
            pos_coms[i][j]=com
    return pos_coms
#### MSD Libs #######

def _autocorr_fft(x):
    """Compute the autocorrelation of 1d signal using FFT."""
    N = x.shape[0]
    # 2*N because of zero-padding to compute non-cyclic correlation
    f = np.fft.fft(x, n=2*N)
    power_spectrum = f * f.conjugate()
    result = np.fft.ifft(power_spectrum)
    # the autocorrelation in the usual convention B
    result = (result[:N]).real
    return result

@njit 
def _msd_fft_compute_s1(r):
    """Compute the non-FFT part of the MSD in the FFT method."""
    N = r.shape[0]
    # Compute the squared distances
    D = np.square(r).sum(axis=1)
    # Apply the recursive algorithm
    s1 = np.zeros(N)
    s1[0] = 2*D.sum()
    for n in range(1, N):
        s1[n] = s1[n-1] - D[n-1] - D[N-n]
    return s1

def self_fft(r):
    """Compute the self MSD from a single particle trajectory using FFT.

    Based on the algorithm outlined in Section 4.2 of the following paper:
    https://doi.org/10.1051/sfn/201112010
    """
    N = r.shape[0]
    # Compute the non-FFT part
    s1 = _msd_fft_compute_s1(r)
    # Compute the FFT-part separately over each position component
    s2 = np.zeros_like(s1)
    for i in range(r.shape[1]):
        s2 += _autocorr_fft(r[:, i])
    return (s1 - 2*s2) / np.arange(N, 0, -1)

#### MSIDs ########
@njit 
def compute_msid(positions, start, stop, num_rings, num_mons):
    
    s = np.arange(1, stop-start)
    msid = np.zeros_like(s, dtype=np.float64)
    upds = np.zeros_like(s, dtype=np.int64)
    
    for n in range(num_rings):
        r = positions[n*num_mons:(n+1)*num_mons]
        for index, ds in enumerate(s):            
            i = start
            j = start + ds
            while j < stop:
                dr = r[j] - r[i]
                msid[index] += dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]
                upds[index] += 1    
                j += 1
                i += 1
    return s, msid/upds
###################################################################

################# Plotting Functions ##############################
def plot_corr_rg(key,xmin,xmax,files,points):
    idx=1
    if key=='tsa':
        plt.ylabel(r'$C^{R_g} [\tau]$')
        x=np.linspace(xmin,xmax,points)
        start_from=0
        to_store='correlations'
    elif key=='rg':
        plt.ylabel(r'$R_g (t) [\sigma]$')
        x=np.linspace(xmin,xmax,points)
        start_from=2
        to_store='timeseries'
    colors = plt.cm.Dark2(np.linspace(0,1,len(files)))
    data_interp=np.zeros(x.shape)
    add_path=os.path.dirname(os.path.abspath(files[0]))+'/'+to_store
    print('Add path:',add_path)
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
    plt.savefig(add_path+'/'+key+'.jpg',dpi=300,bbox_inches='tight')
    np.savetxt(add_path+'/average_%s.avg'%key,np.stack((x[:],data_interp[:]/(idx-1)),axis=-1))
    #plt.show()
    plt.clf()

def plot_eigen(xmin,xmax,rgsq,files,points):
    start_from=2
    colors = plt.cm.Dark2(np.linspace(0,1,len(files)))
    x=np.linspace(xmin,xmax,points)
    plot_avg=np.zeros((x.shape[0],3))
    add_path=os.path.dirname(os.path.abspath(files[0]))
    print('Add path:',add_path)
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
        plt.savefig(add_path+'/timeseries/lambda_%d.jpg'%i,dpi=300,bbox_inches='tight')
        np.savetxt(add_path+'/timeseries/lambda_%d_average.avg'%i,np.stack((x[:],data_interp[:]/(idx-1)),axis=-1))
        #plt.show()
        plt.clf()
    colors = plt.cm.Set1(np.linspace(0,1,plot_avg.shape[1]))
    for i in range (plot_avg.shape[1]):
        plt.plot(x[start_from:],plot_avg[start_from:,i]/rgsq,color=colors[i],label=r'$\lambda_%d$'%(i+1))
    plt.xlabel(r'$\tau[\tau_0]$')
    plt.ylabel(r'$\dfrac{\lambda_i (t)}{\left<R_g^2\right>} $')
    plt.xscale('log')
    plt.legend(frameon=False)
    plt.savefig(add_path+'/timeseries/eigen_values_avg.jpg',dpi=300,bbox_inches='tight')
    plt.clf()

def rg_dist(xmin,xmax,files_rg,points):
    x=np.linspace(xmin,xmax,points)
    start_from=1
    #files_rg=glob.glob('*rg*.dat')
    data_interp_rg=np.zeros(x.shape)
    idx=0
    add_path=os.path.dirname(os.path.abspath(files[0]))
    print('Add path:',add_path)
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
    plt.savefig(add_path+'/distributions/histograms_rg.jpg',dpi=300,bbox_inches='tight')
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
    plt.savefig(add_path+'/distributions/kde_rg.jpg',dpi=300,bbox_inches='tight')
    np.savetxt(add_path+'/distributions/avg_kde_rg.avg',np.stack((x_grid,pdf),axis=-1))
    plt.clf()
    return mean**2

def eigen_dist(rgsq,files):
    #files_rg=glob.glob('./timeseries/lambda_*.avg')
    idx=1
    colors = plt.cm.Set1(np.linspace(0,1,len(files_rg)))
    add_path=os.path.dirname(os.path.dirname(os.path.abspath(files[0])))
    print('Add path:',add_path)
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
        np.savetxt('.os/distributions/avg_kde_lambda_%d.avg'%idx,np.stack((x_grid,pdf),axis=-1))
        idx+=1
    plt.legend(frameon=False)
    plt.ylabel(r'$P(\dfrac{\lambda_i (t)}{\left<R_g^2\right>}) $')
    plt.xlabel(r'$\dfrac{\lambda_i (t)}{\left<R_g^2\right>} $')        
    plt.savefig('./distributions/eigen_distr.jpg',dpi=300,bbox_inches='tight')
    plt.clf()
    idx+=1

def shape_avg(files,rgfile):
    add_path=os.path.dirname(os.path.dirname(os.path.abspath(rgfiles[0])))
    print('Add path:',add_path)
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
    plt.savefig(add_path+'/timeseries/anisotropy.jpg',dpi=300,bbox_inches='tight')
    plt.clf()
    np.savetxt(add_path+'/timeseries/anisotropy.avg',np.stack((time,anis),axis=-1))
    prolat=((3*lambda_1-i_1)*(3*lambda_2-i_1)*(3*lambda_3-i_1))/i_1**3 ## Shape prolateness 
    plt.plot(time,prolat)
    plt.xscale('log')
    plt.ylabel(r'S*')
    plt.xlabel(r'$\tau [\tau_0]$')
    plt.savefig(add_path+'/timeseries/prolateness.jpg',dpi=300,bbox_inches='tight')
    plt.clf()
    np.savetxt(add_path+'/timeseries/prolateness.avg',np.stack((time,prolat),axis=-1))
    rgs=np.genfromtxt(rgfile[0])[:,1]
    # print(rgs.shape);print(prolat.shape)
    # sns.kdeplot(x=rgs, y=prolat, cmap="Reds", shade=False, bw_adjust=np.array(bws).mean(),thresh=0)
    # #plt.xlabel()
    # #plt.ylabel(r)
    plt.hist2d(rgs,prolat,bins=50,density=True,cmap='jet')
    plt.colorbar()
    plt.savefig(add_path+'/2d_distributions/prolat_rg.jpg',dpi=300,bbox_inches='tight')
    plt.show()
##################################################################