import numpy as np 
#import MDAnalysis as mda
from matplotlib import pyplot as plt 
import matplotlib as mpl
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
from utils import *
import scipy.stats as st

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

biden_color='#2986cc'
trump_color='#cc0000'
#########################################################################################################

class analysis_equil():
    def __init__(self,path,frames,warm_frames,run_steps,nchains,DP,L,arch_key):
        self.path=path
        self.frames=frames
        self.warm_frames=warm_frames
        self.run_steps=run_steps
        self.nchains=nchains
        self.DP=DP
        self.L=L
        self.arch_key=arch_key
        ######################
        os.chdir(self.path)
        self.directory='msds_results_pass'
        if not os.path.exists(self.directory):
            os.makedirs(self.directory)
        self.gyr_directory='gyration_results_pass'
        if not os.path.exists(self.gyr_directory):
            os.makedirs(self.gyr_directory)
        self.msid_dir='msid_results_pass'
        if not os.path.exists(self.msid_dir):
            os.makedirs(self.msid_dir)
        # etoe_directory='endtoend_results_pass'
        # if not os.path.exists(etoe_directory):
        #     os.makedirs(etoe_directory)
        self.npart=int(nchains*DP)
        if arch_key=='ring':
            self.bonds=self.npart
        elif arch_key=='linear':
            self.bonds=self.npart-1
        self.angles=0
        self.files=glob.glob("RUN_*/polymer.vtf")
        print(self.files,len(self.files))

    def get_frames(self,file):
        print("Working on file",file)
        f=open(file,'r')
        self.pos=np.zeros((int(frames-skipframes),(npart),3))
        for ii in range (int(self.npart)+self.bonds+self.angles+2):
                f.readline()
        for ii in range (self.warm_frames):
            data=read_traj(f,self.npart)
        self.count_frames=0

        pbar = tqdm(total = frames-skipframes,desc='Creating Pos array')
        while 1: 
            try: 
                data=read_traj(f,npart)
                self.pos[count_frames]=data[:,1:]
                self.count_frames+=1
                pbar.update(1)
            except Exception as e:
                print("Except:",e)
                print('stopped at',self.count_frames)
                self.pos=self.pos[:count_frames,:,:]
                break
        pbar.close()
        self.time = np.zeros(self.pos.shape[0])
        count=0
        for tt in tqdm(range(self.pos.shape[0]),desc='Time array'):
            self.time[tt]=count*self.run_steps
            count+=1
        f.close()

    def calc_g1(self):
        msds = np.zeros((self.npart, self.pos.shape[0]))
        for ii in tqdm(range(self.pos.shape[1]),desc='Monomer MSDs'):
            msds[ii] = self_fft(self.pos[:,ii])
        np.savetxt(self.directory+'/'+'_g1_RUN_%d.dat'%(self.run_count),(np.stack((self.time[:],msds.mean(axis=0)[:]),axis=-1)))
        
    def calc_g3(self):
        pos_coms=np.zeros(self.pos.shape[0],self.nchains,3))
        pos_coms=coms_times(self.pos,pos_coms,self.DP,self.nchains,self.pos.shape[0])
        g3s = np.zeros((self.nchains,self.pos.shape[0]))
        for ii in tqdm(range(pos_coms.shape[1]),desc='CoM MSDs'):
            g3s[ii] = self_fft(pos_coms[:,ii])
        np.savetxt(self.directory+'/'+'_g3_RUN_%d.dat'%(self.run_count),(np.stack((self.time[:],g3s[0]),axis=-1)))

    def calc_gyration_tensor(self):
        rgs_tens= np.zeros((nchains, frames-skipframes,3,3))
        for ii in tqdm(range (self.pos.shape[0]),desc='Gyration Tensor Calc & Storage'):
            for j in range (self.nchains):
                    temp_data=self.pos[ii,j*self.DP:(j+1)*self.DP]
                    rgs_tens[j,ii]=rg_tens(temp_data)
        radii=[]
        radii_all=np.zeros((self.pos.shape[0],self.nchains))
        rg_tensor_eigen=np.zeros((self.pos.shape[0],3))
        for ii in tqdm(range (self.pos.shape[0]),desc='Calculating Rg and shape'):
            tocalc=np.diag(rgs_tens[:,ii][0][0:3])

            tocalc=tocalc[tocalc.argsort()]
            temp_rg=0
            temp_rs=[]
            temp_eigen=np.array([0,0,0])
            for j in range (0,self.nchains):
                temp=(tocalc[0])+(tocalc[1])+(tocalc[2])
                temp_rg+=temp
                temp_eigen=temp_eigen+tocalc
                radii_all[ii,j]=temp**0.5
            radii.append(temp_rg**0.5)
            rg_tensor_eigen[ii]=temp_eigen
        if j==0:
            j=1
        radii=np.asarray(radii)/j
        rg_tensor_eigen=rg_tensor_eigen/j
        corr=np.stack((self.time,radii),axis=-1)
        _M =corr[:].shape[0]
        acf, ci = sm.tsa.acf(radii[skipframes:], alpha=0.05,nlags=_M,fft=True)
        tacf=np.arange(_M)*1000
        np.savetxt(self.gyr_directory+'/'+'tsa_RUN_%d.dat'%(self.run_count),(np.stack((tacf,acf),axis=-1)))
        np.savetxt(self.gyr_directory+'/'+'_rg_RUN_%d.dat'%(self.run_count),(np.stack((self.time[:],radii),axis=-1)))
        np.savetxt(self.gyr_directory+'/'+'_eigen_RUN_%d.dat'%(self.run_count),(np.column_stack([self.time[:],rg_tensor_eigen])))

    def calc_msid(self):
            msid_dist=np.zeros((np.arange(0,int(self.pos.shape[1]/2)).size,self.count_frames+1))
            #print(msid_dist.shape)
            for ii in range (self.count_frames):
                msid_dist[:-1,0],msid_dist[:-1,ii]=compute_msid(self.pos[ii,:], 0, int(self.pos.shape[1]/2), self.nch, int(self.pos.shape[1]))
                #print('Frame for MSID:',ii)
            msid_dist[:-1,-1]=np.mean(msid_dist[:-1,1:-1],axis=1)
            np.savetxt(self.msid_dir+'/'+'msids_RUN_%d.dat'%(self.run_count),msid_dist[:-1,:])

    def run_equil(self):
        self.run_count=1
        for j in range(len(self.paths)):
            self.file=self.paths[j]
            self.get_frames(self.file)
            self.calc_g1()
            self.calc_g3()
            self.calc_gyration_tensor()
            self.calc_msid()

################################################################################################################
class analysis_dipol():
    def __init__(self,path,frames,warm_frames,run_steps,time_step,nchains,DP,L,how_many,dipol,arch_key):
        self.path=path
        self.frames=frames
        self.warm_frames=warm_frames
        self.run_steps=run_steps
        self.nchains=nchains
        self.DP=DP
        self.L=L
        self.how_many=how_many
        self.arch_key=arch_key
        self.part_dict=['total','active','passive']
        self.dipol=dipol
        self.time_step=time_step
        ######################
        os.chdir(self.path)
        self.directory='msds_results'
        if not os.path.exists(self.directory):
            os.makedirs(self.directory)
        self.gyr_directory='gyration_results'
        if not os.path.exists(self.gyr_directory):
            os.makedirs(self.gyr_directory)
        self.msid_dir='msid_results'
        if not os.path.exists(self.msid_dir):
            os.makedirs(self.msid_dir)
        # etoe_directory='endtoend_results_pass'
        # if not os.path.exists(etoe_directory):
        #     os.makedirs(etoe_directory)
        self.npart=int(nchains*DP)
        if arch_key=='ring':
            self.bonds=self.npart
        elif arch_key=='linear':
            self.bonds=self.npart-1
        self.angles=0
        self.files=glob.glob("RUN_*/lambda_%d/polymer_magnetic.vtf"%(self.dipol))
        print(self.files,len(self.files))

    def get_frames(file,act_key):
        self.act_key=act_key
        print("Working on file",file)
        f=open(file,'r')
        if self.act_key=='active': 
            self.npart=npart*how_many
            self.pos=np.zeros((int(self.frames),int(self.npart),3))
        elif self.act_key=='passive':
            self.npart=npart*(1-how_many)
            self.pos=np.zeros((int(self.frames),int(self.npart),3))
        else:
            self.pos=np.zeros((int(self.frames),(self.npart),3))
        for ii in range (int(self.npart)+self.bonds+self.angles+2):
                f.readline()
        for ii in range (self.warm_frames):
            data=read_traj(f,self.npart)
        self.count_frames=0
        pbar = tqdm(total = self.frames-self.warm_frames,desc='Creating Pos array')
        while 1: 
            try: 
                data=read_traj(f,npart)
                if self.act_key='total':
                    pos[count_frames]=data[:,1:]
                elif self.act_key=='active':
                    pos[count_frames]=data[:100,1:]
                elif self.act_key=='passive':
                    pos[count_frames]=data[100:,1:]
                self.pos[count_frames]=data[:,:]
                self.count_frames+=1
                pbar.update(1)
            except Exception as e:
                print("Except:",e)
                self.pos=self.pos[:count_frames,:,:]
                print('stopped at',self.count_frames)
                break
        self.time = np.zeros(self.pos.shape[0])
        count=0
        for tt in tqdm(range(self.pos.shape[0]),desc='Time array'):
            self.time[tt]=count*self.run_steps*self.time_step
            count+=1
        pbar.close()

    def calc_g1(self):
        msds = np.zeros((self.npart, self.pos.shape[0]))
        for ii in tqdm(range(self.pos.shape[1]),desc='Monomer MSDs'):
            msds[ii] = self_fft(self.pos[:,ii])
        np.savetxt(self.directory+'/'+self.act_key+'_g1_RUN_%d_lambda_%d.dat'%(self.run_count,self.dipol),(np.stack((self.time[:],msds.mean(axis=0)[:]),axis=-1)))
        
    def calc_g3(self):
        self.pos_coms=np.zeros(self.pos.shape[0],self.nchains,3))
        self.pos_coms=coms_times(self.pos,pos_coms,self.DP,self.nchains,self.pos.shape[0])
        g3s = np.zeros((self.nchains,self.pos.shape[0]))
        for ii in tqdm(range(pos_coms.shape[1]),desc='CoM MSDs'):
            g3s[ii] = self_fft(pos_coms[:,ii])
        np.savetxt(self.directory+'/'+self.act_key+'_g3_RUN_%d_lambda_%d.dat'%(self.run_count,self.dipol),(np.stack((self.time[:],g3s[0]),axis=-1)))
    
    def traj_com(self):


    def calc_gyration_tensor(self):
        rgs_tens= np.zeros((nchains, frames-skipframes,3,3))
        for ii in tqdm(range (self.pos.shape[0]),desc='Gyration Tensor Calc & Storage'):
            for j in range (self.nchains):
                    temp_data=self.pos[ii,j*self.DP:(j+1)*self.DP]
                    rgs_tens[j,ii]=rg_tens(temp_data)
        radii=[]
        radii_all=np.zeros((self.pos.shape[0],self.nchains))
        rg_tensor_eigen=np.zeros((self.pos.shape[0],3))
        for ii in tqdm(range (self.pos.shape[0]),desc='Calculating Rg and shape'):
            tocalc=np.diag(rgs_tens[:,ii][0][0:3])

            tocalc=tocalc[tocalc.argsort()]
            temp_rg=0
            temp_rs=[]
            temp_eigen=np.array([0,0,0])
            for j in range (0,self.nchains):
                temp=(tocalc[0])+(tocalc[1])+(tocalc[2])
                temp_rg+=temp
                temp_eigen=temp_eigen+tocalc
                radii_all[ii,j]=temp**0.5
            radii.append(temp_rg**0.5)
            rg_tensor_eigen[ii]=temp_eigen
        if j==0:
            j=1
        radii=np.asarray(radii)/j
        rg_tensor_eigen=rg_tensor_eigen/j
        corr=np.stack((self.time,radii),axis=-1)
        _M =corr[:].shape[0]
        acf, ci = sm.tsa.acf(radii[:], alpha=0.05,nlags=_M,fft=True)
        tacf=np.arange(_M)*1000
        np.savetxt(self.gyr_directory+'/'+self.act_key+'_tsa_RUN_%d_lambda_%d.dat'%(self.run_count,%self.dipol),(np.stack((tacf,acf),axis=-1)))
        np.savetxt(self.gyr_directory+'/'+self.act_key+'_rg_RUN_%d_lambda_%d.dat'%(self.run_count,%self.dipol),(np.stack((self.time[:],radii),axis=-1)))
        np.savetxt(self.gyr_directory+'/'+self.act_key+'_eigen_RUN_%d.dat'%(self.run_count),(np.column_stack([self.time[:],rg_tensor_eigen])))

    def calc_msid(self):
            msid_dist=np.zeros((np.arange(0,int(self.pos.shape[1]/2)).size,self.count_frames+1))
            #print(msid_dist.shape)
            for ii in range (self.count_frames):
                msid_dist[:-1,0],msid_dist[:-1,ii]=compute_msid(self.pos[ii,:], 0, int(self.pos.shape[1]/2), self.nch, int(self.pos.shape[1]))
                #print('Frame for MSID:',ii)
            msid_dist[:-1,-1]=np.mean(msid_dist[:-1,1:-1],axis=1)
            np.savetxt(self.msid_dir+'/'+self.act_key+'_msids_RUN_%d.dat'%(self.run_count),msid_dist[:-1,:])

    def analyze(self):
        self.run_count=1
        for j in range(len(self.paths)):
            self.file=self.paths[j]
            self.get_frames(self.file)
            self.calc_g1()
            self.calc_g3()
            self.calc_gyration_tensor()
            self.calc_msid()
            self.run_count+=1

    def analyze_msds(self):
        self.run_count=1
        for j in range(len(self.paths)):
            self.file=self.paths[j]
            self.get_frames(self.file)
            self.calc_g1()
            self.calc_g3()
            self.run_count+=1

    def analyze_shape(self):
        self.run_count=1
        for j in range(len(self.paths)):
            self.file=self.paths[j]
            self.calc_gyration_tensor()
            self.run_count+=1
    
    def analyze_msid(self):
        self.run_count=1
        for j in range(len(self.paths)):
            self.file=self.paths[j]
            self.get_frames(self.file)
            self.calc_msid()
            self.run_count+=1


################################################################################################################
class plot_equil():





################################################################################################################

class plot_dipol():
    def __init__(self,path,xmin,xmax,act_key,dipol,points,how_many,nchains,DP)
        self.path=path
        self.xmin=xmin
        self.xmax=xmax
        self.act_key=act_key
        self.dipol=dipol
        self.points=points
        self.how_many=how_many
        self.nchains=nchains
        self.DP=DP
        os.chdir(self.path)
        os.mkdir(self.act_key)
        os.popen('cp *'+self.act_key+'*.dat ./'+self.act_key) 
        self.dir_names=['timeseries','distributions','2d_distributions','correlations']
        self.x=np.linspace(self.xmin,self.xmax,self.points)
        for nam in self.dir_names:
            if os.path.isdir(self.act_key+'/'+nam)==False:
                print('Make Dir:',self.act_key+'/'+nam)
                os.mkdir(self.act_key+'/'+nam)
        self.gyr_path='./gyration_results'
        self.msds_path='./msds_results'
        self.msid_path='./msid_results'
        self.gyr_t_files=glob.glob(self.gyr_path+'/'+self.act_key+'/*%s*.dat'%key)
        self.gyr_dist_files=glob.glob(self.gyr_path+'/'+self.act_key+'/*rg*.dat')
        self.eigen_t_files=glob.glob(self.gyr_path+'/'+self.act_key+'/timeseries/lambda_*.avg')
        self.eigen_dist_files=glob.glob(self.gyr_path+'/'+self.act_key+'/*_eigen*.dat')
        self.shape_files=glob.glob(self.act_key+'/timeseries/lambda*.avg')
        self.shape_files_rg=glob.glob(self.act_key+'/timeseries/*rg.avg')
        self.g1_files=glob.glob(self.msds_path+'/'+self.act_key+'_g1_RUN_*_lambda_%d.dat'%self.dipol)
        self.g3_files=glob.glob(self.msds_path+'/'+self.act_key+'_g1_RUN_*_lambda_%d.dat'%self.dipol)
        self.g1_files_three_sp=glob.glob(self.msds_path+'/avg'+'*_g1_lambda_%d.avg'%self.dipol)
        self.g3_files_three_sp=glob.glob(self.msds_path+'/avg'+'*_g3_lambda_%d.avg'%self.dipol)
        self.g1_files_lambda=glob.glob(self.msds_path+'/avg'+self.act_key+'_g1_lambda_*.avg')
        self.g3_files_lambda=glob.glob(self.msds_path+'/avg'+self.act_key+'_g3_lambda_*.avg')

    def plot_gyr(self):
        print('Files to plot')
        print(self.gyr_t_files)
        plot_corr_rg('tsa',self.xmin,self.xmax,self.gyr_t_files,self.points)
        plot_corr_rg('rg',self.xmin,self.xmax,self.gyr_t_files,self.points)
        print('Done with gyration timeseries and correlation')
        self.rgsq=rg_dist(self.xmin,self.xmax,self.points)
        print('Done with gyration distribution')
        plot_eigen(self.xmin,self.xmax,self.rgsq,self.eigen_t_files,self.points)
        eigen_dist(self.rgsq)
        print('Done with eigen values distribution')
        print(self.shape_files)
        print(self.shape_files_rg)
        shape_avg(self.shape_files,self.shape_files_rg)

    def plot_msds_avg(self):
        #files=glob.glob(stored_path+'/'+d+'_g1_RUN_*_lambda_%d.dat'%(i))
        print(self.g1_files)
        g1_avg=np.zeros((self.points,2))
        g3_avg=np.zeros((self.points,2))
       
        add_path=os.path.dirname(os.path.abspath(self.g1_files[0]))
        colors = pl.cm.Greens(np.linspace(0,0.8,len(self.g1_files)))

        if self.act_key=='total':
            npart=self.nchains*self*DP
        elif self.act_key=='active': 
            npart=self.nchains*self*DP*how_many
        elif self.act_key=='passive': 
            npart=self.nchains*self*DP*(1-how_many)
        count=0
        for j in self.g1_files: 
            data=np.genfromtxt(j)
            data_interp+=np.interp(x,data[:,0],data[:,1])
            plt.loglog(self.x,data_interp[:,1]/npart,lw=1,linestyle='--',color=colors[count])
            g1_avg[:,1]+=data_interp[:,1]
            count+=1
        g1_avg[:,1]=g1_avg[:,1]/len(files)
        plt.loglog(self.x,g1_avg[:,1]/npart,lw=1,color=trump_color)
        plt.title('$\lambda$=%d/%s'%(self.dipol,self.act_key))
        plt.xlabel(r'Time [$\tau_0$]')
        plt.ylabel(r'$\dfrac{g^1 (t)}{N} [\sigma^2]$')  
        plt.savefig(add_path+'/'+self.act_key+'_g1_lambda_%d.jpg'%(self.dipol),dpi=300,bbox_inches='tight')
#        plt.show()
        plt.clf()
        g1_avg[:,0]=self.x
        np.savetxt(add_path+'/avg_'+self.act_key+'_g1_lambda_%d.avg'%(self.dipol),g1_avg)

        for j in self.g3_files: 
            data=np.genfromtxt(j)
            data_interp+=np.interp(x,data[:,0],data[:,1])
            plt.loglog(self.x,data_interp[:,1]/self.nchains,lw=1,linestyle='--')
            g3_avg[:,1]+=data_interp[:,1]
        g3_avg[:,1]=g3_avg[:,1]/len(files)
        plt.loglog(self.x,g3_avg[:,1]/self.nchains,lw=1,color=trump_color)
        plt.title('$\lambda$=%d/%s'%(self.dipol,self.act_key))
        plt.xlabel(r'Time [$\tau_0$]')
        plt.ylabel(r'$\dfrac{g^3 (t)}{M} [\sigma^2]$')  
        plt.savefig(add_path+'/'+self.act_key+'_g3_lambda_%d.jpg'%(self.dipol),dpi=300,bbox_inches='tight')
#        plt.show()
        plt.clf()
        g3_avg[:,0]=self.x
        np.savetxt(add_path+'/avg_'+self.act_key+'_g3_lambda_%d.avg'%(self.dipol),g3_avg)
    
    def plot_msds_lambda(self):
        count=0
        add_path=os.path.dirname(os.path.abspath(self.g1_files_lambda[0]))
        colors = pl.cm.cividis(np.linspace(0,1,len(self.g1_files_lambda)))
        print('g1 files - lambda:',self.g1_files_lambda)
        for i in self.g1_files_lambda: 
            data=np.genfromtxt(i)
            plt.loglog(data[:,0],data[:,1]/npart,label='$\lambda$=%d'%(count+1),color=colors[count])
            plt.legend(frameon=False,ncol=2)
            plt.title(d)
            plt.xlabel(r'Time [$\tau_0$]')
            plt.ylabel(r'$\dfrac{g^1 (t)}{N} [\sigma^2]$')
            plt.savefig(add_path+'/'+'g1_lambda_dep_%s.png'%self.act_key,dpi=300,bbox_inches='tight')
            count+=1
        #plt.show()
        plt.clf()
        count=0
        print('g3 files - lambda:',self.g3_files_lambda)
        for i in self.g3_files_lambda: 
            data=np.genfromtxt(i)
            plt.loglog(data[:,0],data[:,1]/self.nchains,label='$\lambda$=%d'%(count+1),color=colors[count])
            plt.legend(frameon=False,ncol=2)
            plt.title(d)
            plt.xlabel(r'Time [$\tau_0$]')
            plt.ylabel(r'$\dfrac{g^3 (t)}{M} [\sigma^2]$')
            plt.savefig(add_path+'/'+'g3_lambda_dep_%s.png'%self.act_key,dpi=300,bbox_inches='tight')
            count+=1
        #plt.show()
        plt.clf()


    def plot_msds_three_sp(self):
        legend_dict=['total','active','passive']
        for nam in legend_dict:
            if nam=='total':
                npart=self.nchains*self*DP
            elif nam=='active': 
                npart=self.nchains*self*DP*how_many
            elif nam=='passive': 
                npart=self.nchains*self*DP*(1-how_many)
            for file in self.g1_files_three_sp:
                if (self.msds_path+'/avg'+nam+'_g1_lambda_%d.avg'%self.dipol)==file:
                    print('##################################')
                    print('Choosing legend input')
                    print(nam,file)
                    print('##################################')
                    leg=nam
                    break
            data=np.genfromtxt(file)
            plt.loglog(data[:,0],data[:,1]/npart,legend=nam)
            plt.legend(frameon=False,ncol=2)
            plt.title(r'$\lambda=%d$'%self.dipol)
            plt.xlabel(r'Time [$\tau_0$]')
            plt.ylabel(r'$\dfrac{g^1 (t)}{N} [\sigma^2]$')


    def plot_msids(self):





################################################################################################################
part_dict=['total','active','passive']
frames=9000
warm_frames=2000
run_steps=100000
time_step=5e-3
nchains=1
DP=200
L=80
how_many=0.5
lambdas=np.arange(1,5)
for d in part_dict:
    print('NOW RUNNING',d)
    for i in lambdas:
        print('lambda:',i)
        paths=glob.glob("RUN_*/lambda_%d/polymer_magnetic.vtf"%(i))
        print(paths,len(paths))
        dipol_analyze=analysis_dipol(paths,frames,warm_frames,run_steps,time_step,nchains,DP,L,how_many,i,d)
        dipol_analyze.analyze()