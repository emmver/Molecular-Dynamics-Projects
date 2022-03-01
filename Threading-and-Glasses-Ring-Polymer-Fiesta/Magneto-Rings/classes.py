import numpy as np 
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
from utils import *


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
    def __init__(self,path,frames,warm_frames,run_steps,nchains,DP,L,how_many,dipol,arch_key):
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
        pbar.close()

    def calc_g1(self):
        msds = np.zeros((self.npart, self.pos.shape[0]))
        for ii in tqdm(range(self.pos.shape[1]),desc='Monomer MSDs'):
            msds[ii] = self_fft(self.pos[:,ii])
        np.savetxt(self.directory+'/'+self.act_key+'_g1_RUN_%d_lambda_%d.dat'%(self.run_count,self.dipol),(np.stack((self.time[:],msds.mean(axis=0)[:]),axis=-1)))
        
    def calc_g3(self):
        pos_coms=np.zeros(self.pos.shape[0],self.nchains,3))
        pos_coms=coms_times(self.pos,pos_coms,self.DP,self.nchains,self.pos.shape[0])
        g3s = np.zeros((self.nchains,self.pos.shape[0]))
        for ii in tqdm(range(pos_coms.shape[1]),desc='CoM MSDs'):
            g3s[ii] = self_fft(pos_coms[:,ii])
        np.savetxt(self.directory+'/'+self.act_key+'_g3_RUN_%d_lambda_%d.dat'%(self.run_count,self.dipol),(np.stack((self.time[:],g3s[0]),axis=-1)))

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
class plot_equil():


################################################################################################################

class plot_dipol():


################################################################################################################
