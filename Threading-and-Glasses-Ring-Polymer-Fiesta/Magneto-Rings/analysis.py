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
Otherwise the provided plostyle.mpl file should be used οr created from scratch. 
'''





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
from sympy import Matrix
from tqdm import tqdm
import statsmodels.api as sm
####################################################################
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


######### Calculate Gyration Tensor #######################3
@njit(fastmath=True)
def rg_tens (r):
    
    N = r.shape[0]
    r2=np.zeros((3,3))
    for i in range (N):
        for j in range (N):
            r2[0,0]+=(r[i,0]-r[j,0])*(r[i,0]-r[j,0])#diag
            r2[1,0]+=(r[i,1]-r[j,1])*(r[i,0]-r[j,0])#ndiag
            r2[2,0]+=(r[i,2]-r[j,2])*(r[i,0]-r[j,0])#ndiag
            r2[0,1]+=(r[i,0]-r[j,0])*(r[i,1]-r[j,1])#ndiag
            r2[0,2]+=(r[i,0]-r[j,0])*(r[i,2]-r[j,2])#ndiag
            r2[1,1]+=(r[i,1]-r[j,1])*(r[i,1]-r[j,1])#diag
            r2[2,1]+=(r[i,2]-r[j,2])*(r[i,1]-r[j,1])#ndiag
            r2[1,2]+=(r[i,1]-r[j,1])*(r[i,2]-r[j,2])#ndiag
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
            #print("Calculating small chain",j)
            workdata=data_small[j*Dpolsm:(j+1)*Dpolsm]
            #print(workdata.shape)
            comx=np.sum(workdata[:,0])
            comy=np.sum(workdata[:,1])
            comz=np.sum(workdata[:,2])
            com = np.array([comx,comy,comz])/(workdata[:,0].size)
            #print('COM is:',com)
            #com=CoM(workdata,workdata[:,0].size)
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
    #print(msid.shape)
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



### Dipole Angles #####

@njit(fastmath=True)
def dipol_angles(dipoles,start,stop,num_rings,num_mons):
    s = np.arange(1, stop-start)
    msid = np.zeros_like(s, dtype=np.float64)
    #print(msid.shape)
    upds = np.zeros_like(s, dtype=np.int64)
    
    for n in range(num_rings):
        r = dipoles[n*num_mons:(n+1)*num_mons]
        for index, ds in enumerate(s):            
            i = start
            j = start + ds
            while j < stop:
                dotp=np.dot(r[j]/np.linalg.norm(r[j]),r[i]/np.linalg.norm(r[i]))
                #dr = r[j] - r[i]
                msid[index] += dotp #np.arccos(dotp)*180/np.pi
                upds[index] += 1    
                j += 1
                i += 1
    return s, (msid/upds)
###################################################################



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
os.chdir(sys.argv[1])
directory='msds_results'
import os
if not os.path.exists(directory):
    os.makedirs(directory)
gyr_directory='gyration_results'

if not os.path.exists(gyr_directory):
    os.makedirs(gyr_directory)

etoe_directory='endtoend_results'

if not os.path.exists(etoe_directory):
    os.makedirs(etoe_directory)
msid_dir='msid_results'
if not os.path.exists(msid_dir):
    os.makedirs(msid_dir)

frames_in=10000
skipframes=3000
run_steps=10000
nchains=1
nch=nchains
DP=int(sys.argv[3])

if DP==100:
    dipol_skiplines1=315
    dipol_skiplines2=315
elif DP==200:
    dipol_skiplines1=614
    dipol_skiplines2=412


Dpol=DP
L=80
box=L
npart=int(nchains*DP)
if sys.argv[2]=='ring':
    bonds=npart
elif sys.argv[2]=='linear':
    bonds=npart-1
angles=0
how_many=0.5 
part_dict=['active','total','passive']
#part_dict=['active']
lambdas_w=[1,2,4,6,9,12,16]
#lambdas_w=[1,4,9,16]

for d in part_dict:
    print('NOW RUNNING',d)
    for i in lambdas_w:
        print('lambda:',i)
        #files_list=[]
        paths=glob.glob("RUN_*/lambda_%d/polymer_magnetic.vtf"%(i))
        print(paths,len(paths))

    ##### Calculate total MSDs #############
        run_count=1
        for j in range(len(paths)):
            print(f"Working on run {run_count} and lambda {i}")
            ### Calculating Pos Array and MSDs for each RUN ####
            file=paths[j]
            print("Working on file",file)
            f=open(file,'r')
            if d==part_dict[0] or d==part_dict[2]:
                    pos=np.zeros((int(frames_in),int(npart*how_many),3))
            else:
                pos=np.zeros((int(frames_in),(npart),3))
            for ii in range (int(npart)+bonds+angles+2):
                    print(f.readline())
            print(f.readline())

            print('skipframes')
            for ii in range (skipframes):
                data=read_traj(f,npart)
            count_frames=0
            start=timer()
            frame_num=0
            gg=open('frames_list.txt','w')
            while 1: 
                try: 
                    g="Frame:%d"%(count_frames+1)
                    sys.stdout.write("\r" + str(g))
                    sys.stdout.flush()
                    data=read_traj(f,npart)
                    if d==part_dict[1]:
                        pos[count_frames]=data[:,1:]
                    elif d==part_dict[0]:
                        pos[count_frames]=data[:int(npart*how_many),1:]
                    elif d==part_dict[2]:
                        pos[count_frames]=data[int(npart*how_many):,1:]
                    count_frames+=1
                except Exception as e:
                    print("Except:",e)
                    print('stopped at',count_frames)
                    pos=pos[:count_frames,:,:]
                    gg.write('%d'%count_frames)
                    #gg.write('\n')
                    break

            frames=count_frames
            mean_pos=pos.mean(axis=1).reshape(-1, 1, 3)
            #pos-=mean_pos
            end=timer()
            print('Wall Time:',end-start)
            time = np.zeros(frames)
            count=0
            for tt in range(frames):
                time[tt]=count*run_steps
                count+=1
            msds = np.zeros((npart, frames))
            #print(msds.shape)
            for ii in range(pos.shape[1]):
                start=timer()
                msds[ii] = self_fft(pos[:,ii])
                end=timer()
                g="Particle %d Wall time %.4f seconds"%(ii, end-start)
            #     sys.stdout.write("\r" + str(g))
            #     sys.stdout.write('\n')
            #     sys.stdout.flush()
            # #       np.savetxt('./'+file[:-4]+'_g1_all.dat',(np.stack((msds),axis=-1)))
            np.savetxt(directory+'/'+d+'_g1_RUN_%d_lambda_%d.dat'%(run_count,i),(np.stack((time[:],msds.mean(axis=0)[:]),axis=-1)))
            nch=int(nchains)
            #pos+=mean_pos
            pos_coms=np.zeros((frames,nch,3))
            pos_coms=coms_times(pos,pos_coms,Dpol,nch,frames)
            np.savetxt('pos_coms.txt',pos_coms[:,0])
            print('Shape POS',pos_coms.shape) 
            #pos=np.zeros((len(files),npart,3))
            g3s = np.zeros((nch,frames))
            for ii in range(pos_coms.shape[1]):
                start=timer()
                g3s[ii] = self_fft(pos_coms[:,ii])
                end=timer()
            #     g="Chain %d Wall time %.4f seconds"%(ii, end-start)
            #     sys.stdout.write("\r" + str(g))
            #     sys.stdout.write('\n')
            #     sys.stdout.flush()
            # #print("================= COMs g3 =================")
            print('G3 shape',g3s[0].shape)
            print('time shape',time.shape)
            np.savetxt(directory+'/'+d+'_g3_RUN_%d_lambda_%d.dat'%(run_count,i),(np.stack((time[:],g3s[0]),axis=-1)))
            


            rgs_tens= np.zeros((nchains, frames,3,3))
            for ii in tqdm(range(0,frames),desc='Gyration Tensor Calculation'):
                    g="Frame:%d"%(ii+1)
                    sys.stdout.write("\r" + str(g))
                    sys.stdout.flush()
                   # sys.stdout.write(" ")  
                    sys.stdout.flush()

                    for j in range (nchains):
                            ##g="Chain:%d"%(j+1)
                            #sys.stdout.write("\r" + str(g))
                            #sys.stdout.flush()
                            temp_data=pos[ii,j*DP:(j+1)*DP]
                            rgs_tens[j,ii]=rg_tens(temp_data)
            #(rgs[0,500,0,0]**2+rgs[0,500,1,1]**2+rgs[0,500,2,2]**2)**0.5
            radii=[]
            radii_all=np.zeros((frames,nchains))
            rg_tensor_eigen=np.zeros((frames,3))
            for ii in tqdm(range(0,frames),desc='Gyration Tensor Diag'):
                M=Matrix(rgs_tens[:,ii][0][0:3])
                P,D=M.diagonalize()
                tocalc=np.sort(np.diag(D))
                temp_rg=0
                temp_rs=[]
                temp_eigen=np.array([0,0,0])
                for j in range (0,nchains):
                    temp=(tocalc[0])+(tocalc[1])+(tocalc[2])
                    temp_rg+=temp
                    temp_eigen=temp_eigen+tocalc
                    radii_all[ii,j]=temp**0.5
                radii.append(temp_rg**0.5)
                rg_tensor_eigen[ii]=temp_eigen
            #print("\n Radii")
            #print(radii)
         
            if j==0:
                j=1
            radii=np.asarray(radii)/j
            rg_tensor_eigen=rg_tensor_eigen/j
            corr=np.stack((time,radii),axis=-1)
            _M =corr[:].shape[0]
            acf, ci = sm.tsa.acf(radii[:], alpha=0.05,nlags=_M,fft=True)
            tacf=np.arange(_M)*1000
            np.savetxt(gyr_directory+'/'+d+'_tsa_RUN_%d_lambda_%d.dat'%(run_count,i),(np.stack((tacf,acf),axis=-1)))
            np.savetxt(gyr_directory+'/'+d+'_rg_RUN_%d_lambda_%d.dat'%(run_count,i),(np.stack((time[:],radii),axis=-1)))
            np.savetxt(gyr_directory+'/'+d+'_eigen_RUN_%d.dat'%(run_count),(np.column_stack([time[:],rg_tensor_eigen])))
            
            ### Internal Distances #####
            #DP=pos.shape[1]
            if d=='total':
                start=0
                stop=int(DP/2)
                eff_DP=DP
                msid_dist=np.zeros((np.arange(0,int(stop-start)).size,frames+1))
                for ii in tqdm(range(0,frames),desc='Internal Distances'):
                    msid_dist[:-1,0],msid_dist[:-1,ii]=compute_msid(pos[ii,:], start, stop, nch, eff_DP)
            elif d=='active':
                start=0
                stop=int(DP*how_many)
                eff_DP=int(DP/2)
                msid_dist=np.zeros((np.arange(0,int(stop-start)).size,frames+1))
                temp_dipol_files=glob.glob(os.path.dirname(file)+'/vtk_equil/part_test_*.vtk')
                #print(dipol_files)
                pop_idx=np.arange(500,len(temp_dipol_files))
                dipol_files=[]
                for pop in pop_idx:
                    dipol_files.append(glob.glob(os.path.dirname(file)+'/vtk_equil/part_test_%d.vtk'%pop)[0])
                print('Files:',len(dipol_files))
               # sys.exit()
                dipol_dist=np.zeros((np.arange(0,int(stop-start)).size,len(dipol_files)+1))
                print(dipol_files[0])
                for ii in tqdm(range(0,frames),desc='Internal Distances'):
                    msid_dist[:-1,0],msid_dist[:-1,ii]=compute_msid(pos[ii,:], start, stop, nch, eff_DP)
                for dip_count in tqdm(range(0,len(dipol_files)),desc='Dipole Angles'):
                    if i==2 or i==6 or i==12:
                        dipoles_in=np.genfromtxt(dipol_files[dip_count],skip_header=dipol_skiplines1)#614)
                    else:
                        dipoles_in=np.genfromtxt(dipol_files[dip_count],skip_header=dipol_skiplines2)#412)
                    del_list=[]
                #    for del_idx, i in enumerate(dipoles_in):
                        #print(i)
                        #if np.linalg.norm(i)==0:
                            #del_list.append(del_idx)
                            #dipoles=np.delete(dipoles_in,del_idx)
                    
                    #if i==2:
                    dip_idx=[]   
                    for dip in range(dipoles_in[:,0].size):
                        if np.linalg.norm(dipoles_in[dip])==0:
                            #print('Deleted:',dip)
                            #print('----------------------')i
                            dip_idx.append(dip)
                    dipoles=np.delete(dipoles_in,dip_idx,axis=0)
                    #print(dipoles.shape)
                        #print(dipoles[dip])    
                        #print(dip,np.linalg.norm(dipoles[dip])**2)
                   # print(dipoles[:eff_DP])
                    dipol_dist[:-1,0],dipol_dist[:-1,dip_count]=dipol_angles(dipoles,start,stop,nch,eff_DP)
                dipol_dist[:-1,-1]=np.mean(dipol_dist[:-1,1:-1],axis=1)
                np.savetxt(msid_dir+'/'+d+'dipol_dist_RUN_%d_lambda_%d.dat'%(run_count,i),dipol_dist[:-1,:])
            elif d=='passive':
                start=0
                stop=int(DP*(1-how_many))
                eff_DP=int(DP/2)
                msid_dist=np.zeros((np.arange(0,int(stop-start)).size,frames+1))
                for ii in tqdm(range(0,frames),desc='Internal Distances'):
                    msid_dist[:-1,0],msid_dist[:-1,ii]=compute_msid(pos[ii,:], start, stop, nch, eff_DP)

            msid_dist[:-1,-1]=np.mean(msid_dist[:-1,1:-1],axis=1)
            
            np.savetxt(msid_dir+'/'+d+'msids_RUN_%d_lambda_%d.dat'%(run_count,i),msid_dist[:-1,:])
            
            run_count+=1
            f.close()






