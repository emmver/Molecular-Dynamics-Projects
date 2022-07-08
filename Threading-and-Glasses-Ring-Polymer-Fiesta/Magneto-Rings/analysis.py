'''
This script finds and reads in the vtf and vtfk files (see variables filename and files). These files are produced during the simulation runs.
VTF files contain particle positions for each frame while VTK files contain both positions and dipole vectors.

Usage: python3 analysis.py /path/to/files DP

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
    '''
    This function reads the vtf file split frame-wise.
    Each time the function is called, and the f file is not closed the next frame is read. 
    The input needed is the file f, which is opened in the main part of the script, and the number of particles npart.
    Returns an array with the position of each particle, thus with dimensions, (npart,3)
    '''
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
    '''
    Simply calculates the gyration tensor. 
    Input needed is the particles positions.
    @njit in the beginning is from the numba module and is used to make the calculation faster
    '''
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
    r2=r2/(2*N**2) #proper normalization
    #rg2=r2[0,0]**2+r2[1,1]**2+r2[2,2]**2
    return r2 

#################### MSDs Lib ############
@njit(fastmath=True)
def coms_times(pos,pos_coms,DP,nch,frames):
   '''
   This function calculates the center of mass position for each frame in the trajectory.
   Input
   ------------
   1) pos: particle positions
   2) pos_coms: np.zero array with proper dimensions. This is not created in the function so njit works. njit does not work with 
      some numpy functions such as np.zeros
   3) DP: degree of polymerization
   4) nch: number of chains
   5) frames: number of frames
   
   Returns: array with the center of mass coordinates of each chain for each frame.
   '''
   
    for i in range (frames):
        data=pos[i]
        for j in range (nch):
            workdata=data[j*DP:(j+1)*DP]
            comx=np.sum(workdata[:,0])
            comy=np.sum(workdata[:,1])
            comz=np.sum(workdata[:,2])
            com = np.array([comx,comy,comz])/(workdata[:,0].size)
            pos_coms[i][j]=com
    return pos_coms
#### MSD Libs #######

def _autocorr_fft(x):
    """Compute the autocorrelation of 1d signal using FFT.
    Input: particle positions
    """
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
    """Compute the non-FFT part of the MSD in the FFT method.
    Input: particle positions 
    """
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
    
    Input: particle positions
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
    '''
    Compute Mean-Squared Internal Distances (MSID) as a function of the contour length s. 
    The MSID is computed as the squared Euclidean distance between two monomers s-units apart. 
    
    Input
    -----------------
    1) positions: particle positions
    2) start: At which monomer to start the calculation 
    3) stop: At which monomer to stop the calculation
    4) num_rings: number of rings
    5) num_mons: number of monomers per ring 
    
    Returns
    ----------------
    1) The contour segment length s
    2) The averaged-per-particle MSID
    '''
    
    s = np.arange(1, stop-start)
    msid = np.zeros_like(s, dtype=np.float64)
    #print(msid.shape)
    upds = np.zeros_like(s, dtype=np.int64)
    
    for n in range(num_rings):
        r = positions[n*num_mons:(n+1)*num_mons]
        for index, ds in enumerate(s):            
            i = start
            j = start + ds # ds monomers apart: Typically ds takes values from 1 to DP/2
            while j < stop:
                dr = r[j] - r[i]
                msid[index] += dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]
                upds[index] += 1    #tracks how many times the calculation was done
                j += 1
                i += 1 #moving along the chain 
    return s, msid/upds #averaging



### Dipole Angles #####

@njit(fastmath=True)
def dipol_angles(dipoles,start,stop,num_rings,num_mons):
     '''
    Compute the cosinus of the angle between dipoles as a function of the contour length s. 
    This is computed from the dot product of the dipole vectors of a pair of particles, 
    separated by s-units along the contour.
    
    Input
    -----------------
    1) dipoles: dipole vectors from VTK files
    2) start: At which monomer to start the calculation 
    3) stop: At which monomer to stop the calculation
    4) num_rings: number of rings
    5) num_mons: number of monomers per ring 
    
    Returns
    ----------------
    1) The contour segment length s
    2) The averaged-per-particle dipole cosinus
    '''
    s = np.arange(1, stop-start)
    cosin = np.zeros_like(s, dtype=np.float64)
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
                cosin[index] += dotp #np.arccos(dotp)*180/np.pi
                upds[index] += 1    
                j += 1
                i += 1
    return s, (cosin/upds)
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


os.chdir(sys.argv[1]) #change working directory
###### creating folders for storing resulsts ######
directory='msds_results'                           
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
##################################################

frames_in=10000 # frames in trajectory 
skipframes=3000 # how many to skip at first 
run_steps=10000 #
nchains=1
nch=nchains #number of chains
DP=int(sys.argv[3]) # Degree of polymerization



########### This will later on be used when reading VTK files ################
if DP==100:
    dipol_skiplines1=315
    dipol_skiplines2=315
elif DP==200:
    dipol_skiplines1=614
    dipol_skiplines2=412
###############################################################################

Dpol=DP
L=80
box=L
npart=int(nchains*DP)

##### number of bonds for the case of linear or ring #########
#### Used while reading the VTF files #####
if sys.argv[2]=='ring':
    bonds=npart
elif sys.argv[2]=='linear':
    bonds=npart-1
angles=0
how_many=0.5 #ratio of magnetic/total particles
part_dict=['active','total','passive'] #dictionary to calculate properties for differnt parts of the chain
lambdas_w=[1,2,4,6,9,12,16] # list of lambdas

''' 
Following loops run for each object in the part_dict list i.e. passive, active or total
for each configuration (up to now 8 independent runs) 
for each frame in each configuration
and does the following:

1) Skip non-needed lines and frames from the file read. 
2) Calculate particle MSD 
3) Calculate Center of Mass (CoM) MSD
4) Calculate gyration tensor
5) Mean squared internal distances 
'''
for d in part_dict: # passive active or total
    print('NOW RUNNING',d)
    for i in lambdas_w: # which lambda 
        print('lambda:',i)
        #files_list=[]
        paths=glob.glob("RUN_*/lambda_%d/polymer_magnetic.vtf"%(i)) #reads files for a specific lambda
        print(paths,len(paths))

    ##### Calculate total MSDs #############
        run_count=1
        for j in range(len(paths)):# which run
            print(f"Working on run {run_count} and lambda {i}")
            ### Calculating Pos Array and MSDs for each RUN ####
            file=paths[j]
            print("Working on file",file)
            f=open(file,'r')
            if d==part_dict[0] or d==part_dict[2]: #choose particles active or passive
                    pos=np.zeros((int(frames_in),int(npart*how_many),3)) #create array to store positions of particles
            else: # or all
                pos=np.zeros((int(frames_in),(npart),3)) #create array to store positions of particles
            for ii in range (int(npart)+bonds+angles+2): #skip structural info in VTF file
                    print(f.readline())
            print(f.readline())

            print('skipframes')
            for ii in range (skipframes): #skipping initial frames
                data=read_traj(f,npart)
            count_frames=0
            start=timer()
            frame_num=0
            gg=open('frames_list.txt','w')
            while 1: 
                try: 
                    g="Frame:%d"%(count_frames+1)
                    sys.stdout.write("\r" + str(g))
                    sys.stdout.flush() ## printing info on frame 
                    data=read_traj(f,npart) #reading frame
                    if d==part_dict[1]: #total
                        pos[count_frames]=data[:,1:]
                    elif d==part_dict[0]: #active
                        pos[count_frames]=data[:int(npart*how_many),1:] #which particles to keep
                    elif d==part_dict[2]: #passive
                        pos[count_frames]=data[int(npart*how_many):,1:] #which particles to keep
                    count_frames+=1 #update frame index
                except Exception as e: #stops reading when there is an exception
                    print("Except:",e)
                    print('stopped at',count_frames)
                    pos=pos[:count_frames,:,:] ## keeps non-zero position vectors.
                    gg.write('%d'%count_frames)
                    #gg.write('\n')
                    break

            frames=count_frames
            #mean_pos=pos.mean(axis=1).reshape(-1, 1, 3) #this is useful when more than 1 ring is in the box. Substracts CoM of the box to correct for drift.
            #pos-=mean_pos 
            end=timer()
            print('Wall Time:',end-start)
            time = np.zeros(frames) # simulation time array 
            count=0
            for tt in range(frames): # fill time array 
                time[tt]=count*run_steps
                count+=1
            msds = np.zeros((npart, frames)) # particle msd array
            #print(msds.shape)
            for ii in range(pos.shape[1]): #calculate msd for each particle
                start=timer()
                msds[ii] = self_fft(pos[:,ii])
                end=timer()
                g="Particle %d Wall time %.4f seconds"%(ii, end-start)
            #     sys.stdout.write("\r" + str(g))
            #     sys.stdout.write('\n')
            #     sys.stdout.flush()
            # #       np.savetxt('./'+file[:-4]+'_g1_all.dat',(np.stack((msds),axis=-1)))
            np.savetxt(directory+'/'+d+'_g1_RUN_%d_lambda_%d.dat'%(run_count,i),(np.stack((time[:],msds.mean(axis=0)[:]),axis=-1))) #save MSD for this run and lambda
            nch=int(nchains)
            #pos+=mean_pos
            pos_coms=np.zeros((frames,nch,3))
            pos_coms=coms_times(pos,pos_coms,Dpol,nch,frames) # create array for polymer CoM positions.
            print('Shape POS',pos_coms.shape) 
            g3s = np.zeros((nch,frames)) # CoM MSD array
            for ii in range(pos_coms.shape[1]):
                start=timer()
                g3s[ii] = self_fft(pos_coms[:,ii]) # calculate CoM MSD
                end=timer()
            #     g="Chain %d Wall time %.4f seconds"%(ii, end-start)
            #     sys.stdout.write("\r" + str(g))
            #     sys.stdout.write('\n')
            #     sys.stdout.flush()
            # #print("================= COMs g3 =================")
            print('G3 shape',g3s[0].shape)
            print('time shape',time.shape)
            np.savetxt(directory+'/'+d+'_g3_RUN_%d_lambda_%d.dat'%(run_count,i),(np.stack((time[:],g3s[0]),axis=-1))) #save MSD for this run and lambda
            


            rgs_tens= np.zeros((nchains, frames,3,3)) # create array for gyratio tensor (3x3) for each chain and frame.
            for ii in tqdm(range(0,frames),desc='Gyration Tensor Calculation'):
                    g="Frame:%d"%(ii+1)
                    sys.stdout.write("\r" + str(g))
                    sys.stdout.flush()
                   # sys.stdout.write(" ")  
                    sys.stdout.flush()

                    for j in range (nchains): #calculate gyration tensor for each chain
                            ##g="Chain:%d"%(j+1)
                            #sys.stdout.write("\r" + str(g))
                            #sys.stdout.flush()
                            temp_data=pos[ii,j*DP:(j+1)*DP]
                            rgs_tens[j,ii]=rg_tens(temp_data)
            radii=[]
            radii_all=np.zeros((frames,nchains)) #stores rg values
            rg_tensor_eigen=np.zeros((frames,3)) #stores eigenvalues 
            for ii in tqdm(range(0,frames),desc='Gyration Tensor Diag'): #Properly diagonalize gyration tensor
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
                radii.append((temp_rg)**0.5) 
                rg_tensor_eigen[ii]=temp_eigen
  
         
            if j==0:
                j=1
            radii=np.asarray(radii)/j #averaged for chains
            rg_tensor_eigen=rg_tensor_eigen/j#averaged for chains 
            corr=np.stack((time,radii),axis=-1)
            _M =corr[:].shape[0]
            acf, ci = sm.tsa.acf(radii[:], alpha=0.05,nlags=_M,fft=True) # autocorrelation of Rg
            tacf=np.arange(_M)*1000 #times for autocorrelation
            
            ### save all results for this lambda and this run####
            np.savetxt(gyr_directory+'/'+d+'_tsa_RUN_%d_lambda_%d.dat'%(run_count,i),(np.stack((tacf,acf),axis=-1)))
            np.savetxt(gyr_directory+'/'+d+'_rg_RUN_%d_lambda_%d.dat'%(run_count,i),(np.stack((time[:],radii),axis=-1)))
            np.savetxt(gyr_directory+'/'+d+'_eigen_RUN_%d.dat'%(run_count),(np.column_stack([time[:],rg_tensor_eigen])))
            
            ### Internal Distances #####
            #DP=pos.shape[1]
            if d=='total': #for the whole chain
                start=0
                stop=int(DP/2)
                eff_DP=DP
                msid_dist=np.zeros((np.arange(0,int(stop-start)).size,frames+1))
                for ii in tqdm(range(0,frames),desc='Internal Distances'):
                    msid_dist[:-1,0],msid_dist[:-1,ii]=compute_msid(pos[ii,:], start, stop, nch, eff_DP)
            elif d=='active': #active part
                start=0
                stop=int(DP*how_many)
                eff_DP=int(DP/2) #effective degree of polymerization
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
            elif d=='passive': #passive part
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






