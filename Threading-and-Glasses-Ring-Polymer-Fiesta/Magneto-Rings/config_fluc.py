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
from tqdm import tqdm
from mpl_toolkits import mplot3d
### My plotting style is inputted here #####
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 150
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
#############################################

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

os.chdir(sys.argv[1])


frames=10000
nchains=1
DP=200
Dpol=DP
L=80
box=L
npart=int(nchains*DP)
if sys.argv[2]=='lin':
    bonds=npart-1
elif sys.argv[2]=='rin':
    bonds=npart
angles=0
how_many=0.5 
part_dict=['total','active','passive']
early_frames=1000
f_step_e=20
f_step_l=100
file_idx=np.random.randint(6)+1
print(file_idx)
store_folder='configurations_se_%d_sl_%d'%(f_step_e,f_step_l)
if os.path.isdir(store_folder)==False:
    os.mkdir(store_folder)

for d in part_dict:
    print('NOW RUNNING',d)
    for ii in range(1,5):
        print('lambda:',ii)
        #files_list=[]
        paths=glob.glob("RUN_%d/lambda_%d/polymer_magnetic.vtf"%(file_idx,ii))
        print(paths)
        file=paths[0]
        f=open(file,'r')
        for hh in range (int(npart)+bonds+angles+2):
            f.readline()
        print(f.readline())
        print(paths,len(paths))
        pos=np.zeros((int(frames),(npart),3))
        count_frames=0

        ##creating positions array
        ggg=open('frames_list.txt','w')
        while 1: 
            try: 
                g="Frame:%d"%(count_frames+1)
                sys.stdout.write("\r" + str(g))
                sys.stdout.flush()
                data=read_traj(f,npart)
                pos[count_frames]=data[:,1:]
                count_frames+=1
            except Exception as e:
                print("Except:",e)
                print('stopped at',count_frames)
                ggg.write('%d'%count_frames)
                ggg.write('\n')
                pos=pos[:count_frames,:,:]
                break


        f_step=f_step_e
        c_act = mpl.cm.Reds(np.linspace(0,0.8,int(early_frames/f_step)+1))
        c_pass=mpl.cm.Blues(np.linspace(0,0.8,int(early_frames/f_step)+1))
        fig = plt.figure(figsize=(9,16))
        ax = fig.add_subplot(111, projection='3d')
        
        for i in range(0,early_frames,f_step): 
            g="Frame:%d"%(count_frames+1)
            sys.stdout.write("\r" + str(g))
            sys.stdout.flush()
            data=pos[i,:,:]
            toplot=data[:,:];#print('Shape:',toplot.shape)
            mean_pos=toplot.mean(axis=0).reshape(1, 3)
            #print('Mean Pos',mean_pos)
            toplot-=mean_pos
            ax.view_init(azim=0, elev=90)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_zticks([])
            ax.grid(False)
            ax.axis("off")
            ax.xaxis.pane.fill = False
            ax.yaxis.pane.fill = False
            ax.zaxis.pane.fill = False
            if d==part_dict[0]:
                ax.plot(toplot[:int(npart*how_many),0],toplot[:int(npart*how_many),1],toplot[:int(npart*how_many),2],color=c_act[int(i/f_step)],lw=1)
                ax.plot(toplot[int(npart*how_many):npart,0],toplot[int(npart*how_many):npart,1],toplot[int(npart*how_many):npart,2],color=c_pass[int(i/f_step)],lw=1)
            elif d==part_dict[1]:
                ax.plot(toplot[:int(npart*how_many),0],toplot[:int(npart*how_many),1],toplot[:int(npart*how_many),2],color=c_act[int(i/f_step)],lw=1)
            elif d==part_dict[2]:
                ax.plot(toplot[int(npart*how_many):npart,0],toplot[int(npart*how_many):npart,1],toplot[int(npart*how_many):npart,2],color=c_pass[int(i/f_step)],lw=1)
            # if i==early_frames-1 or i==early_frames-2:
            #     print('early last')
            #     if sys.argv[2]=='rin':
            #         xyz=np.array([[toplot[0,0],toplot[0,1],toplot[0,2]],
            #         [toplot[-1,0],toplot[-1,1],toplot[-1,2]]])
            #         ax.plot(xyz[:,0],xyz[:,1],xyz[:,2],'k-',lw=1)
            #     ax.plot(toplot[:,0],toplot[:,1],toplot[:,2],'k-',lw=1)
           ## elif i==0:
             #   ax.plot(toplot[:,0],toplot[:,1],toplot[:,2],color='#949392',lw=1)
            count_frames+=1
        data=pos[early_frames-1,:,:]
        toplot=data[:,:]
        mean_pos=toplot.mean(axis=0).reshape(1, 3)
        toplot-=mean_pos
        ax.plot(toplot[:,0],toplot[:,1],toplot[:,2],'k-',lw=1)
        if sys.argv[2]=='rin':
            xyz=np.array([[toplot[0,0],toplot[0,1],toplot[0,2]],
            [toplot[-1,0],toplot[-1,1],toplot[-1,2]]])
            ax.plot(xyz[:,0],xyz[:,1],xyz[:,2],'k-',lw=1)
        plt.close(fig)
        fig.savefig(store_folder+'/early_frames_%s_lambda_%d.jpg'%(d,ii),dpi=300,bbox_inches='tight')
        plt.clf()
        ax.clear()
        f_step=f_step_l
        c_act = mpl.cm.Reds(np.linspace(0,0.8,int(frames/(2*f_step)+1)))
        c_pass=mpl.cm.Blues(np.linspace(0,0.8,int(frames/(2*f_step)+1)))
        for i in range(0,int(frames/2),f_step): 
            g="Frame:%d"%(i)
            sys.stdout.write("\r" + str(g))
            sys.stdout.flush()
            data=pos[i+early_frames,:,:]
            toplot=data[:,:]
            mean_pos=toplot.mean(axis=0).reshape(1, 3)
            #print('Mean Pos',mean_pos)
            toplot-=mean_pos
            ax.view_init(azim=0, elev=90)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_zticks([])
            ax.grid(False)
            ax.axis("off")
            ax.xaxis.pane.fill = False
            ax.yaxis.pane.fill = False
            ax.zaxis.pane.fill = False
            ax.set_title(r'$\lambda=%d$'%ii,fontsize=25)
            if d==part_dict[0]:
                ax.plot(toplot[:int(npart*how_many),0],toplot[:int(npart*how_many),1],toplot[:int(npart*how_many),2],color=c_act[int(i/f_step)],lw=1)
                ax.plot(toplot[int(npart*how_many):npart,0],toplot[int(npart*how_many):npart,1],toplot[int(npart*how_many):npart,2],color=c_pass[int(i/f_step)],lw=1)
            elif d==part_dict[1]:
                ax.plot(toplot[:int(npart*how_many),0],toplot[:int(npart*how_many),1],toplot[:int(npart*how_many),2],color=c_act[int(i/f_step)],lw=1)
            elif d==part_dict[2]:
                ax.plot(toplot[int(npart*how_many):npart,0],toplot[int(npart*how_many):npart,1],toplot[int(npart*how_many):npart,2],color=c_pass[int(i/f_step)],lw=1)
            # if i==int(frames/2)-1 or i==int(frames/2)-2:
            #     print('last',i+early_frames)
            #     ax.plot(toplot[:,0],toplot[:,1],toplot[:,2],'k-',lw=1)
            #     if sys.argv[2]=='rin':
            #         xyz=np.array([[toplot[0,0],toplot[0,1],toplot[0,2]],
            #         [toplot[-1,0],toplot[-1,1],toplot[-1,2]]])
            #         ax.plot(xyz[:,0],xyz[:,1],xyz[:,2],'k-',lw=1)
            #elif i==0:
             #   ax.plot(toplot[:,0],toplot[:,1],toplot[:,2],color='#949392',lw=1)
            count_frames+=1
        data=pos[int(frames/2)+early_frames-1,:,:]
        toplot=data[:,:]
        mean_pos=toplot.mean(axis=0).reshape(1, 3)
        toplot-=mean_pos
        ax.plot(toplot[:,0],toplot[:,1],toplot[:,2],'k-',lw=1)
        if sys.argv[2]=='rin':
            xyz=np.array([[toplot[0,0],toplot[0,1],toplot[0,2]],
            [toplot[-1,0],toplot[-1,1],toplot[-1,2]]])
            ax.plot(xyz[:,0],xyz[:,1],xyz[:,2],'k-',lw=1)
            #elif i==0:
        plt.close(fig)
        fig.savefig(store_folder+'/late_frames_%s_lambda_%d.jpg'%(d,ii),dpi=300,bbox_inches='tight')
        plt.clf()
