import numpy as np # Library for multi-dimensional arrays and math operations on them 
import pandas as pd # Library for data manipulation and analysis
import os # Interfacing python with the operating system
import matplotlib.pyplot as plt #Plotting library
import time
from matplotlib.ticker import FuncFormatter
import glob
import matplotlib.pylab as pl
import sys
### My plotting style is inputted here #####
import matplotlib as mpl
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
#############################################

biden_color='#2986cc'
trump_color='#cc0000'

os.chdir(sys.argv[1])
directory='plots'
import os
if not os.path.exists(directory):
    os.makedirs(directory)
stored_path='msds_results'
part_dict=['total','active','passive']
shape_params=['b','c','κ']
shape_names=['Asphericity','Acylindricity','Shape' +' '+'Anisotropy']



### g1 total #####
for d in part_dict:
### g1 total #####
    for i in range(1,5):
        files=glob.glob(stored_path+'/'+d+'_g1_RUN_*_lambda_%d.dat'%(i))
        print(files)
        temp=np.genfromtxt(files[0])
        g1_avg=np.zeros(temp[1:].shape)
        time=temp[1:,0]*5e-3
        for j in files: 
            data=np.genfromtxt(j)
            plt.loglog(time,data[1:,1],lw=1,linestyle='--')
            g1_avg[:,1]+=data[1:,1]
        g1_avg[:,1]=g1_avg[:,1]/len(files)
        plt.loglog(time,g1_avg[:,1],lw=1,color=trump_color)
        plt.title('$\lambda$=%d/%s'%(i,d))
        plt.xlabel(r'Time [$\tau_0$]')
        plt.ylabel('$g^1 (t) [\sigma^2]$')  
        plt.savefig(directory+'/'+d+'_g1_lambda_%d.jpg'%(i),dpi=300,bbox_inches='tight')
#        plt.show()
        plt.clf()
        g1_avg[:,0]=time
        np.savetxt(stored_path+'/avg_'+d+'_g1_lambda_%d.dat'%(i),g1_avg)
 
### lambda dep ####
lambdas=np.arange(4)+1
n=len(lambdas)
colors = pl.cm.gist_earth(np.linspace(0,0.8,n))
npart=[200,100,100]
count_ind=0

for d in part_dict:
        files=glob.glob(stored_path+'/avg_'+d+'_g1_lambda_*.dat')
        print(files)
        #count=1
        count=0
        for i in files: 
            data=np.genfromtxt(i)
            plt.loglog(data[:,0],data[:,1]/npart[count_ind],label='$\lambda$=%d'%(count+1),color=colors[count])
            plt.legend(frameon=False,ncol=2)
            plt.title(d)
            plt.xlabel(r'Time [$\tau_0$]')
            plt.ylabel('$g^1 (t)/N [\sigma^2]$')
            count+=1
        plt.savefig(directory+'/'+'g1_lambda_dep_%s.png'%d,dpi=300,bbox_inches='tight')
        plt.clf()
            #count+=1
        #plt.show()
        #plt.clf()
        count_ind+=1
## g1 avg total/passive/active
lambdas=np.arange(4)+1
npart=[200,100,100]
n=len(npart)
colors = pl.cm.Dark2(np.linspace(0,1,n))
for i in lambdas:
    count=0
    for d in part_dict:
            files=glob.glob(stored_path+'/avg_'+d+'_g1_lambda_%d.dat'%(i))
            print(files)
            data=np.genfromtxt(files[0])
            plt.loglog(data[:,0],data[:,1]/npart[count],label='%s'%d,color=colors[count])
            plt.legend(frameon=False,ncol=2)
            plt.title(r'$\lambda=%d$'%i)
            plt.xlabel(r'Time [$\tau_0$]')
            plt.ylabel('$g^1 (t)/N [\sigma^2]$')
            count+=1
    plt.savefig(directory+'/'+'g1_three_species_%d.png'%i,dpi=300,bbox_inches='tight')

            #count+=1
    plt.clf()
    
for d in part_dict:
  
    ### g3 total #####

    for i in range(1,5):
        files=glob.glob(stored_path+'/'+d+'_g3_RUN_*_lambda_%d.dat'%(i))
        #print(files)
        temp=np.genfromtxt(files[0])
        g3_avg=np.zeros(temp[1:].shape)
        time=temp[1:,0]*5e-3
        for j in files: 
            data=np.genfromtxt(j)
            print(j)
            plt.loglog(time,data[1:,1],lw=3,linestyle='--',label=j)
            plt.title(d)
            g3_avg[:,1]+=data[1:,1]
        g3_avg[:,1]=g3_avg[:,1]/len(files)
        plt.loglog(time,g3_avg[:,1],lw=1,color=trump_color)
        plt.title('$\lambda$=%d/%s'%(i,d))
        plt.ylabel(r'$g^3(t) [\sigma^2]$')
        plt.xlabel(r'$Time [\tau_0]$')
        plt.savefig(directory+'/'+d+'_g3_lambda_%d.jpg'%(i),dpi=300,bbox_inches='tight')
        #plt.show()
        plt.clf()
        g3_avg[:,0]=time
        np.savetxt(stored_path+'/avg'+d+'_g3_lambda_%d.dat'%(i),g3_avg)
    

    
    
    
### lambda dep ####
n=len(lambdas)
colors = pl.cm.gist_earth(np.linspace(0,0.8,n))
npart=[200,100,100]
count_ind=0

for d in part_dict:
        files=glob.glob(stored_path+'/avg'+d+'_g3_lambda_*.dat')
        print(files)
        #count=1
        count=0
        for i in files: 
            data=np.genfromtxt(i)
            plt.loglog(data[:,0],data[:,1],label='$\lambda$=%d'%(count+1),color=colors[count])
            plt.legend(frameon=False,ncol=2)
            plt.title(d)
            plt.xlabel(r'Time [$\tau_0$]')
            plt.ylabel('$g^1 (t)/N [\sigma^2]$')
            count+=1
        plt.savefig(directory+'/'+'g3_lambda_dep_%s.png'%d,dpi=300,bbox_inches='tight')
        plt.clf()

            #count+=1
 #       plt.show()
        count_ind+=1
## g3 avg total/passive/active
lambdas=np.arange(4)+1
npart=[200,100,100]
n=len(npart)
colors = pl.cm.Dark2(np.linspace(0,1,n))
stored_path='msds_results'
for i in lambdas:
    count=0
    for d in part_dict:
            files=glob.glob(stored_path+'/avg'+d+'_g3_lambda_%d.dat'%(i))
            print(files)
            data=np.genfromtxt(files[0])
            plt.loglog(data[:,0],data[:,1],label='%s'%d,color=colors[count])
            plt.legend(frameon=False,ncol=2)
            plt.title(r'$\lambda=%d$'%i)
            plt.xlabel(r'Time [$\tau_0$]')
            plt.ylabel('$g^1 (t)/N [\sigma^2]$')
            count+=1
    plt.savefig(directory+'/'+'g3_three_species_%d.png'%i,dpi=300,bbox_inches='tight')

            #count+=1
#    plt.show()
    plt.clf()    



#### gyration tensor total averaging #####
stored_path='gyration_results'
for d in part_dict:
    for i in range(1,5):
        files=glob.glob(stored_path+'/'+d+'_eigen_RUN_*_lambda_%d.dat'%(i))
        print(files)
        temp=np.genfromtxt(files[0])
        g3_avg=np.zeros(temp[1:].shape)
        time=temp[1:,0]*5e-3
        for j in files:
            data=np.genfromtxt(j)
            g3_avg[:,1:]+=data[1:,1:]
        g3_avg[:,1:]=g3_avg[:,1:]/len(files)
        g3_avg[:,0]=time
        np.savetxt(stored_path+'/avg'+d+'_eigen_lambda_%d.dat'%(i),g3_avg)

### gyration tensor average plotting #####
lambdas=np.arange(4)+1
npart=[200,100,100]
n=len(npart)
colors = pl.cm.Set2(np.linspace(0,1,n))
stored_path='gyration_results'

for idx in range (3):
    avg_sh=np.zeros((lambdas.size,3))
    for i in lambdas:
        count=0
        for d in part_dict:
                files=glob.glob(stored_path+'/avg'+d+'_eigen_lambda_%d.dat'%(i))
                print(files)
                data=np.genfromtxt(files[0])
                data[:,1:]=np.sort(data[:,1:],axis=1)
                if idx==0:
                    toplot=data[:,3]-(data[:,1]+data[:,2])/2
                    avg_sh[i-1,count]=toplot.mean()
                elif idx==1:
                    toplot=data[:,2]-data[:,1]
                    avg_sh[i-1,count]=toplot.mean()

                elif idx==2:
                    toplot=-0.5+1.5*(np.sum(data[:,1:]**2,axis=1)/(np.sum(data[:,1:],axis=1)**2))
                    avg_sh[i-1,count]=toplot.mean()

                    
                    
                plt.plot(data[:,0],toplot,label='%s'%d,color=colors[count])
                plt.xscale('log')
                plt.legend(frameon=False,ncol=1)
                plt.title(r'$\lambda=%d$'%i)
                plt.xlabel(r'Time [$\tau_0$]')
                plt.ylabel('$%s$'%(shape_params[idx]))
                plt.savefig(directory+'/'+'_%s_shape_%d.png'%(shape_params[idx],i),dpi=300,bbox_inches='tight')

                count+=1
     #   plt.show()
        plt.clf()
    np.savetxt(stored_path+'/lambda_dep_shape_%s.txt'%(shape_names[idx]),np.column_stack((lambdas,avg_sh)))
    for jj in range (3):
        plt.plot(lambdas,avg_sh[:,jj],marker='o',markerfacecolor='none',markersize=9,color=colors[jj],label='%s'%(part_dict[jj]))
    plt.legend(frameon=False,ncol=1)
    plt.title(r'%s'%shape_names[idx])
    plt.xlabel(r'$\lambda$')
    plt.ylabel('$%s$'%(shape_params[idx]))
    plt.savefig(directory+'/'+'_%s_shape.png'%(shape_names[idx]),dpi=300,bbox_inches='tight')
    #plt.show()
    plt.clf()

### eigen values total #####
stored_path='gyration_results'
count=0
lambdas=np.arange(4)+1
avg_store=np.zeros((lambdas.size,3))
for d in part_dict:
    for i in range(1,5):
        files=glob.glob(stored_path+'/'+d+'_eigen_RUN_*_lambda_%d.dat'%(i))
        print(files)
        temp=np.genfromtxt(files[0])
        g3_avg=np.zeros(temp[1:].shape)
        time=temp[1:,0]*5e-3
        for j in files: 
            data=np.genfromtxt(j)
            plt.loglog(time,data[1:,1],lw=3,linestyle='--',label=j)
            plt.title(d)
            g3_avg[:,1]+=data[1:,1]
        g3_avg[:,1]=g3_avg[:,1]/len(files)
        plt.loglog(time,g3_avg[:,1],lw=1,color=trump_color)
        plt.title(r'$\lambda=%d/%s$'%(i,d))
        plt.xlabel(r'Time [$\tau_0$]')
        plt.ylabel('$R_g (t) [\sigma]$')

        plt.savefig(directory+'/'+d+'_eigen_lambda_%d.jpg'%(i),dpi=300,bbox_inches='tight')
        #plt.show()
        plt.clf()
        g3_avg[:,0]=time
        np.savetxt(stored_path+'/avg'+d+'_eigen_lambda_%d.dat'%(i),g3_avg)
        avg_store[i-1,count]=g3_avg[-100:,1].mean()
    count+=1
np.savetxt(stored_path+'/lambda_dep_rg.txt',np.column_stack((lambdas,avg_store)))


### Radius of gyration vs lambda ###
toplot=np.genfromtxt(stored_path+'/lambda_dep_rg.txt')
for jj in range (3):
    plt.plot(lambdas,toplot[:,jj+1],marker='o',markerfacecolor='none',markersize=9,color=colors[jj],label='%s'%(part_dict[jj]))
    plt.legend(frameon=False,ncol=1)
    plt.xlabel(r'$\lambda$')
    plt.ylabel('$R_g [\sigma]$')
plt.legend(frameon=False)
plt.savefig(directory+'/'+'lambda_dep_rg.png',dpi=300,bbox_inches='tight')
#plt.show()



