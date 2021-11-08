#### Script for calculating MSD and Rg for each species in a polymer blend #####
## Usage is: python3 msd.py directory-containing-data number-of-polymers DP-small DP-Large##############


import numpy as np 
from matplotlib import pyplot as plt 
import glob 
import os 
import sys
from numba import njit
################## Defining Functions to use while calculating ###############
def CoM (out,npart): #### Calculates CoM of particles takes as input unwrapped data, box and #of particles ############
    comx=np.sum(out[:,0])
    comy=np.sum(out[:,1])
    comz=np.sum(out[:,2])
    com = np.array([comx,comy,comz])/(npart)
    return com

def read_h (files):
    #### Read Header File, getting number of particles, box size, timestep #######
    f=open(files[0],'r')
    trash = f.readline()
    trash = f.readline()
    elems = str.split(trash.strip()," ") 
    #  print (str.isdigit(elems[0]))
    if str.isdigit(elems[0])==True:
        #print ("yay")
        time1 = float(elems[0])/1
        trash = f.readline() 
        trash = f.readline()
        elems = str.split(trash," ") 
        npart = int(elems[0]);
        trash = f.readline()
        trash = f.readline()
        trash = f.readline()
        trash = f.readline()
        elems = str.split(trash," ")
        mybox = float(elems[1])-float(elems[0])
        trash = f.readline()
    #### Read Header File of second file, getting number of particles, box size, timestep #######
    f=open(files[1],'r')
    trash = f.readline()
    trash = f.readline()
    elems = str.split(trash.strip()," ") 
    #  print (str.isdigit(elems[0]))
    if str.isdigit(elems[0])==True:
        #print ("yay")
        time2 = float(elems[0])/1
        trash = f.readline() 
        trash = f.readline()
        elems = str.split(trash," ") 
        npart = int(elems[0]);
        trash = f.readline()
        trash = f.readline()
        trash = f.readline()
        trash = f.readline()
        elems = str.split(trash," ")
        mybox = float(elems[1])-float(elems[0])
        trash = f.readline()
    framestep=int(time2-time1)
    return framestep, npart, mybox


def save():
    file = request.files['xx']
    extension = os.path.splitext(file.filename)[1]

    xx = generate_filename(extension)

    file.save(os.path.join(app.config['UPLOAD_FOLDER'], xx))

def generate_filename(extension):
    xx = str(uuid.uuid4()) + extension
    if os.path.isfile(os.path.join(app.config['UPLOAD_FOLDER'], xx)):
        return generate_filename(extension)
    return xx
@njit(fastmath=True,)
def g1 (data_t,data_0,box):
    tosum=0
    for k in range (data_t[:,0].size):
        xx=(data_t[k,2]+data_t[k,5]*box)-(data_0[k,2]+data_0[k,5]*box)
        yy=(data_t[k,3]+data_t[k,6]*box)-(data_0[k,3]+data_0[k,6]*box)
        zz=(data_t[k,4]+data_t[k,7]*box)-(data_0[k,4]+data_0[k,7]*box)
        magn=np.sqrt(xx**2+yy**2+zz**2)
        tosum+=magn**2
    tosum=tosum/data_t[:,0].size
    return tosum

@njit(fastmath=True,)      
def compute_rg(data,com):
    temp=0
    for i in range (data[:,0].size):
        xx=np.abs(data[i,0]-com[0])
        yy=np.abs(data[i,1]-com[1])
        zz=np.abs(data[i,2]-com[2])
        temp=temp+xx*2+yy**2+zz**2
    temp=temp/data[:,0].size
    #print(temp)
    return temp
################################################# End of Functions Part #####################################################

############### Defining initial parameters for script ###########################
wdir=sys.argv[1]

filename=wdir+"/test_pos*"
files=glob.glob(filename)
files.sort(key=os.path.getmtime)
print("\n".join(files))

nchains=int(sys.argv[2])

DPsm=int(sys.argv[3]) ### Degree of Polymerization of small. input from user
DPlarge=int(sys.argv[4]) ### Degree of Polymerization of large, input from user
framestep,npart,box=read_h(files)
print("Framestep:",framestep)
print("Npart:",npart)
print("Box:",box)

count=1

data=np.loadtxt(files[0],skiprows=9)
data=sorted(data,key= lambda part: part[0] ) ## Sorting data with increasing particle number
data=np.asarray(data)#### Converts data into an array, because it is initially loaded and sorted as a list
data_small=data[np.where(data[:,1]==1)]#### Finds where particle type is 1 (small) ###########
data_large=data[np.where(data[:,1]==2)]#### Finds where particle type is 2 (large) ###########
comp_small=data_small[:,0].size/(data_small[:,0].size+data_large[:,0].size)
nsmall=int(npart*comp_small/DPsm)
nlarge=int((npart-npart*comp_small)/DPlarge)
print('Small:',nsmall,'Large:',nlarge)
g123=np.ones((len(files)-1,5))
g123_s=np.copy(g123)
g123_l=np.copy(g123)
print("Composition small:",100*comp_small)
#####################################################################################################
################### calculating positions of CoM, at time zero ######################################
CoM_s=[]
CoM_l=[]

#### construct lits with CoM coords for small and large rings at t=0 #########################
for j in range (nsmall):
    #print("Calculating small chain",j)
    workdata=data_small[j*DPsm:(j+1)*DPsm,2:5]+data_small[j*DPsm:(j+1)*DPsm,5:8]*box
    com=CoM(workdata,workdata[:,0].size)
    CoM_s.append(com)
for j in range (nlarge):
    #print("Calculating large chain",j)
    workdata=data_large[j*DPlarge:(j+1)*DPlarge,2:5]+data_large[j*DPlarge:(j+1)*DPlarge,5:8]*box
    com=CoM(workdata,workdata[:,0].size)
    CoM_l.append(com)
#####################################################################################################
rg_list=[]
### Starting Calculation of MSD and Rg for small and large rings separately#######
for i in files[1:]:
    g123[count-1,0]=count*framestep
    g123_s[count-1,0]=count*framestep
    g123_l[count-1,0]=count*framestep
    data_1=np.genfromtxt(i,skip_header=9)
    data_1=sorted(data_1,key= lambda part: part[0] ) ## Sorting data with increasing particle number
    data_1=np.asarray(data_1)
    data_small_1=data_1[np.where(data[:,1]==1)]
    data_large_1=data_1[np.where(data[:,1]==2)]
    g3_ts=0
    counter=0
    g3_tl=0
    g1_ts=0
    g1_tl=0
    rg_s=0
    rg_l=0
    print("Working on:",i)
    if nsmall==0:
        print("No Small")
    else:
        for j in range (nsmall):
            workdata=data_small_1[j*DPsm:(j+1)*DPsm,2:5]+data_small_1[j*DPsm:(j+1)*DPsm,5:8]*box
            com=CoM(workdata,workdata[:,0].size)
            g3_ts+=(np.linalg.norm(com-CoM_s[j]))**2
            counter+=1
            temp_rg=(compute_rg(workdata,com))**0.5
            rg_s+=temp_rg
            if count>500:
                rg_list.append(temp_rg)
        g1_ts=g1(data_small_1,data_small,box)

        #np.mean(((data_small_1[:,2:5]+data_small_1[:,5:8]*box)-(data_small[:,2:5]+data_small[:,5:8]*box))**2)
        g3_ts=g3_ts/(counter)
        rg_s=rg_s/counter
        #g2_ts=g1_ts-g3_ts
        g123_s[count-1,1]=g1_ts
        g123_s[count-1,2]=g3_ts
        g123_s[count-1,3]=rg_s
    if nlarge==0:
            print("No Large")
    else:
        counter=0
        for j in range (nlarge):
            #print(" Calculating large chain",j)
            workdata=data_large_1[j*DPlarge:(j+1)*DPlarge,2:5]+data_large_1[j*DPlarge:(j+1)*DPlarge,5:8]*box
            com=CoM(workdata,workdata[:,0].size)
            g3_tl+=(np.linalg.norm(com-CoM_l[j]))**2
            temp_rg=(compute_rg(workdata,com))**0.5
            rg_l+=temp_rg
            #for k in range (DPlarge):
             #   g1_tl+=np.linalg.norm(((data_large_1[k*j,2:5]+data_large_1[k*j,5:8]*box)-(data_large[k*j,2:5]+data_large[k*j,5:8]*box)))**2
            #g1_tl+=np.linalg.norm(((data_large_1[j*DPlarge:(j+1)*DPlarge,2:5]+data_small_1[j*DPlarge:(j+1)*DPlarge,5:8]*box)-(data_large[j*DPlarge:(j+1)*DPlarge,2:5]+data_large[j*DPlarge:(j+1)*DPlarge,5:8]*box)))**2
            #rg_l+=(np.mean((workdata-com)**2))**0.5
            if count>500:
                rg_list.append(temp_rg)
            counter+=1
        g1_tl=g1(data_large_1,data_large,box)
        g3_tl=g3_tl/((counter))
        rg_l=rg_l/counter
        g123_l[count-1,1]=g1_tl
        g123_l[count-1,2]=g3_tl
        g123_l[count-1,3]=rg_l
    count+=1
         #g2_ts=g1_ts-g3_ts
   # g123_l[count,1:4]=np.asarray(tostore_list)
####################################################################3
#### outputting files for small and large rings #######
import os
import time 
t = time.localtime()
timestamp = time.strftime('%b-%d-%Y_%H%M', t)
if not os.path.exists('results'):
    os.makedirs('results')
if nsmall==0:
    print("No Small")

else:
    filename=wdir+"g13_rg_SMALL_%d.dat"%(DPsm)
    if os.path.isfile(filename)==True:
        filename=filename+timestamp
    np.savetxt(filename,g123_s)
if nlarge==0:
    print("No Large")
else:
    if os.path.isfile(filename)==True:
        filename=filename+timestamp
    filename=wdir+"g13_rg_LARGE_%d.dat"%(DPlarge)
    np.savetxt(filename,g123_l)
rg_list=np.asarray(rg_list)
np.savetxt("rg_list.dat",rg_list)