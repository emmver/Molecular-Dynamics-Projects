import numpy as np
import sys
import os 
import string
import glob

sys.path.insert(1,"./ANALYSIS/")

from utilities import *


class mypolymer:

    def __init__(self, length):
        self.length = length
        self.RoG = 0.
        self.pol_shape = np.zeros(2)
        self.Ree2 = []
        self.Ree = []

        self.bondcorr = np.zeros(length-1,dtype = float)
        self.persist_Ree = np.zeros(length-1,dtype = float)
        self.blenavg = []

        self.REcorr = np.zeros(1,dtype = float)
        self.RErep = np.zeros(1,dtype = float)
        self.bondhist = np.zeros((50,2), dtype = float)
	
        self.rhist = np.zeros((100,2),dtype = float) 
        self.thhist = np.zeros((100,2),dtype = float)

        self.lp = []


if len(sys.argv) != 6:
    print("usage is", sys.argv[0], "folder gen Lbackbone lp_density nskip")
    exit(1)

myfold = sys.argv[1]
mygen = int(sys.argv[2])
myLbb = int(sys.argv[3])
lp_density = float(sys.argv[4])
nskip = int(sys.argv[5])

myNpol = 1
myLpol = int(myLbb*(2.+(2.**(mygen+1.)-1.)))

os.chdir(myfold)
myfold=os.getcwd()
mypath = os.path.abspath(myfold)
filelist=list()

temppath = myfold #os.path.join(dirname, subdirname)
filelist.append(sorted(glob.glob(temppath+'/G{}*.trj'.format(gen))))


flat_list = np.asarray([item for sublist in filelist for item in sublist])

newlist = np.empty(len(flat_list), dtype = int)
k=0
for item in flat_list:
    elems = item.strip(".trj").split("_"); 
    #elems = item.strip(".lammpstrj").split(".")
    newlist[k] = int(elems[-1]); k+=1

idlist = newlist[newlist.argsort()]
Deltat = 0.01 ##float(elems[2])

times = []

polymer = mypolymer(myLbb)
Nbins = polymer.rhist.shape[0]

jj = 0; _z = 0
totrep = 0
for ff in flat_list[newlist.argsort()]:

    print(ff)

    if jj < nskip:
        jj += 1
        continue

    fz = open(ff,"r")
        
    _time, mybox, nparticles, data = read_lammpstrj(fz)	
    fz.close()

    if totrep == 0:
        time0 = _time

    mypolymers = data[:,2:5] + mybox*data[:,5:8]
    temp = data[np.where(data[:,1] == 1)]
    mybackbone = temp[:,2:5] + mybox*temp[:,5:8]
    temp = data[np.where(data[:,1] == 3)]       
    dendrimer = temp[:,2:5] + mybox*temp[:,5:8]
    times.extend([(_time-time0)*Deltat])

    #CHAINS
    myCoM = np.empty((myNpol,3),dtype = 'float')
        
    _z = idlist[jj];
    
    for i in range(myNpol):
        poltemp = mypolymers[i*myLpol:(i+1)*myLpol,:]; 
        CoM = center_of_mass(poltemp, mybox, 'nopbc'); 
       
        rhist, thhist = radial_density_cylinder(dendrimer, mybackbone, int(lp_density), mygen, Nbins)   

        polymer.rhist[:,0] = rhist[0];  polymer.rhist[:,1] += rhist[1]
        polymer.thhist[:,0] = thhist[0]; polymer.thhist[:,1] += thhist[1]

    totrep = totrep+1
    jj = jj+1

print("finished")

norm = np.trapz(polymer.rhist[:,1], polymer.rhist[:,0]); polymer.rhist[:,1] /= norm; 
norm = np.trapz(polymer.thhist[:,1], polymer.thhist[:,0]); polymer.thhist[:,1] /= norm;

np.savetxt(mypath+"/cylindrical_radial_density_%d.dat" % (int(lp_density)), polymer.rhist, fmt='%.8e') 
np.savetxt(mypath+"/cylindrical_angle_density_%d.dat" % (int(lp_density)), polymer.thhist, fmt='%.8e')

print("all done")

