import numpy as np
import sys
#import string
import ctypes
import os
import glob

#sys.path.insert(1,"./ANALYSIS/")

from utilities import print_data, random_on_sphere
from mykgrid import *
from utilities import read_lammpstrj, center_of_mass, avg_bending_angle


_sk = ctypes.CDLL('./libsk.so')
_sk.calc_sk.argtypes = (ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_int, ctypes.c_int)
#_sk.calc_sk_2.argtypes = (ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_int, ctypes.c_int)

_mypath = sys.argv[1]
##_file = _mypath+'/trj2.lammpstrj' ##sys.argv[1]
mygen = int(sys.argv[2])
Lbackbone = int(sys.argv[3])
mybox = float(sys.argv[4])
nskip = int(sys.argv[5])
print("Run in Folder",_mypath)
files=_mypath + '/G{}{}_*.lammpstrj'.format(mygen,Lbackbone)
#print("Files:",files)
filelist=sorted(glob.glob(files))
#filelist.extend((glob.glob(_mypath+'/G{}{}_*.lammpstrj'.format(mygen,Lbackbone))))
#print(filelist)

#flat_list = np.asarray([item for item in filelist])
flat_list=filelist

newlist = np.empty(len(flat_list), dtype = np.int_)
k=0
for item in flat_list:
        elems = str.split(item.strip(".lammpstrj"), "_");
        newlist[k] = int(elems[-1]); k+=1

newlist = np.asarray(newlist)
#print(newlist)

nk = int(np.log(10./0.001)/np.log(1.05)); print ("How many k",nk)
kgrid, knorm, nkbins= generate_kgrid_rand(nk, 0.001, 1.05)

Sskarr = np.zeros((nkbins),dtype = np.double)
Ssktrans = np.zeros((nkbins),dtype = np.double)

Sskbb = np.zeros((nkbins),dtype = np.double)

dr = 0.1; Lmax = 200.; nbins=int(Lmax/dr)
c_r = np.zeros((nbins), dtype = float)
cr_bb = np.zeros((nbins), dtype = float)

Nth= 100
btheta = np.zeros((Nth))
#bcount = np.zeros((Lbackbone-1))

j = 0
k = 0

##frames=1000 #684  ##runs/framestep
jj = -1
totrep = 0
#newlist.argsort()
for ff in flat_list[:]:

    print ("Print FF:",ff)


    _f = open(ff,'r')
    time, mybox, nparticles, data = read_lammpstrj(_f)
    _f.close()
    #nparticles=Lbackbone*2**(gen+1)+1
    unwrap = data[:,2:5] #+ mybox*data[:,-3:]
    CoM = center_of_mass(unwrap, mybox, 'nopbc')

    unwrap2 = np.copy(unwrap); unwrap2[:] -= CoM
    
    _sk.calc_sk (ctypes.c_void_p(kgrid.ctypes.data), ctypes.c_void_p(unwrap.ctypes.data), ctypes.c_void_p(Sskarr.ctypes.data), ctypes.c_int(nparticles), ctypes.c_int(nkbins))

    ##calc_molecule_density(unwrap, CoM, c_r, dr, mybox, nbins)

    temp = data[np.where(data[:,1] == 1)]
    backbone = temp[:,2:5]# + mybox*temp[:,5:8]
    ##backbone = unwrap2[::(2+(2**(mygen+1)-1))] 

    _sk.calc_sk (ctypes.c_void_p(kgrid.ctypes.data), ctypes.c_void_p(backbone.ctypes.data), ctypes.c_void_p(Sskbb.ctypes.data), ctypes.c_int(Lbackbone), ctypes.c_int(nkbins))

    ##calc_molecule_density(backbone, CoM, cr_bb, dr, mybox, nbins)

    bonds = backbone[1:] - backbone[:-1]    
    blength = np.sum(bonds[:]**2, axis=1)    

    avg_bending_angle(bonds, blength, btheta)
    

    jj += 1
    #print("ADD")
    totrep += 1
    
#print("Totrep: ",totrep,"\n\n")
for i in range(nkbins):
    Ssktrans[i] = (1./(1.*totrep)) * Sskarr[i] / float(nparticles)**2
    
    
_mypath="./ANALYSIS/G{}".format(mygen)
print(_mypath)
toprint = normsk_rand(knorm, Ssktrans)
filename=_mypath+"/G{}{}".format(mygen,Lbackbone)
print("To Store:",filename)
np.savetxt(filename, toprint, fmt = '%.8e')

for i in range(nkbins):
    Ssktrans[i] = (1./(1.*totrep)) * Sskbb[i] / float(Lbackbone)**2
    
    
    
toprint = normsk_rand(knorm, Ssktrans)
filename=_mypath+"/G{}{}_lin".format(mygen,Lbackbone)
print("To Store:",filename)
np.savetxt(filename, toprint, fmt = '%.8e')


##############################

toprint = np.empty((Nth,2))

temp = np.linspace(-1,1,Nth+1)
toprint[:,0] = 0.5*(temp[1:]+temp[:-1])
norm = np.trapz(btheta, toprint[:,0])
toprint[:,1] = btheta/norm

np.savetxt(_mypath+"/G{}{}_distro_theta.dat".format(mygen,Lbackbone), toprint, fmt="%.5e")

'''
########################
toprint = np.empty((nbins,2))

toprint[:,0] = np.linspace(dr/2., Lmax, nbins, endpoint=False)
toprint[:,1] = c_r/(4.*np.pi*toprint[:,0]**2)

np.savetxt(_mypath+"/c_r.dat", toprint, fmt="%.5e")

dq = 0.001; qmax = 10; qbins = int(np.log(qmax/dq)/np.log(1.05))

sq,qarr =  calc_fouriergr(c_r, qbins, dq, nbins, dr, Lmax)

toprint =  np.empty((qbins,2))
toprint[:,0] = qarr
toprint[:,1] = sq

np.savetxt(_mypath+"/s_q.dat", toprint, fmt="%.5e")
######################################
toprint = np.empty((nbins,2))

toprint[:,0] = np.linspace(dr/2., Lmax, nbins, endpoint=False)
toprint[:,1] = cr_bb/(4.*np.pi*toprint[:,0]**2)

np.savetxt(_mypath+"/cr_bb.dat", toprint, fmt="%.5e")

dq = 0.001; qmax = 10; qbins = int(np.log(qmax/dq)/np.log(1.05))

sq,qarr =  calc_fouriergr(cr_bb, qbins, dq, nbins, dr, Lmax)

toprint =  np.empty((qbins,2))
toprint[:,0] = qarr
toprint[:,1] = sq

np.savetxt(_mypath+"/sq_bb.dat", toprint, fmt="%.5e")
'''

