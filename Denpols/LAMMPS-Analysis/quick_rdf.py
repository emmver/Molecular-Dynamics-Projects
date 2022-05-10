import numpy as np
import sys
import os 

from library import print_data,center_of_mass,calc_rdf,calc_sq,gyration
from library import read_lammpstrj

myfile = sys.argv[1]
myL = float(sys.argv[2])
skip = int(sys.argv[3])

dr = myL/(2.*nbins)
nbins = int(myL/(2.*dr))

rdf = np.zeros((nbins),dtype = float)

f = open(myfile,"r")

path=os.path.realpath(f.name)
rpath=os.path.dirname(path)
#filename=input("Filename")

jj = 0
Rg = 0
#runs=1E9
framestep=float(sys.argv[4])
frames=int(sys.argv[5])
#time=0
time, box, nparticles1, data = read_lammpstrj(f)
t0=time+1
f.close()
f = open(myfile,"r")

for jj in range (0,frames):
    
    #print (jj,(time*1/framestep))
    print("yay")
    time, box, nparticles, data = read_lammpstrj(f)
    if jj > int((time/t0)*framestep):
        print ("Exit",jj,framestep*time/t0)
        break    
    
    #print (time)

    data1 = data[::skip,2:5] + box*data[::skip,5:]
    
    CoM = center_of_mass(data1, box, 'nopbc');
    Rg += gyration(data1, CoM, box)
    
    calc_rdf(data1, rdf, dr, myL, nbins, 1.)
    jj +=1
    print (jj)


f.close()

print (Rg/float(jj))
module_path = os.path.dirname(os.path.realpath(myfile));print(module_path)
file=sys.argv[6]
#file = os.path.join(module_path, file)

toprint = print_data(rdf, [nbins, dr, 4./3.*np.pi*(nparticles1/skip), (jj), file], 'rdf')

#toprint = np.empty((nbins,2))
#toprint[:,0] = np.linspace(dr/2.,nbins+dr/2.,nbins)
#toprint[:,1] = rdf/float(jj)

#np.savetxt('rdf.dat', toprint, fmt='%.8e')

#dq = (50.)/2000
sq,qmin,qmax = calc_sq(toprint[:,1], 2000, dr, nbins,myL)
module_path = os.path.dirname(os.path.realpath(myfile));print(module_path)
file=sys.argv[7]
#file = os.path.join(module_path, file)
rho = 1 ##(nparticles/skip)/box[0]**3
print_data(sq, [qmin, qmax, 2000,  4.*np.pi*dr*rho, rho, myL/2., file], 'formfactor')

#outfiles=os.path.join(rpath,filename+".txt")
#np.savetxt(outfiles, sq, fmt='%.5e',delimiter=' ')  

