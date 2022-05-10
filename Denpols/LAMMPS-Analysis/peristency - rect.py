import numpy as np 
from matplotlib import pyplot as plt 
##from numba import njit
from mpl_toolkits.mplot3d import Axes3D
from utilities import read_lammpstrj
import sys
import os.path



_file=sys.argv[1]
f=open (_file , "r")
path=os.path.realpath(f.name)
rpath=os.path.dirname(path)
Nbb=int(input("Backbone monomers?"))
gen=int(input("Generation?"))
frames=int(input("Frames?"))
filename=input("Filename")
lp_proj=[]
bonds=[]

for i in range (0,Nbb-1):
	bonds.extend([(i+1)/Nbb])

for i in range (frames):
    print (i)
    out,time,mybox,nparticles = read_lammpstrj(f)
    unwrap=np.zeros(np.shape(out[:,2:5]))
    print ('Got trajectory',mybox)
    #data1 = data[:,2:5] + box*data[:,5:]
    unwrap[:,0] = out[:,2]+out[:,5]*mybox/1.85714285714
    unwrap[:,1] = out[:,3]+out[:,6]*mybox/1.85714285714
    unwrap[:,2] = out[:,4]+out[:,7]*mybox


    endtoend=unwrap[(Nbb-1)*(1+2**(gen+1))]-unwrap[0];
    tmp=[]

    for j in range (Nbb-1):
       # print(j*(1+2**(gen+1)))
        bondvect=unwrap[(j+1)*(1+2**(gen+1))]-unwrap[j*(1+2**(gen+1))]
        bondlength=np.linalg.norm(bondvect)
        calc=np.dot(bondvect,endtoend)/bondlength**2
        tmp.extend([calc])

    lp_proj.append(tmp)
    #plt.plot(bonds,lp_proj[i])




       
   # jj +=1
    #print (jj)
   

    
f.close

plt.show()
kk=0
lp_proj_f=lp_proj[0]
#print(len(lp_proj[0]))
for i in range (1,frames-1):
    lp_proj_f[:]=[sum(x) for x in zip(lp_proj_f,lp_proj[i])]
    kk+=1
    #print ("Length",len(lp_proj_f))
lp_proj_f[:]=[x/kk for x in lp_proj_f]
plt.plot(bonds,lp_proj_f,"ko-")
plt.show()

outfiles=os.path.join(rpath,filename+".txt")

g=open(outfiles,"w")
for i in range (len(bonds)):
    g.write("%.5f\t%.5f\n" % (bonds[i],lp_proj_f[i]))
g.close()

