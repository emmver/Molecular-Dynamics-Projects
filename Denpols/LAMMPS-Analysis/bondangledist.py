import numpy as np
import sys
import os 
from library import kde_scipy,read_lammpstrj
#from utilities import 
from matplotlib import pyplot as plt

####################### I n p u t P a r a m e t e r s ###############################
myfile = sys.argv[1]
f = open(myfile,"r")
Nbb=int(sys.argv[2])


x_grid = np.linspace(-1, 1, 1000)
skipframes=int(sys.argv[3])
#frames=int(sys.argv[4])
gen=int(sys.argv[4])

path=os.path.realpath(__file__)
rpath=os.path.dirname(path)
rpath=os.path.join(rpath,"Bond_Dist/G{}".format(gen))
####################### S k i p   F r a m e s ###############################
for i in range (skipframes):
	time, box, nparticles, data = read_lammpstrj(f)
	print("Skipping...",i)



####################### L o o p   t h r o u g h   F r a m e s ###############################
cosavg=np.zeros((Nbb-2))
count=0
while 1:
	try:
		print("yay")
		j=1
		costheta=[]
		time, box, nparticles, data=read_lammpstrj(f)
		unwrap = data[:,2:5]+data[:,5:]*box

		while j<Nbb-1:
			bondvect1=unwrap[(j)*(1+2**(gen+1))]-unwrap[(j-1)*(1+2**(gen+1))]
			#bondvect1=bondvect1/np.linalg.norm(bondvect1)
			bondvect2=unwrap[(j+1)*(1+2**(gen+1))]-unwrap[j*(1+2**(gen+1))]
			#bondvect2=bondvect2/np.linalg.norm(bondvect2)
			dot_product=np.dot(bondvect1,bondvect2)
			#cosinus=np.arccos(dot_product)*(180/np.pi)
			magn=np.linalg.norm(bondvect1)*np.linalg.norm(bondvect2)
			cosinus=dot_product/magn

			costheta.append(cosinus)
			j+=1
			print("Done",j)
		
		count+=1
		costheta=np.array(costheta)
		#cosavg=np.array(cosavg)
		cosavg=cosavg+costheta
		#cosavg=[sum(x) for x in zip(cosavg, costheta)];print(cosavg)
		#print(cosavg)
		print("Calculating..",count)
	except Exception as e:
		print(e)
		break
myfile='avg_bond_dist_g{}{}'.format(int(gen),int(Nbb))
cosavg=np.asarray(cosavg)
cosavg=cosavg/count;print(cosavg)
pdf=kde_scipy(cosavg,x_grid,0.01)
plt.plot(x_grid,pdf)
#plt.yscale('log')
#plt.ylim(1e-4,100)
plt.show()
#plt.hist(cosavg,bins=5)
#plt.show()
outfiles=os.path.join(rpath,myfile+".txt")
np.savetxt(outfiles,np.stack((x_grid,pdf),axis=1))