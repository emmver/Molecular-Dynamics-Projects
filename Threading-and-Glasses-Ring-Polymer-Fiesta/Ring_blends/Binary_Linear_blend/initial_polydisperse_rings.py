import numpy as np
import math
from scipy.stats import special_ortho_group
import os
import time

np.random.seed(int(time.time()))#22378)
for nca in [10]:#[10,30,50,70,90,100,110,130,150,170,190]:
	for dirnr in range(1):
		mpclist = [200,200]#[400,800]#
		outpath = './initial_data/initialN'+str(mpclist[0])+'N'+str(mpclist[1])+'nca'+str(nca)+'data/dir'+str(dirnr)+'/'
		os.system('mkdir -p '+outpath)
		maxmpc = max(mpclist)
		nc = 10

		#lattice constant
		la = 3+maxmpc/3.14

		n_rows = int(math.ceil(np.power(nc,1.0/3.0)))
		print(n_rows)
		L = n_rows*la
		print(L, L , L) 

		Nlist = np.array([mpclist[0]]*int(nca) + [mpclist[1]]*int(nc-nca))
		#Nlist = np.random.choice(Nlist,len(Nlist),replace=False) #random order of big and small rings
		print(Nlist) 
		np.savetxt(outpath+'chain_lengths.dat',Nlist,fmt='%d')
		totnp = sum(Nlist)
		print(totnp)

		ring_primitive = []
		for mpc in mpclist:
			r = mpc/np.pi/2.0
			ring_coords = []
			for i in range(mpc):
				phi = 2.0*np.pi*i/mpc
				x = r*np.cos(phi)
				y = r*np.sin(phi)
				z = 0.0
				ring_coords.append([x,y,z])
			ring_primitive.append(ring_coords)

		print(len(ring_primitive[0]))
		print(len(ring_primitive[1]))

		data = []
		data2 = []
		bonds = []
		angles=[]
		ch=0
		pid=0

		alist = np.random.choice(range(n_rows*n_rows*n_rows),n_rows*n_rows*n_rows,replace=False)
		chck = []
		print(alist)
		for u in alist:
			i=u%n_rows #int(math.floor(u/(n_rows*n_rows)))
			j=((u-i)/n_rows)%n_rows #int(math.floor((u-i*(n_rows*n_rows))/n_rows))
			k=(u-j*n_rows-i)/(n_rows*n_rows) #u%n_rows
			if [i,j,k] in chck:
				print("UAAAAAAAAA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
			else:
				chck.append([i,j,k])
			if ch>=nc:
				break
			R = np.array([-L/2.0 + la/2.0 + i*la,-L/2.0 + la/2.0 + j*la, -L/2.0 + la/2.0 + k*la	])				
			M = special_ortho_group.rvs(3) #random rotational matrix
			rp = [elem for elem in ring_primitive if len(elem)==Nlist[ch]][0]
			for ii in range(Nlist[ch]):
				rvec = R + np.dot(M,rp[ii])
				pid+=1

				data.append([pid,ch+1, 1, rvec[0],rvec[1],rvec[2]])
				data2.append([pid, 1, rvec[0],rvec[1],rvec[2]])

				if ii==0:
					bonds.append([pid, 1, pid+Nlist[ch]-1,pid])
				else:
					bonds.append([pid, 1, pid-1,pid])

				if ii==0:
					angles.append([pid, 1, pid+Nlist[ch]-2,pid+Nlist[ch]-1,pid])
				elif ii==1:
					angles.append([pid, 1, pid+Nlist[ch]-2,pid-1,pid])
				else:
					angles.append([pid, 1, pid-2,pid-1,pid])

			ch+=1

		with open(outpath+'head','w') as of:
			of.write('LAMMPS\n')
			of.write('\n')
			of.write(str(totnp)+' atoms\n')
			of.write(str(totnp)+' bonds\n')
			of.write(str(totnp)+' angles\n')
			of.write('\n')
			of.write('1 atom types\n')
			of.write('1 bond types\n')
			of.write('1 angle types\n')
			of.write('\n')
			of.write(str(-L/2.0)+' '+str(L/2.0)+' xlo xhi\n')
			of.write(str(-L/2.0)+' '+str(L/2.0)+' ylo yhi\n')
			of.write(str(-L/2.0)+' '+str(L/2.0)+' zlo zhi\n')
			of.write('\n')
			of.write('Atoms\n')
			of.write('\n')

		with open(outpath+'head2','w') as of:
			of.write('ITEM: TIMESTEP\n')
			of.write('0\n')
			of.write('ITEM: NUMBER OF ATOMS\n')
			of.write(str(totnp)+'\n')
			of.write('ITEM: BOX BOUNDS pp pp pp\n')
			of.write(str(-L/2.0)+' '+str(L/2.0)+'\n')
			of.write(str(-L/2.0)+' '+str(L/2.0)+'\n')
			of.write(str(-L/2.0)+' '+str(L/2.0)+'\n')
			of.write('ITEM: ATOMS id type xu yu zu\n')
				
		np.savetxt(outpath+'coords.dat',data,fmt = '%d %d %d %f.5 %f.5 %f.5')
		np.savetxt(outpath+'coords2.dat',data2,fmt = '%d %d %f.5 %f.5 %f.5')
		np.savetxt(outpath+'bonds.dat',bonds,fmt = '%d %d %d %d')				
		np.savetxt(outpath+'angles.dat',angles,fmt = '%d %d %d %d %d')				

		os.system('cat '+outpath+'head > '+outpath+'pos.lammps')
		os.system('cat '+outpath+'coords.dat >> '+outpath+'pos.lammps')
		os.system('echo "\n"Bonds"\n" >> '+outpath+'pos.lammps')
		os.system('cat '+outpath+'bonds.dat >> '+outpath+'pos.lammps')
		os.system('echo "\n"Angles"\n" >> '+outpath+'pos.lammps')
		os.system('cat '+outpath+'angles.dat >> '+outpath+'pos.lammps')

		os.system('cat '+outpath+'head2 > '+outpath+'pos2.lammpstrj')
		os.system('cat '+outpath+'coords2.dat >> '+outpath+'pos2.lammpstrj')









