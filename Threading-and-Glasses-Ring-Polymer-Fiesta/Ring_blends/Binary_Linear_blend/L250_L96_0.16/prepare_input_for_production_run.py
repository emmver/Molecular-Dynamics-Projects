import numpy as np
import os

nmin=96#400
nmax=250#800
nca = 88#int(400000/nmin)
nc = 88+6#nca+1000

totnp = nmin*nca + nmax*(nc-nca)
L = np.power(totnp/0.85,1./3.)
if nmax==250:
	stepmax=1100000
if nmax==400:
	stepmax=100000
if nmax==800:
	stepmax=200000
elif nmax==1600:
	stepmax=400000
tp = 'N'+str(nmin)+'N'+str(nmax)+'nca'+str(nca)+'data'


for dirnr in range(1):
	path = './data'+'/'
	print (dirnr), 
	outpath = path
	#os.system('mkdir -p '+outpath)

	chlens = [nmin]*(nca)+[nmax]*(nc-nca)# map(int,np.loadtxt(path+'dir'+str(dirnr)+'/chain_lengths.dat'))
	data = np.loadtxt(path+'test_pos'+str(stepmax)+'.dat',skiprows=9)
	data = sorted(data,key= lambda part: part[0] )
	datalow = []
	datahigh = []
	chlsum = 0
	pidlow = 1
	pidhigh = 1+nmin*(nca)
	for chl in chlens:
		for par in data[chlsum:chlsum+chl]:
			if chl==nmin:
				par[0] = pidlow
				datalow.append(par)
				pidlow+=1
			elif chl==nmax:
				par[0] = pidhigh
				datahigh.append(par)
				pidhigh+=1
		chlsum+=chl

	np.savetxt('./dl.dat',datalow,fmt='%d %d %f %f %f')
	#np.savetxt('./dh.dat',datahigh,fmt='%d %d %f %f %f')

	fl = './final.dat'
	os.system('cat ./dl.dat > '+fl)
	os.system('cat ./dh.dat >> '+fl)

	os.system('rm dl.dat dh.dat')


	newdata = np.loadtxt(fl)
	Nlist = [nmin]*(nca)+[nmax]*(nc-nca)

	data = []
	data2=[]
	bonds = []
	angles = []
	chlensum = 0
	for ch in range(nc):
			for ii in range(Nlist[ch]):
				pid = newdata[chlensum+ii,0]
				rvec = newdata[chlensum+ii,2:5]

				data.append([pid,ch+1, 1, rvec[0],rvec[1],rvec[2]])
				data2.append([pid, 1, rvec[0],rvec[1],rvec[2]])

				if ii==0:
					#bonds.append([pid, 1, pid+Nlist[ch]-1,pid])
					continue
				else:
					bonds.append([pid, 1, pid-1,pid])

				if ii==0:
					#angles.append([pid, 1, pid+Nlist[ch]-2,pid+Nlist[ch]-1,pid])
					continue
				elif ii==1:
					#angles.append([pid, 1, pid+Nlist[ch]-2,pid-1,pid])
					continue
				else:
					angles.append([pid, 1, pid-2,pid-1,pid])

			chlensum+=Nlist[ch]


	with open(outpath+'head','w') as of:
		of.write('LAMMPS\n')
		of.write('\n')
		of.write(str(totnp)+' atoms\n')
		of.write(str(nc*(nmax-1))+' bonds\n')
		of.write(str(nc*(nmax-2))+' angles\n')
		of.write('\n')
		of.write('1 atom types\n')
		of.write('1 bond types\n')
		of.write('1 angle types\n')
		of.write('\n')
		of.write(str(-L/2.)+' '+str(L/2.)+' xlo xhi\n')
		of.write(str(-L/2.)+' '+str(L/2.)+' ylo yhi\n')
		of.write(str(-L/2.)+' '+str(L/2.)+' zlo zhi\n')
		of.write('\n')
		of.write('Atoms\n')
		of.write('\n')

	with open(outpath+'head2','w') as of:
		of.write('ITEM: TIMESTEP\n')
		of.write('0\n')
		of.write('ITEM: NUMBER OF ATOMS\n')
		of.write(str(totnp)+'\n')
		of.write('ITEM: BOX BOUNDS pp pp pp\n')
		of.write(str(-L/2.)+' '+str(L/2.)+'\n')
		of.write(str(-L/2.)+' '+str(L/2.)+'\n')
		of.write(str(-L/2.)+' '+str(L/2.)+'\n')
		of.write('ITEM: ATOMS id type xu yu zu\n')
				
	np.savetxt(outpath+'coords.dat',data,fmt = '%d %d %d %f %f %f')
	np.savetxt(outpath+'coords2.dat',data2,fmt = '%d %d %f %f %f')
	np.savetxt(outpath+'bonds.dat',bonds,fmt = '%d %d %d %d')				
	np.savetxt(outpath+'angles.dat',angles,fmt = '%d %d %d %d %d')				
	np.savetxt(outpath+'chain_lengths.dat',Nlist,fmt = '%d')				


	os.system('cat '+outpath+'head > '+outpath+'indata.lammps')
	os.system('cat '+outpath+'coords.dat >> '+outpath+'indata.lammps')
	os.system('echo "\n"Bonds"\n" >> '+outpath+'indata.lammps')
	os.system('cat '+outpath+'bonds.dat >> '+outpath+'indata.lammps')
	os.system('echo "\n"Angles"\n" >> '+outpath+'indata.lammps')
	os.system('cat '+outpath+'angles.dat >> '+outpath+'indata.lammps')

	os.system('cat '+outpath+'head2 > '+outpath+'pos.lammpstrj')
	os.system('cat '+outpath+'coords2.dat >> '+outpath+'pos.lammpstrj')

	os.system('rm final.dat')


