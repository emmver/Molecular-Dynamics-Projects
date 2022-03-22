import hoomd
import hoomd.md
import numpy as np 
import sys 
import os 
######################################################################################################################################################################################################################
def make_dendron(gen):

    blobs=2**(gen+1)-1
    dendron=np.zeros((blobs,3))
    a=np.zeros((1,3))
    theta0 = 60
    b=np.array([1*np.cos(theta0*np.pi/180),1*np.sin(theta0*np.pi/180),1.0])+np.array([0,0,1])
    c=np.array([1*np.cos((180-theta0)*np.pi/180),1*np.sin((180-theta0)*np.pi/180),-1.0])+np.array([0,0,-1])
    dendron_1 = [a, [b,c]]
    dendron_2 = [a, [b,c]]
    dendron_3 = [a, [b,c]]
    outmap = [[0,1],[0,2]]
    idx = [0, [1,2]]
    
    ww=3
    thetas = [20,20,20,20,20,20,20]
# xunit=np.array([1,0,0])
# zunit=np.array([0,0,1])
# yunit=np.array([0,1,0])
    for i in range (2,gen+1):
        blobs=2**(gen+1)-1
        dendron=np.zeros((blobs,3))
        a=np.zeros((1,3))
        theta0 = 60
        b=np.array([1*np.cos(theta0*np.pi/180),1*np.sin(theta0*np.pi/180),0])+np.array([0,0.1,0.1])
        c=np.array([1*np.cos((180-theta0)*np.pi/180),1*np.sin((180-theta0)*np.pi/180),0])+np.array([0,0.1,-0.1])
        dendron_1 = [a, [b,c]]
        dendron_2 = [a, [b,c]]
        #print ("b is =",b,"c is =",c)
        outmap = [[0,1],[0,2]]
        #unitvect=[[0],[1,2]]
        idx = [0, [1,2]]
    
        ww=3
        thetas = [20,20,20,20,20,20,20]
        for i in range (2,gen+1):
            dendron_1.append([])
            dendron_2.append([])
            #dendron_3.append([])
            idx.append([])
            #print("generation", i)
            curr_g = dendron_2[i-1]
            curridx = idx[i-1]

            rot=np.array([np.cos(thetas[i]*np.pi/180),np.sin(thetas[i]*np.pi/180),0,-np.sin(thetas[i]*np.pi/180),np.cos(thetas[i]*np.pi/180),0,0,0,1])
            rot=rot.reshape(3,3)
            rot1=np.copy(rot) 
            rot1[0,1] = -rot1[0,1]; rot1[1,0] = -rot1[1,0] #Changing first and second row?
#         rotx1=np.array([1.,0.,0.,0.,np.cos(phi*np.pi/180),-np.sin(phi*np.pi/180),0,np.sin(phi*np.pi/180),np.cos(phi*np.pi/180)])
#         rotx1=rotx1.reshape(3,3)
#         rotx2=np.copy(rotx1) 
#         rotx2[2,1] = -rotx1[2,1]; rotx2[1,2] = -rotx1[1,2]
        
            for j in range(2**(i-1)):
                    v = dendron_1[i-1][j]
                    ref_idx = curridx[j]
                    for k in [0,1]:
                        
                        if k == 0:
                            dendron_1[i].extend([np.dot(rot,v)])
                            #+np.array([0.5,0.2,0.04])
                            dendron_2[i].extend([np.dot(rot,v)+curr_g[j]])
                        #dendron_3[i].extend(np.dot(([np.dot(rot,v)+curr_g[j]]),xunit))
                            idx[i].extend([ww]); outmap.append([ref_idx,ww]); ww+=1
                        else:
                            #+np.array([-0.5,-0.2,-0.04])
                            dendron_1[i].extend([np.dot(rot1,v)])
                            dendron_2[i].extend([np.dot(rot1,v)+curr_g[j]])
                            idx[i].extend([ww]); outmap.append([ref_idx,ww]); ww+=1
                            
        kk=1
        for i in range(1,gen+1):
            kstart=kk
            for j in dendron_2[i]:
                dendron[kk,:] = j
                kk+=1

    return dendron, outmap




def make_small(filename,Nmon,gen,scale_factor):
	frame=np.random.randint(70)
	snapshot=hoomd.data.gsd_snapshot(filename=filename,frame=frame)
	#print(snapshot.box)
	ldendron = 2**(gen+1)-1
	npart=Nmon* (2**(gen+1)+1)
	stop=int(snapshot.particles.N/scale_factor)
	#stop=int(snapshot.particles.N)
	system=hoomd.data.make_snapshot(N=npart,box=snapshot.box,particle_types=['A','B','C'],bond_types=['backbone','linker','dendron'],angle_types=['backbone','linker'])
	print(system.particles.position[:,0].size)
	for i in range(stop):
		#system.particles.add(i)
		system.particles.position[i]=snapshot.particles.position[i]
		system.particles.orientation[i]=snapshot.particles.orientation[i]
		system.particles.velocity[i]=snapshot.particles.velocity[i]

	kk=0
	dendron, bondmap= make_dendron(gen)
	system.bonds.resize(Nmon-1+2*Nmon+Nmon*len(bondmap))
	system.angles.resize(Nmon-2+Nmon)

	for k in range (0,Nmon-1):
		system.bonds.group[kk]=[k*(2+(2**(gen+1)-1)),(k+1)*(2+(2**(gen+1)-1))]
		#print("Backbone:",(k)*(2+(2**(gen+1)-1)),(k+1)*(2+(2**(gen+1)-1)))
		system.bonds.typeid[kk]=0
		kk+=1
	for k in range (Nmon):
		system.bonds.group[kk]=[k*(2+(2**(gen+1)-1))+1-1, k*(2+(2**(gen+1)-1))+2-1]
		#print("Linker 1")
		#print(k*(2+(2**(gen+1)-1))+1-1, k*(2+(2**(gen+1)-1))+2-1)
		snapshot.bonds.typeid[kk]=1
		kk+=1
		system.bonds.group[kk]=[k*(2+(2**(gen+1)-1))+2-1, k*(2+(2**(gen+1)-1))+3-1]
	# print("Linker 2")
		#print(k*(2+(2**(gen+1)-1))+2-1, k*(2+(2**(gen+1)-1))+3-1)
		system.bonds.typeid[kk]=1
		kk+=1
	for k in range (Nmon):
		for w,l in bondmap:
			system.bonds.group[kk]=[k*(2+(2**(gen+1)-1))+3+w-1,k*(2+(2**(gen+1)-1))+3+l-1];
	#       print("Dendron")
	#      print(k*(2+(2**(gen+1)-1))+3+w-1)
	#     print(k*(2+(2**(gen+1)-1))+3+l-1)
			system.bonds.typeid[kk]=2;
			kk+=1
	#print("Bonds Added:",kk)

	#### Angles #####
	kk=0

	for k in range (Nmon-2):
		system.angles.group[kk]=[k*(2+(2**(gen+1)-1)),(k+1)*(2+(2**(gen+1)-1)),(k+2)*(2+(2**(gen+1)-1))];
		system.angles.typeid[kk]=0;kk+=1
	for k in range(Nmon):
		system.angles.group[kk]=[k*(2+(2**(gen+1)-1)),k*(2+(2**(gen+1)-1))+1,k*(2+(2**(gen+1)-1))+2]
		system.angles.typeid[kk]=1;kk+=1

	return system
################################################################################################################################################################


hoomd.context.initialize("");
Nmon =int(sys.argv[1]) # of monomers
gen = int(sys.argv[2])   # generation
scale_factor=2
xbox =1.5*Nmon##float(sys.argv[3])
ybox = xbox #float(sys.argv[4])
zbox =xbox #float(sys.argv[5])
filename=sys.argv[3]
system=make_small(filename,Nmon,gen,scale_factor)
hoomd.init.read_snapshot(system)
print(hoomd.group.all())
#path_to_store='G%d_%d'%(gen,Nmon)
path_to_store=sys.argv[4]
# top_path=os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(sys.argv[3]))))#os.path.dirname(os.path.dirname(sys.argv[3]))))
# os.chdir(top_path)
# if os.path.isdir(path_to_store)==False:
# 	os.mkdir(path_to_store)
#print('Path',curr_path)
#sys.exit()
#print(snapshot.particles.id)
print('################################################################')
print('########################## F I N A L #############################')

print('Particles after cut',system.particles.N)
print('Bonds after cut',system.bonds.N)
print('Angles after cut',system.angles.N)


print('################################################################')



#hoomd.dump.gsd(filename='cut_half.gsd', period=None,group=hoomd.group.all(),overwrite=False)

#sys.exit()
nl=hoomd.md.nlist.tree()
### Pair Potential ##########
wca=hoomd.md.pair.force_shifted_lj(r_cut=2**(1/6),nlist=nl)
wca.pair_coeff.set('A','A',epsilon=1.0,sigma=1.0,r_cut=2**(1/6))
wca.pair_coeff.set('B','B',epsilon=1.0,sigma=1.0,r_cut=2**(1/6))
wca.pair_coeff.set('A','B',epsilon=1.0,sigma=1.0,r_cut=2**(1/6))
wca.pair_coeff.set('A','C',epsilon=1.0,sigma=0.9,r_cut=0.9*2**(1/6))
wca.pair_coeff.set('B','C',epsilon=1.0,sigma=0.9,r_cut=0.9*2**(1/6))
wca.pair_coeff.set('C','C',epsilon=1.0,sigma=0.8,r_cut=0.8*2**(1/6))

#gauss = hoomd.md.pair.gauss(r_cut=2**(1/6.), nlist = nl)
#gauss.pair_coeff.set('A','A',epsilon=200.0,sigma=0.25,r_cut=2**(1/6))
#gauss.pair_coeff.set('B','B',epsilon=200.0,sigma=0.25,r_cut=2**(1/6))
#gauss.pair_coeff.set('A','B',epsilon=200.0,sigma=0.25,r_cut=2**(1/6))
#gauss.pair_coeff.set('A','C',epsilon=200.0,sigma=0.18,r_cut=0.9*2**(1/6))
#gauss.pair_coeff.set('B','C',epsilon=200.0,sigma=0.18,r_cut=0.9*2**(1/6))
#gauss.pair_coeff.set('C','C',epsilon=200.0,sigma=0.16,r_cut=0.8*2**(1/6))

#morse = hoomd.md.pair.morse(r_cut=2**(1/6.), nlist = nl)
#morse.pair_coeff.set('A','A',D0=1.,alpha=10.,r0=1.2,r_cut=2**(1/6))
#morse.pair_coeff.set('B','B',D0=1.,alpha=10.,r0=1.2,r_cut=2**(1/6))
#morse.pair_coeff.set('A','B',D0=1.,alpha=10.,r0=1.2,r_cut=2**(1/6))
#morse.pair_coeff.set('A','C',D0=1.,alpha=10.,r0=1.08,r_cut=0.9*2**(1/6))
#morse.pair_coeff.set('B','C',D0=1.,alpha=10.,r0=1.08,r_cut=0.9*2**(1/6))
#morse.pair_coeff.set('C','C',D0=1.,alpha=10.,r0=0.97,r_cut=0.8*2**(1/6))

### Bonding Potential #######
fene=hoomd.md.bond.fene()
fene.bond_coeff.set('backbone',k=30.0,r0=1.5,sigma=1.0,epsilon=1.0)
fene.bond_coeff.set('linker',k=30.0,r0=1.5,sigma=1.0,epsilon=1.0)
fene.bond_coeff.set('dendron',k=30.0,r0=1.5,sigma=0.9,epsilon=1.0)

#harm = hoomd.md.bond.harmonic()
#harm.bond_coeff.set('backbone',k=100.0,r0=1.0)
#harm.bond_coeff.set('linker',k=100.0,r0=1.0,)
#harm.bond_coeff.set('dendron',k=100.0,r0=1.0)

###Bending Potential####
angular=hoomd.md.angle.harmonic()
angular.angle_coeff.set('backbone',k=6.0,t0=np.pi)
angular.angle_coeff.set('linker',k=240.0,t0=np.pi)
itg=hoomd.md.integrate.langevin(group=hoomd.group.all(),kT=1.0,seed=np.random.randint(1e6)+1)
#itg.randomize_velocities(kT=1.0,s)
itg.set_gamma('A',gamma=1.0)
itg.set_gamma('B',gamma=1.0)
itg.set_gamma('C',gamma=1.0)
traj_fname='traj.gsd'
tidx=0
while os.path.isfile(path_to_store+'/'+traj_fname)==True:
	tidx+=1
	traj_fname='traj_%d.gsd'%(tidx)
hoomd.md.integrate.mode_standard(dt=0.005);
to_run=1e9
frames=1000
hoomd.dump.gsd(filename=path_to_store+'/'+traj_fname,period=int(to_run/frames), group=hoomd.group.all(),overwrite=False)
hoomd.analyze.log(filename=path_to_store+'/'+"log-output_%d.log"%(tidx),quantities=['potential_energy', 'temperature'],period=int(to_run/1000) ,overwrite=True)

hoomd.run (to_run)

