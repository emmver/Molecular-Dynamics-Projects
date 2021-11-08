import hoomd
import hoomd.md
import numpy as np 
import sys
import os 
from utils_rebuild import read_lammpstrj, mybroadcast_to, C_pd_L, rotation_matrix, center_of_mass 
#from myanalysis import center_of_mass

hoomd.context.initialize("");

def make_dendron(gen):

    outmap = [[0,1],[0,2]]
    idx = [0, [1,2]]
    
    ww=3
    thetas = [20,20,20,20,20,20,20]
    
    for i in range (2,gen+1):
        blobs=2**(gen+1)-1
        dendron=np.zeros((blobs,3))
        a=np.zeros((1,3))
        theta0 = 60
        b=np.array([1*np.cos(theta0*np.pi/180),1*np.sin(theta0*np.pi/180),0])+np.array([0,0.1,0.1])
        c=np.array([1*np.cos((180-theta0)*np.pi/180),1*np.sin((180-theta0)*np.pi/180),0])+np.array([0,0.1,-0.1])
       	dendron_1 = [a, [b,c]]
        dendron_2 = [a, [b,c]]
        
        outmap = [[0,1],[0,2]]
        idx = [0, [1,2]]
    
        ww=3
        thetas = [20,20,20,20,20,20,20]
        for i in range (2,gen+1):
            dendron_1.append([])
            dendron_2.append([])
            idx.append([])
            curr_g = dendron_2[i-1]
            curridx = idx[i-1]

            rot=np.array([np.cos(thetas[i]*np.pi/180),np.sin(thetas[i]*np.pi/180),0,-np.sin(thetas[i]*np.pi/180),np.cos(thetas[i]*np.pi/180),0,0,0,1])
            rot=rot.reshape(3,3)
            rot1=np.copy(rot) 
            rot1[0,1] = -rot1[0,1]; rot1[1,0] = -rot1[1,0] #Changing first and second row?
        
            for j in range(2**(i-1)):
                v = dendron_1[i-1][j]
                ref_idx = curridx[j]
                for k in [0,1]:
                    if k == 0:
                        dendron_1[i].extend([np.dot(rot,v)])
                        dendron_2[i].extend([np.dot(rot,v)+curr_g[j]])
                        idx[i].extend([ww]); outmap.append([ref_idx,ww]); ww+=1
                    else:	
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



def assemble_unit(monomer, side_mon, dendron): #assembles monomer linker and dendron

    out = np.empty((dendron.shape[0]+1+1,3)) #dendron.shape[0] gives number of rows 

    out[0,:] = monomer
    out[1,:] = side_mon
    _dir = side_mon - monomer; _dir /= np.linalg.norm(_dir)
    CdM = center_of_mass(dendron[1:],(1,1,1), 'nopbc')
    _dir2 = CdM - dendron[0]; _dir2 /= np.linalg.norm(_dir2)
    nn, nn, mat = rotation_matrix(_dir2, _dir)
    oldd = np.copy(dendron);
    for i in range(oldd.shape[0]):
        dendron[i] = np.dot(mat,oldd[i])
    out[2:,:] = dendron + mybroadcast_to(monomer + 2*_dir, (out.shape[0]-2,3))

    return out


Nmon =int(sys.argv[1]) #t(sys.argv[1]) # of monomers
gen = int(sys.argv[2]) ##int(sys.argv[2])   # generation
xbox =1.5*Nmon##float(sys.argv[3])
ybox = xbox #float(sys.argv[4])
zbox =xbox #float(sys.argv[5])
path_for_traj= sys.argv[3]
newgen = gen
ldendron = 2**(gen+1)-1

allparticles = np.empty((Nmon*(2+(2**(gen+1)-1)),3), dtype = float)
dendron, bondmap= make_dendron(newgen)
path_for_store="G%d_%d"%(gen,Nmon)
didx=0

if os.path.isdir(path_for_store)==False:
	os.mkdir(path_for_store)
while os.path.isdir(path_for_store)==True:
	didx+=1
	path_for_store="G%d_%d_run%d"%(gen,Nmon,didx)

os.mkdir(path_for_store)
snapshot0 = hoomd.data.gsd_snapshot(filename=path_for_traj+"/traj.gsd",frame=-1) ##hoomd.init.read_gsd(filename="traj0.gsd",frame=-1)

snapshot = hoomd.data.make_snapshot(N=allparticles[:,0].size,box=hoomd.data.boxdim(Lx=xbox, Ly=ybox, Lz=zbox),particle_types=['A', 'B','C'],bond_types=['backbone','linker','dendron'],angle_types=['backbone','linker'])
### Set Positions ####
kk = 0
for i in range(Nmon):
    
    monomer = snapshot0.particles.position[i*(2+2**(1+1)-1)]; side_m = snapshot0.particles.position[i*(2+2**(1+1)-1)+1]
    dendronized = assemble_unit(monomer, side_m, dendron)
     
    #allparticles[i*(2+(2**(gen+1)-1)):(i+1)*(2+(2**(gen+1)-1)),:] =np.copy(dendronized)
    #dendronized=dendronized-np.array([xbox/2,xbox/2,xbox/2])
    
    snapshot.particles.position[kk]=[dendronized[0,0],dendronized[0,1],dendronized[0,2]];
    snapshot.particles.typeid[kk]=0;
    kk+=1
    snapshot.particles.position[kk]=[dendronized[1,0],dendronized[1,1],dendronized[1,2]];
    snapshot.particles.typeid[kk]=1;
    kk+=1
    for j in range((2**(gen+1)-1)):
        snapshot.particles.typeid[kk]=2;
        snapshot.particles.position[kk]=[dendronized[j+2,0],dendronized[j+2,1],dendronized[j+2,2]];kk+=1

### Randomize Velocities #####

#snapshot.thermalize_particle_momenta(filter=hoomd.filter.All, kT=1.0)
snapshot.bonds.resize(Nmon-1+2*Nmon+Nmon*len(bondmap))
kk=0
for k in range (0,Nmon-1):
    snapshot.bonds.group[kk]=[k*(2+(2**(gen+1)-1)),(k+1)*(2+(2**(gen+1)-1))]
    snapshot.bonds.typeid[kk]=0
    kk+=1
for k in range (Nmon):
    snapshot.bonds.group[kk]=[k*(2+(2**(gen+1)-1))+1-1, k*(2+(2**(gen+1)-1))+2-1]
    snapshot.bonds.typeid[kk]=1
    kk+=1
    snapshot.bonds.group[kk]=[k*(2+(2**(gen+1)-1))+2-1, k*(2+(2**(gen+1)-1))+3-1]
    snapshot.bonds.typeid[kk]=1
    kk+=1
for k in range (Nmon):
    for w,l in bondmap:
        snapshot.bonds.group[kk]=[k*(2+(2**(gen+1)-1))+3+w-1,k*(2+(2**(gen+1)-1))+3+l-1];
        snapshot.bonds.typeid[kk]=2;
        kk+=1

#### Angles #####
kk=0
snapshot.angles.resize(Nmon-2+Nmon)

for k in range (Nmon-2):
	snapshot.angles.group[kk]=[k*(2+(2**(gen+1)-1)),(k+1)*(2+(2**(gen+1)-1)),(k+2)*(2+(2**(gen+1)-1))];
	snapshot.angles.typeid[kk]=0;kk+=1
for k in range(Nmon):
	snapshot.angles.group[kk]=[k*(2+(2**(gen+1)-1)),k*(2+(2**(gen+1)-1))+1,k*(2+(2**(gen+1)-1))+2]
	snapshot.angles.typeid[kk]=1;kk+=1

hoomd.init.read_snapshot(snapshot)

nl=hoomd.md.nlist.tree()
### Pair Potential ##########
#wca=hoomd.md.pair.force_shifted_lj(r_cut=2**(1/6),nlist=nl)
#wca.pair_coeff.set('A','A',epsilon=1.0,sigma=1.0,r_cut=2**(1/6))
#wca.pair_coeff.set('B','B',epsilon=1.0,sigma=1.0,r_cut=2**(1/6))
#wca.pair_coeff.set('A','B',epsilon=1.0,sigma=1.0,r_cut=2**(1/6))
#wca.pair_coeff.set('A','C',epsilon=1.0,sigma=0.9,r_cut=0.9*2**(1/6))
#wca.pair_coeff.set('B','C',epsilon=1.0,sigma=0.9,r_cut=0.9*2**(1/6))
#wca.pair_coeff.set('C','C',epsilon=1.0,sigma=0.8,r_cut=0.8*2**(1/6))

gauss = hoomd.md.pair.gauss(r_cut=2**(1/6.), nlist = nl)
gauss.pair_coeff.set('A','A',epsilon=10.0,sigma=0.25,r_cut=2**(1/6))
gauss.pair_coeff.set('B','B',epsilon=10.0,sigma=0.25,r_cut=2**(1/6))
gauss.pair_coeff.set('A','B',epsilon=10.0,sigma=0.25,r_cut=2**(1/6))
gauss.pair_coeff.set('A','C',epsilon=10.0,sigma=0.18,r_cut=0.9*2**(1/6))
gauss.pair_coeff.set('B','C',epsilon=10.0,sigma=0.18,r_cut=0.9*2**(1/6))
gauss.pair_coeff.set('C','C',epsilon=10.0,sigma=0.16,r_cut=0.8*2**(1/6))

#morse = hoomd.md.pair.morse(r_cut=2**(1/6.), nlist = nl)
#morse.pair_coeff.set('A','A',D0=0.05,alpha=10.,r0=1.2,r_cut=2**(1/6))
#morse.pair_coeff.set('B','B',D0=0.05,alpha=10.,r0=1.2,r_cut=2**(1/6))
#morse.pair_coeff.set('A','B',D0=0.05,alpha=10.,r0=1.2,r_cut=2**(1/6))
#morse.pair_coeff.set('A','C',D0=0.05,alpha=10.,r0=1.08,r_cut=0.9*2**(1/6))
#morse.pair_coeff.set('B','C',D0=0.05,alpha=10.,r0=1.08,r_cut=0.9*2**(1/6))
#morse.pair_coeff.set('C','C',D0=0.05,alpha=10.,r0=0.97,r_cut=0.8*2**(1/6))

### Bonding Potential #######
#fene=hoomd.md.bond.fene()
#fene.bond_coeff.set('backbone',k=30.0,r0=1.5,sigma=1.0,epsilon=1.0)
#fene.bond_coeff.set('linker',k=30.0,r0=1.5,sigma=1.0,epsilon=1.0)
#fene.bond_coeff.set('dendron',k=30.0,r0=1.5,sigma=0.9,epsilon=1.0)

harm = hoomd.md.bond.harmonic()
harm.bond_coeff.set('backbone',k=100.0,r0=1.0)
harm.bond_coeff.set('linker',k=100.0,r0=1.0,)
harm.bond_coeff.set('dendron',k=100.0,r0=1.0)

###Bending Potential####
angular=hoomd.md.angle.harmonic()
angular.angle_coeff.set('backbone',k=6.0,t0=np.pi)
angular.angle_coeff.set('linker',k=240.0,t0=np.pi)
itg=hoomd.md.integrate.langevin(group=hoomd.group.all(),kT=1.0,seed=np.random.randint(1e6)+1)
#itg.randomize_velocities(kT=1.0,s)
itg.set_gamma('A',gamma=10.0)
itg.set_gamma('B',gamma=10.0)
itg.set_gamma('C',gamma=10.0)

hoomd.md.integrate.mode_standard(dt=0.0001);
to_run=1e8
frames=1000
hoomd.dump.gsd(filename=path_for_store+'/traj.gsd',period=int(to_run/frames), group=hoomd.group.all(),overwrite=True)
hoomd.analyze.log(filename=path_for_store+"/log-output.log",quantities=['potential_energy', 'temperature'],period=int(to_run/1000) ,overwrite=True)

hoomd.run (to_run)

