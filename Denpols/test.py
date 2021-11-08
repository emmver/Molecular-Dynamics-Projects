import hoomd
import hoomd.md
import numpy as np 
import os 
import sys
hoomd.context.initialize("");

def print_lammpstrj(myinput, thisL, myname):

        f = open(myname,"w")
        howmany = myinput.shape[0]

        f.write("ITEM: TIMESTEP\n%ld\n" % 0);
        f.write("ITEM: NUMBER OF ATOMS\n%ld\n" % howmany);
        f.write("ITEM: BOX BOUNDS pp pp pp\n%.0f %.1f\n%.0f %.1f\n%.0f %.1f\n" % (0.,thisL,0.,thisL,0.,thisL));
        f.write("ITEM: ATOMS id type x y z\n");
        for i in range(howmany):
                f.write("%d %d %.5f %.5f %.5f 0 0 0\n" % (i+1,0, myinput[i,0], myinput[i,1], myinput[i,2]));

        f.close()
        
        

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

def assemble_unit(monomer, dendron): #assembles monomer linker and dendron
    
    out = np.empty((dendron.shape[0]+1+1,3)) #dendron.shape[0] gives number of rows 
    
    out[0,:] = monomer #writes monomer on first row of array out
    out[1,:] = monomer + [0,1,0] #writes monomer on second row of array out adding a y-translation
    #np.broadcast_to turns monomer+ y-translation to an array with shape similar to out but with two less
    #rows, basically this translated the dendron array to y=2 since in the function creating the dendron
    #the dendron starts from zero and here the backbone monomer starts at zero.
    out[2:,:] = dendron + np.broadcast_to(monomer + [0,2,0], (out.shape[0]-2,3)) 
    
    return out
    
################# Functions loaded ##########################################3    
Nmon = 500 #t(sys.argv[1]) # of monomers
gen = 1 ##int(sys.argv[2])   # generation
xbox =1.5*Nmon##float(sys.argv[3])
ybox = xbox #float(sys.argv[4])
zbox =xbox #float(sys.argv[5])
if os.path.isdir("G%d_%d"%(gen,Nmon))==False:
	os.mkdir("G%d_%d"%(gen,Nmon))

path_to_store="G%d_%d"%(gen,Nmon)
ldendron = 2**(gen+1)-1

nbonds_dendr = np.sum(np.asarray([2**i for i in range(1,gen+1)]))

allparticles = np.empty((Nmon*(2+(2**(gen+1)-1)),3), dtype = float) 

centerx = xbox/2.; centery = ybox/2.; centerz = zbox/2.- Nmon/2

dendron, bondmap= make_dendron(gen)
kk = 0
radius0 = 1.1
snapshot = hoomd.data.make_snapshot(N=allparticles[:,0].size,box=hoomd.data.boxdim(Lx=xbox, Ly=ybox, Lz=zbox),particle_types=['A', 'B','C'],bond_types=['backbone','linker','dendron'],angle_types=['backbone','linker'])
### Set Positions ####
for i in range(Nmon):
    radius = radius0
    myx = centerx 
    myy = centery 
    myz = centerz + radius*i #grows the chain each time by one radius
    monomer = np.asarray([myx,myy,myz]) #gives coordinates
    dendronized = assemble_unit(monomer, dendron) #calls function with monomer coordinates and dendron coordinates 
    #uses allparticles array's corrseponding row to copy the dendronized array 
    allparticles[i*(2+(2**(gen+1)-1)):(i+1)*(2+(2**(gen+1)-1)),:] =np.copy(dendronized)
    dendronized=dendronized-np.array([xbox/2,xbox/2,xbox/2])
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
#print (len(bondmap))
snapshot.bonds.resize(Nmon-1+2*Nmon+Nmon*len(bondmap))
#print("Bonds are: %d"%(Nmon-1+2*Nmon+Nmon*len(bondmap)))
kk=0
for k in range (0,Nmon-1):
    snapshot.bonds.group[kk]=[k*(2+(2**(gen+1)-1)),(k+1)*(2+(2**(gen+1)-1))]
    #print("Backbone:",(k)*(2+(2**(gen+1)-1)),(k+1)*(2+(2**(gen+1)-1)))
    snapshot.bonds.typeid[kk]=0
    kk+=1
for k in range (Nmon):
    snapshot.bonds.group[kk]=[k*(2+(2**(gen+1)-1))+1-1, k*(2+(2**(gen+1)-1))+2-1]
    #print("Linker 1")
    #print(k*(2+(2**(gen+1)-1))+1-1, k*(2+(2**(gen+1)-1))+2-1)
    snapshot.bonds.typeid[kk]=1
    kk+=1
    snapshot.bonds.group[kk]=[k*(2+(2**(gen+1)-1))+2-1, k*(2+(2**(gen+1)-1))+3-1]
   # print("Linker 2")
    #print(k*(2+(2**(gen+1)-1))+2-1, k*(2+(2**(gen+1)-1))+3-1)
    snapshot.bonds.typeid[kk]=1
    kk+=1
for k in range (Nmon):
    for w,l in bondmap:
        snapshot.bonds.group[kk]=[k*(2+(2**(gen+1)-1))+3+w-1,k*(2+(2**(gen+1)-1))+3+l-1];
 #       print("Dendron")
  #      print(k*(2+(2**(gen+1)-1))+3+w-1)
   #     print(k*(2+(2**(gen+1)-1))+3+l-1)
        snapshot.bonds.typeid[kk]=2;
        kk+=1
#print("Bonds Added:",kk)

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

#nl = hoomd.md.nlist.cell();
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

hoomd.md.integrate.mode_standard(dt=0.001);
to_run = 1e9
frames=100
hoomd.dump.gsd(filename=path_to_store+'/traj.gsd',period=int(to_run/frames), group=hoomd.group.all(),overwrite=True)
hoomd.analyze.log(filename=path_to_store+"/log-output.log",quantities=['potential_energy', 'temperature'],period=int(to_run/100) ,overwrite=True)

hoomd.run (to_run)

