import hoomd
import hoomd.md
import numpy as np 
import sys 
import os 
hoomd.context.initialize("");


Nmon =int(sys.argv[1]) # of monomers
gen = int(sys.argv[2])   # generation
xbox =1.5*Nmon##float(sys.argv[3])
ybox = xbox #float(sys.argv[4])
zbox =xbox #float(sys.argv[5])

ldendron = 2**(gen+1)-1

hoomd.init.read_gsd(filename=sys.argv[3],frame=-1)

nl=hoomd.md.nlist.tree()
### Pair Potential ##########
# wca=hoomd.md.pair.force_shifted_lj(r_cut=2**(1/6),nlist=nl)
# wca.pair_coeff.set('A','A',epsilon=1.0,sigma=1.0,r_cut=2**(1/6))
# wca.pair_coeff.set('B','B',epsilon=1.0,sigma=1.0,r_cut=2**(1/6))
# wca.pair_coeff.set('A','B',epsilon=1.0,sigma=1.0,r_cut=2**(1/6))
# wca.pair_coeff.set('A','C',epsilon=1.0,sigma=0.9,r_cut=0.9*2**(1/6))
# wca.pair_coeff.set('B','C',epsilon=1.0,sigma=0.9,r_cut=0.9*2**(1/6))
# wca.pair_coeff.set('C','C',epsilon=1.0,sigma=0.8,r_cut=0.8*2**(1/6))

# gauss = hoomd.md.pair.gauss(r_cut=2**(1/6.), nlist = nl)
# gauss.pair_coeff.set('A','A',epsilon=200.0,sigma=0.25,r_cut=2**(1/6))
# gauss.pair_coeff.set('B','B',epsilon=200.0,sigma=0.25,r_cut=2**(1/6))
# gauss.pair_coeff.set('A','B',epsilon=200.0,sigma=0.25,r_cut=2**(1/6))
# gauss.pair_coeff.set('A','C',epsilon=200.0,sigma=0.18,r_cut=0.9*2**(1/6))
# gauss.pair_coeff.set('B','C',epsilon=200.0,sigma=0.18,r_cut=0.9*2**(1/6))
# gauss.pair_coeff.set('C','C',epsilon=200.0,sigma=0.16,r_cut=0.8*2**(1/6))

morse = hoomd.md.pair.morse(r_cut=2**(1/6.), nlist = nl)
morse.pair_coeff.set('A','A',D0=1.,alpha=10.,r0=1.2,r_cut=2**(1/6))
morse.pair_coeff.set('B','B',D0=1.,alpha=10.,r0=1.2,r_cut=2**(1/6))
morse.pair_coeff.set('A','B',D0=1.,alpha=10.,r0=1.2,r_cut=2**(1/6))
morse.pair_coeff.set('A','C',D0=1.,alpha=10.,r0=1.08,r_cut=0.9*2**(1/6))
morse.pair_coeff.set('B','C',D0=1.,alpha=10.,r0=1.08,r_cut=0.9*2**(1/6))
morse.pair_coeff.set('C','C',D0=1.,alpha=10.,r0=0.97,r_cut=0.8*2**(1/6))

### Bonding Potential #######
# fene=hoomd.md.bond.fene()
# fene.bond_coeff.set('backbone',k=30.0,r0=1.5,sigma=1.0,epsilon=1.0)
# fene.bond_coeff.set('linker',k=30.0,r0=1.5,sigma=1.0,epsilon=1.0)
# fene.bond_coeff.set('dendron',k=30.0,r0=1.5,sigma=0.9,epsilon=1.0)

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
itg.set_gamma('A',gamma=1.0)
itg.set_gamma('B',gamma=1.0)
itg.set_gamma('C',gamma=1.0)
path_to_store='G%d_%d'%(gen,Nmon)
traj_fname='traj.gsd'
tidx=0
while os.path.isfile(path_to_store+'/'+traj_fname)==True:
	tidx+=1
	traj_fname='traj_%d.gsd'%(tidx)
hoomd.md.integrate.mode_standard(dt=0.005);
to_run=1e9
frames=100
hoomd.dump.gsd(filename=path_to_store+'/'+traj_fname,period=int(to_run/frames), group=hoomd.group.all(),overwrite=False)
hoomd.analyze.log(filename="log-output_%d.log"%(tidx),quantities=['potential_energy', 'temperature'],period=int(to_run/1000) ,overwrite=True)

hoomd.run (to_run)
#
