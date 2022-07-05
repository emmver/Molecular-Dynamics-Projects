
'''
This script can creates an initial configuration of a given number of rings, with all monomers of a single ring postioned on a circle.
The centers of the rings are put on simple cubic lattice points, under dilute conditions, to avoid overlaps.
Then bonds are added, except the closing bond, and hence, you get a linear polymer. 


Usage: /path/to/pypresso initial_poly.py /path/to/working_directory

The simulation consists of two steps:
1) Warm up: The force is capped and slowly ramped up over time. This is done because the circular conformation may lead 
to frustrated bonds that may break during the first steps of integration. This way, by force capping, bonds are not broken
and the confirmation slowly changes to a stable one. 

Here, the force cap, gamma and minimal distances are tracked and plotted as well as the energies.

2) Integration: Here the force cap has been removed and the simulation is a typical MD simulation
with a Langevin thermostat.

Here we keep track of the energies 

Important Variables
1) DP: Degree of polymerization
2) nchains: Number of chains
'''
import espressomd
espressomd.assert_features('WCA','DIPOLES')
from espressomd import thermostat
from espressomd import interactions
from espressomd import polymer
from espressomd.io.writer import vtf  # pylint: disable=import-error
import numpy as np
import logging
from espressomd import checkpointing
import espressomd.magnetostatics as magnetostatics
import sys
import warnings
from matplotlib import pyplot as plt
import matplotlib as mpl
import os 
import math
from scipy.stats import special_ortho_group
import time
from espressomd.observables import ParticlePositions
from espressomd.accumulators import Correlator

mpl.rcParams['figure.dpi'] = 400
plt.rcParams["font.family"] = "Ubuntu"
#plt.style.use('C:/Users/Toumba/Documents/plotstyle.mplstyle')
plt.style.use('~/plotstyle.mplstyle')

plt.rcParams['axes.linewidth'] = 4
plt.rcParams['xtick.major.size'] = 8 
plt.rcParams['ytick.major.size'] = 8 
plt.rcParams['xtick.labelsize']=20
plt.rcParams['ytick.labelsize']=20
plt.rcParams['xtick.minor.visible']=True
plt.rcParams['ytick.minor.visible']=True
plt.rcParams['xtick.minor.size'] = 5 
plt.rcParams['ytick.minor.size'] = 5
plt.rcParams['xtick.minor.width'] = 1.5
plt.rcParams['ytick.minor.width'] = 1.5
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['ytick.major.width'] = 2
plt.rcParams['xtick.major.pad']='8'
plt.rcParams['ytick.major.pad']='8'

np. set_printoptions(threshold=np. inf)


warnings.filterwarnings('ignore')

np.seterr(all="ignore")





def plot_min_dist_force(key):
    '''
    This function plots the minimal distance, force cap and gamma  , during the warmup step.
    key denotes the warmup or integration step
    '''
        plt.plot(steps,min_dist,'ro-',markerfacecolor=None)
        plt.ylabel("min_dist",fontsize=35)
        plt.yscale('log')
        #plt.xscale('log')
        #plt.legend(loc=(0.0,0.5),frameon=False)
        plt.xlabel("steps",fontsize=35)
        plt.savefig('./min_dist_{}.png'.format(key),dpi=300,bbox_inches='tight')
        plt.clf()
        plt.plot(steps,fcap,'ro-',markerfacecolor=None)
        plt.ylabel("force_cap",fontsize=35)
        plt.yscale('log')
        #plt.xscale('log')
        #plt.legend(loc=(0.0,0.5),frameon=False)
        plt.xlabel("steps",fontsize=35)
        plt.savefig('./force_cap__{}.png'.format(key),dpi=300,bbox_inches='tight')
        plt.clf()
        plt.plot(steps,gammas,'ro-',markerfacecolor=None)
        plt.ylabel("gamma",fontsize=35)
        plt.yscale('log')
        #plt.xscale('log')
        #plt.legend(loc=(0.0,0.5),frameon=False)
        plt.xlabel("steps",fontsize=35)
        plt.savefig('./gamma__{}.png'.format(key),dpi=300,bbox_inches='tight')
        plt.clf()

def plot_energies(key,plotlist,key_2):
    '''
    This function plots the list given through the plotlist variable.
    key: which energy is plotted e.g. bonded
    key_2: for which part. e.g active (magnetic) or passive
    '''
        plt.plot(steps,plotlist,'ro-',markerfacecolor=None)
        plt.ylabel("{} energy".format(key),fontsize=35)
       # plt.yscale('log')
        #plt.xscale('log')
        #plt.legend(loc=(0.0,0.5),frameon=False)
        plt.xlabel("steps",fontsize=35)
        plt.savefig(tostore_plots+'/energy_{}_{}.png'.format(key,key_2),dpi=300,bbox_inches='tight')
        plt.clf()


################################################################################################################

# System parameters
#############################################################
os.chdir(sys.argv[1])
log = open("myprog.log", "w")
sys.stdout = log

DP=200
chains=1
mpclist = [DP,DP]
maxmpc = max(mpclist)
nc = chains

#lattice constant
la = 3+maxmpc/3.14
n_rows = int(math.ceil(np.power(nc,1.0/3.0)))
print(n_rows)
L = 1.1*n_rows*la
print(L, L , L) 

'''
Here the system variables are set (box size, timestep, skin of neighbor list etc.)

See espressomd documentation for details: https://espressomd.github.io/doc4.1.4/system_setup.html
'''
system = espressomd.System(box_l=[L, L, L])
system.set_random_state_PRNG()
np.random.seed(seed=system.seed)

system.time_step = 0.01
system.min_global_cut=1.0
system.cell_system.skin = 0.4
system.thermostat.set_langevin(kT=1.0, gamma=1.0, seed=np.random.randint(1e6))
system.cell_system.set_domain_decomposition(use_verlet_lists=True)

#Pair potentials are set

system.non_bonded_inter[0,0].wca.set_params(
    epsilon=1, sigma=1)

system.non_bonded_inter[1,1].wca.set_params(
    epsilon=1, sigma=1)

system.non_bonded_inter[0,1].wca.set_params(
    epsilon=1, sigma=1)


fene = interactions.FeneBond(k=30, d_r_max=1.5)
system.bonded_inter.add(fene)


np.random.seed(int(time.time()))

'''
This double loop constructs the ring circular conformation and places the rings centers of mass on lattice points (on a simple cubic lattice).
This script was provided by Jan Smrek and was built to accomodate polydispersity in ring sizes. 
The variable mpclist is supposed to accomodate the different degrees of polymerization.
'''
for nca in [chains]:
    for dirnr in range(1):
        mpclist = [DP,DP]
        maxmpc = max(mpclist)
        nc = chains

        #lattice constant
        la = 3+maxmpc/3.14

        n_rows = int(math.ceil(np.power(nc,1.0/3.0)))
        print(n_rows)
        L = n_rows*la
        print(L, L , L) 

        Nlist = np.array([mpclist[0]]*int(nca) + [mpclist[1]]*int(nc-nca))
        print(Nlist) 
        totnp = sum(Nlist)
        print(totnp)
        '''
        This loop writes each monomer coordinate, for each ring (in the case of more than 1 rings), in the ring_primitive list. 
        The x and y coordinates follow a circle of radius r. 
        '''
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

        ch=0
        pid=0

        alist = np.random.choice(range(n_rows*n_rows*n_rows),n_rows*n_rows*n_rows,replace=False)
        chck = []
        print('alist',alist)
        '''
        This loop places the rings on the lattice points.
        '''
        for u in alist:
            i=u%n_rows 
            j=((u-i)/n_rows)%n_rows 
            k=(u-j*n_rows-i)/(n_rows*n_rows) 
            if [i,j,k] in chck:
                print("UAAAAAAAAA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            else:
                chck.append([i,j,k])
            if ch>=nc:
                
                break
            R = np.array([-L/2.0 + la/2.0 + i*la,-L/2.0 + la/2.0 + j*la, -L/2.0 + la/2.0 + k*la    ])                
            M = special_ortho_group.rvs(3) #random rotational matrix
            rp = [elem for elem in ring_primitive if len(elem)==Nlist[ch]][0]
            for ii in range(Nlist[ch]):
                rvec = R + np.dot(M,rp[ii])
                print("Adding Particle:",pid)
                system.part.add(id=pid, pos=[rvec[0],rvec[1],rvec[2]]);

                if ii==0:
                    pass # no closing bond	   
                else:
                    system.part[pid].add_bond((fene, pid-1)) # adding bonds
                pid+=1
            ch+=1


print("============= BONDS ==============")
outfile = open('polymer.vtf', 'w') #vtf files are used in VMD for visualization. Here the polymer.vtf file is created.

print("Number of particles:",DP*chains)
print("Example position",system.part[0].pos)
print("{} chains with {} monomers each".format(chains, DP))
print(system.part[:].id)
np.savetxt("pid.dat",system.part[:].id,fmt="%d")
vtf.writevsf(system, outfile) #vtf.writevsf writes in the polymer.vtf file the system structure i.e. number of particles, type of each particle, bonds between particles, angles between particles etc. 
vtf.writevcf(system, outfile) #vtf.writevtf writes in the particle coordinates x,y,z (unwrapped) for each frame.

npart=system.part[:].pos[:,0].size



#############################################################
#      Warmup                                               #
#############################################################

warm_steps = 10 # steps for each warm-up integration step (i.e. each )
wca_cap = 1 #force cap
system.force_cap = wca_cap
i = 0
act_min_dist = system.analysis.min_dist() # minimal distance between particles
tostore_plots="./passive_warmup"
if os.path.isdir(tostore_plots)==False:
    os.mkdir(tostore_plots)
gamma=1.0
# slowly ramp un up the cap
steps=[]
min_dist=[]
fcap=[]
count_steps=0
gammas=[]
while (act_min_dist < 0.95):
    energy=system.analysis.energy()
    gammas.append(gamma)
    vtf.writevcf(system, outfile)
    print("min_dist: {} \t force cap: {}".format(act_min_dist, wca_cap))
    system.part[:].v=[0,0,0]
    print("--------------------- Energies --------------------------")
    print(energy)
    system.integrator.run(warm_steps) # each time the following loop runs the integration moves forward by warm_steps #
    steps.append(count_steps)
    min_dist.append(system.analysis.min_dist())
    fcap.append(wca_cap)
    count_steps+=1
    act_min_dist = system.analysis.min_dist()
    wca_cap = wca_cap * 1.01
    system.force_cap = wca_cap
    plot_min_dist_force('warmup')

# remove force cap
wca_cap = 0
system.force_cap = wca_cap
system.integrator.run(warm_steps * 10)

log = open("myprog_afterwarmup.log", "w") #change logging file
sys.stdout = log

system.integrator.run(warm_steps * 10)
print("Finished warmup")

log = open("myprog_equil.log", "w") #change logging file
sys.stdout = log
#############################################################
#      Equilibration                                          #
#############################################################
tostore_plots="./passive_equil"
if os.path.isdir(tostore_plots)==False:
    os.mkdir(tostore_plots)  #folder to store plots for equilibration

# lists to store energies

steps=[]
total_en=[]
kin_en=[]
dipolar_en=[]
pair_en=[]
bond_en=[]
print("simulating...")
check_path="mycheck_passive"
checkpoint=checkpointing.Checkpoint(checkpoint_id=check_path,checkpoint_path='.')
checkpoint.register("system")
t_steps = 10000
check_idx=0
count=0
key_2='passive'
warm_steps=1000

for t in range(t_steps):
    energy=system.analysis.energy()
    print("step {} of {}".format(t, t_steps))    
    print("--------------------- Energies --------------------------")
    print(energy)
    system.integrator.run(warm_steps)
    count+=1
    steps.append(count*warm_steps)
    vtf.writevcf(system, outfile)
    total_en.append(energy['total']/npart)
    plot_energies('total',total_en,key_2)
    kin_en.append(2*energy['kinetic']/(3*npart))
    plot_energies('kinetic',kin_en,key_2)
    dipolar_en.append(energy['dipolar']/npart)
    plot_energies('dipolar',dipolar_en,key_2)
    pair_en.append(energy['non_bonded']/npart)
    plot_energies('non-bonded',pair_en,key_2)
    bond_en.append(energy['bonded']/npart)
    plot_energies('bonded',bond_en,key_2)
    if t%100==0:
        checkpoint.save(checkpoint_index=check_idx)  #save checkpoint every 100 integration steps
    check_idx+=1

outfile.close()
