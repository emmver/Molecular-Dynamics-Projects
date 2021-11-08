#
# Copyright (C) 2013-2019 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
"""
This sample sets up a polymer.
"""
import espressomd
espressomd.assert_features(["WCA"])
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


mpl.rcParams['figure.dpi'] = 400
plt.rcParams["font.family"] = "Ubuntu"
#plt.style.use('C:/Users/Toumba/Documents/plotstyle.mplstyle')
plt.style.use('/home/manosver/plotstyle.mplstyle')

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
#######################################




warnings.filterwarnings('ignore')

np.seterr(all="ignore")


log = open("myprog.log", "w")
sys.stdout = log


def checkIfDuplicates_1(listOfElems):
    ''' Check if given list contains any duplicates '''
    if len(listOfElems) == len(set(listOfElems)):
        return False
    else:
        return True


def generate_random_list(times,max_n):
    out=[]
    times=int(times)
    max_n=int(max_n)
    for i in range(times):
        out.append(np.random.randint(max_n))
    result = checkIfDuplicates_1(out)
    while result==True:
        print('Yes, list contains duplicates')
        out=[]
        for i in range(times):
            out.append(np.random.randint(max_n))
    
    print('No duplicates found in list')
    return out


def add_activity(npart,N,key,percent,moment): 
    nchains=int(npart/N)
    if key=='block':
        for j in range (nchains):
            idx=0
            for i in range (int(j*N),int((j+1)*N)):
                if idx<int(percent*N): 
                    rand_vect=np.random.rand(3)
                    dip_hat=rand_vect/np.linalg.norm(rand_vect)
                    system.part[i].mol_id=int(np.floor(i/npart))
                    system.part[i].type=1
                    system.part[i].dip=dip_hat
                    system.part[i].dipm=moment
                    system.part[i].rotation=[1,1,1]
                    idx+=1
                else:
                    system.part[i].mol_id=int(np.floor(i/npart))
                    #system.part[i].type=1
                    system.part[i].rotation=[1,1,1]

    elif key=='copolymer':
        for i in range (npart,2):
            rand_vect=np.random.rand(3)
            dip_hat=rand_vect/np.linalg.norm(rand_vect)
            system.part[i].mol_id=int(np.floor(i/npart))
            system.part[i].type=1
            system.part[i].dip=dip_hat
            system.part[j].dipm=moment
            system.part[j].rotation=[1,1,1]
        for i in ragne (npart):
            system.part[i].mol_id=int(np.floor(i/npart))
            system.part[i].rotation=[1,1,1]

    elif key=='random':
        nchains=int(npart/N)
        for i in range(nchains):
            random_list=generate_random_list(percent*N,i*N)
            for j in random_list:
                rand_vect=np.random.rand(3)
                dip_hat=rand_vect/np.linalg.norm(rand_vect)
                system.part[j].mol_id=int(np.floor(i/npart))
                system.part[j].type=1
                system.part[j].dip=dip_hat
                system.part[j].dipm=moment
        for i in ragne (npart):
            system.part[i].mol_id=int(np.floor(i/npart))
            system.part[i].rotation=[1,1,1]

def plot_min_dist_force(key):
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
        plt.plot(steps,plotlist,'ro-',markerfacecolor=None)
        plt.ylabel("{} energy".format(key),fontsize=35)
        plt.yscale('log')
        #plt.xscale('log')
        #plt.legend(loc=(0.0,0.5),frameon=False)
        plt.xlabel("steps",fontsize=35)
        plt.savefig('./energy_{}_{}.png'.format(key,key_2),dpi=300,bbox_inches='tight')
        plt.clf()

# System parameters
#############################################################

system = espressomd.System(box_l=[80, 80, 80])
system.set_random_state_PRNG()
#system.seed = system.cell_system.get_state()['n_nodes'] * [1234]
np.random.seed(seed=system.seed)

system.time_step = 0.01
system.min_global_cut=1.0
system.cell_system.skin = 0.4
system.thermostat.set_langevin(kT=1.0, gamma=1.0, seed=42)
system.cell_system.set_domain_decomposition(use_verlet_lists=True)

outfile = open('polymer.vtf', 'w')


system.non_bonded_inter[0,0].wca.set_params(
    epsilon=1, sigma=1)

system.non_bonded_inter[1,1].wca.set_params(
    epsilon=1, sigma=1)

system.non_bonded_inter[0,1].wca.set_params(
    epsilon=1, sigma=1)



fene = interactions.FeneBond(k=30, d_r_max=1.5,r_0=1.0)
system.bonded_inter.add(fene)


energy=system.analysis.energy()
positions = polymer.positions(n_polymers=5,
                              beads_per_chain=200,
                              bond_length=1.0,
                              seed=3210)
print("Number of particles:",positions[:,:,0].size)
print("Example position",positions[0,:,0].shape[0])
print("Example position 2",positions[:,0].shape[0])
print("example 3", positions[0,2,:])
nchains=positions[:,0].shape[0]
monomers=positions[0].shape[0]
DP=monomers

print("{} chains with {} monomers each".format(nchains, monomers))
for j in range (nchains): #chains
    idx=0
    for i in range (int(j*monomers),int((j+1)*monomers)): #monomers per chain
#        print("Particle ID:",int(i))
        system.part.add(id=i, pos=positions[j,idx,:]);idx+=1

for j in range (nchains):
    for i in range (int(j*monomers),int((j+1)*monomers)-1):
        system.part[i].add_bond((fene, i+1))


vtf.writevsf(system, outfile)
npart=system.part[:].pos[:,0].size


#############################################################
#      Warmup                                               #
#############################################################

warm_steps = 10
wca_cap = 1
system.force_cap = wca_cap
i = 0
act_min_dist = system.analysis.min_dist()

# warmup with zero temperature to remove overlaps
system.thermostat.set_langevin(kT=0.0, gamma=1.0)
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
    #print("total energy: {} \t kinetic energy: {} \t bonded energy: {} \t non-bonded energy: {} \t".format(energy["total"],energy["kinetic"],energy["bonded"],energy["non_bonded"]))
    system.integrator.run(warm_steps)
    #system.part[:].v = [0, 0, 0]
    # Warmup criterion
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

log = open("myprog_afterwarmup.log", "w")
sys.stdout = log

system.thermostat.set_langevin(kT=1.0, gamma=1.0)
print("#### Temp Test ######")
print("Temperatures=\n {}".format(system.part[:].temp))
system.integrator.run(warm_steps * 10)
print("Finished warmup")

log = open("myprog_equil.log", "w")
sys.stdout = log
#############################################################
#      Equilibration                                          #
#############################################################
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
t_steps = 100
check_idx=0
count=0
key_2='passive'
for t in range(t_steps):
    energy=system.analysis.energy()
    print("step {} of {}".format(t, t_steps))
    #print("total energy: {} \t kinetic energy: {} \t bonded energy: {} \t non-bonded energy: {} \t".format(energy["total"],energy["kinetic"],energy["bonded"],energy["non_bonded"]))
    
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
    if t%10==0:
        checkpoint.save(checkpoint_index=check_idx)
    check_idx+=1

outfile.close()
