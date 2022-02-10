
"""
This reads the checkpoint of an equilibrated passive system and adds magnetic dipoles in diffent forms, block, random or alternating copolymer"
"""
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
                   # system.part[i].dipm=moment
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
            #system.part[j].dipm=moment
            system.part[j].rotation=[1,1,1]
        for i in ragne (npart):
            system.part[i].mol_id=int(np.floor(i/npart))
            system.part[i].rotation=[1,1,1]

    elif key=='random':
        nchains=int(npart/N)
        for i in range(nchains):
            random_list=generate_random_list(percent*N,i*N)
            for j in random_list:
                dip_phi = np.random.random(()) *2. * np.pi
                dip_cos_theta = 2*np.random.random(()) -1
                dip_sin_theta = np.sin(np.arccos(dip_cos_theta))
                rand_vect=np.array([dip_sin_theta *np.sin(dip_phi),dip_sin_theta *np.cos(dip_phi),dip_cos_theta])
                dip_hat=rand_vect/np.linalg.norm(rand_vect)
                system.part[j].mol_id=int(np.floor(i/npart))
                system.part[j].type=1
                system.part[j].dip=dip_hat
                #system.part[j].dipm=moment
        for i in range (npart):
            system.part[i].mol_id=int(np.floor(i/npart))
            system.part[i].rotation=[1,1,1]



    system.non_bonded_inter[0,0].wca.set_params(
    epsilon=epsilon, sigma=sigma)

    system.non_bonded_inter[1,1].wca.set_params(
    epsilon=epsilon, sigma=sigma)

    system.non_bonded_inter[0,1].wca.set_params(
    epsilon=epsilon, sigma=sigma)



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

DP=10
chains=1
mpclist = [DP,DP]#[400,800]#
maxmpc = max(mpclist)
nc = chains

#lattice constant
la = 3+maxmpc/3.14
n_rows = int(math.ceil(np.power(nc,1.0/3.0)))
print(n_rows)
L = 1.1*n_rows*la
print(L, L , L) 


system = espressomd.System(box_l=[L, L, L])
system.set_random_state_PRNG()
#system.seed = system.cell_system.get_state()['n_nodes'] * [1234]
np.random.seed(seed=system.seed)

system.time_step = 0.01
system.min_global_cut=1.0
system.cell_system.skin = 0.4
system.thermostat.set_langevin(kT=1.0, gamma=1.0, seed=np.random.randint(1e6))
system.cell_system.set_domain_decomposition(use_verlet_lists=True)

system.non_bonded_inter[0,0].wca.set_params(
    epsilon=1, sigma=1)

system.non_bonded_inter[1,1].wca.set_params(
    epsilon=1, sigma=1)

system.non_bonded_inter[0,1].wca.set_params(
    epsilon=1, sigma=1)


fene = interactions.FeneBond(k=30, d_r_max=1.5)
system.bonded_inter.add(fene)


np.random.seed(int(time.time()))#22378)
for nca in [chains]:#[10,30,50,70,90,100,110,130,150,170,190]:
    for dirnr in range(1):
        mpclist = [DP,DP]#[400,800]#
        maxmpc = max(mpclist)
        nc = chains

        #lattice constant
        la = 3+maxmpc/3.14

        n_rows = int(math.ceil(np.power(nc,1.0/3.0)))
        print(n_rows)
        L = n_rows*la
        print(L, L , L) 

        Nlist = np.array([mpclist[0]]*int(nca) + [mpclist[1]]*int(nc-nca))
        #Nlist = np.random.choice(Nlist,len(Nlist),replace=False) #random order of big and small rings
        print(Nlist) 
        #np.savetxt(outpath+'chain_lengths.dat',Nlist,fmt='%d')
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
        print('alist',alist)
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
            R = np.array([-L/2.0 + la/2.0 + i*la,-L/2.0 + la/2.0 + j*la, -L/2.0 + la/2.0 + k*la    ])                
            M = special_ortho_group.rvs(3) #random rotational matrix
            rp = [elem for elem in ring_primitive if len(elem)==Nlist[ch]][0]
            for ii in range(Nlist[ch]):
                rvec = R + np.dot(M,rp[ii])
            #    data.append([pid,ch+1, 1, rvec[0],rvec[1],rvec[2]])
            #    data2.append([pid, 1, rvec[0],rvec[1],rvec[2]])
                print("Adding Particle:",pid)
                system.part.add(id=pid, pos=[rvec[0],rvec[1],rvec[2]]);

                if ii==0:
                    bonds.append([pid,pid+Nlist[ch-1]-1])
                    #bonds.append([pid, 1, pid+Nlist[ch]-1,pid])
                    system.part[pid].add_bond((fene, pid+Nlist[ch]-1))
                else:
                    #bonds.append([pid, 1, pid-1,pid])
                    bonds.append([pid,pid+1])
                    system.part[pid].add_bond((fene, pid-1))
               
                pid+=1
            ch+=1


print("============= BONDS ==============")
print(bonds)
outfile = open('polymer.vtf', 'w')

#energy=system.analysis.energy()
print("Number of particles:",DP*chains)
print("Example position",system.part[0].pos)
#print("Example position 2",system.part[240].pos)
#print("example 3", system.part[991].pos)
print("{} chains with {} monomers each".format(chains, DP))
print(system.part[:].id)
np.savetxt("pid.dat",system.part[:].id,fmt="%d")
vtf.writevsf(system, outfile)
vtf.writevcf(system, outfile)

npart=system.part[:].pos[:,0].size



#############################################################
#      Warmup                                               #
#############################################################

warm_steps = 10
wca_cap = 1
system.force_cap = wca_cap
i = 0
act_min_dist = system.analysis.min_dist()
tostore_plots="./passive_warmup"
if os.path.isdir(tostore_plots)==False:
    os.mkdir(tostore_plots)
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
tostore_plots="./passive_equil"
if os.path.isdir(tostore_plots)==False:
    os.mkdir(tostore_plots)
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
    if t%100==0:
        checkpoint.save(checkpoint_index=check_idx)
    check_idx+=1

outfile.close()
