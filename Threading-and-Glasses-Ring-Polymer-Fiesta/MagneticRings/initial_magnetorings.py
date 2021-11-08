
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
        for i in ragne (npart):
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

DP=200
chains=5
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
				#pid+=1

			#	data.append([pid,ch+1, 1, rvec[0],rvec[1],rvec[2]])
			#	data2.append([pid, 1, rvec[0],rvec[1],rvec[2]])
                system.part.add(id=pid, pos=[rvec[0],rvec[1],rvec[2]]);

				if ii==0:
					bonds.append([pid, 1, pid+Nlist[ch]-1,pid])
                    system.part[ pid+Nlist[ch]-1].add_bond((fene, pid))
				else:
					#bonds.append([pid, 1, pid-1,pid])
                    system.part[pid].add_bond((fene, pid+1))
               
                pid+=1
				


# System parameters
#############################################################

system = espressomd.System(box_l=[L, L, L])
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


fene = interactions.FeneBond(k=30, d_r_max=1.5)
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


print("{} chains with {} monomers each".format(nchains, monomers))



vtf.writevsf(system, outfile)
npart=system.part[:].pos[:,0].size