'''
This reads the checkpoint of an equilibrated passive system and adds magnetic dipoles in diffent forms, block, random or alternating copolymer.
It works for both linear and ring polymers. 
Use is: /path/to/pypresso restart_min_poly.py checkpoint_folder_name value_of_lambda working_directory

important variables:
1) DP : Degree of Polymerization
2) dip_moment: is actually the lambda value (proportional do dipole_moment^2)
3) how_many: How many particles will be magnetic. Expressed as a percentage of the total number of particles.
'''

import espressomd
espressomd.assert_features('WCA','DIPOLES')
from espressomd import thermostat
from espressomd import interactions
from espressomd import polymer
from espressomd.io.writer import vtf  # pylint: disable=import-error
import numpy as np
from espressomd import checkpointing
import espressomd.magnetostatics as magnetostatics
import sys
import warnings
from matplotlib import pyplot as plt
import matplotlib as mpl
import os 
import glob


##########################################################
mpl.rcParams['figure.dpi'] = 400
plt.rcParams["font.family"] = "Ubuntu"
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


def generate_random_list(times,max_n
    ''' Generates a random list of numbers '''
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
     '''
     This function adds dipole moments to a percentage of the particles.
     Needs as input the number of particles and the degree of polymerization, here N. 
     The key, which will dictate the way in which the dipole moments will be added. Options are
     'block', to create a diblock copolymer, 'copolymer' to create an alternating copolymer, 'random' to create 
     a random copolymer.
     '''
    nchains=int(npart/N)
    if key=='block':
        for j in range (nchains):
            idx=0
            for i in range (int(j*N),int((j+1)*N)):
                if idx<int(percent*N): 
                    rand_vect=np.random.rand(3) # random vector for the orientation of the dipole
                    dip_hat=rand_vect/np.linalg.norm(rand_vect) # unit vector for the orientation of the dipole
                    system.part[i].mol_id=int(np.floor(i/npart)) # id of molecule. useful in the case of more than one polymers.
                    system.part[i].type=1 # magnetic particle type 
                    system.part[i].dip=dip_hat # adding the unitary dipole vector to the system
                    system.part[i].dipm=moment**0.5 #adding a dipole moment value to the system, dictated by the lambda given initially.
                    system.part[i].rotation=[1,1,1] # activating rotation of each particle, it is turned off by default in espressomd.
                    idx+=1
                else:
                    system.part[i].mol_id=int(np.floor(i/npart)) # id of molecule. useful in the case of more than one polymers.
                    #system.part[i].type=1
                    system.part[i].rotation=[1,1,1]  # activating rotation of each particle, it is turned off by default in espressomd.

    elif key=='copolymer':
        for i in range (npart,2): 
            rand_vect=np.random.rand(3)  # random vector for the orientation of the dipole
            dip_hat=rand_vect/np.linalg.norm(rand_vect) # unit vector for the orientation of the dipole
            system.part[i].mol_id=int(np.floor(i/npart)) # id of molecule. useful in the case of more than one polymers.
            system.part[i].type=1 # magnetic particle type
            system.part[i].dip=dip_hat # adding the unitary dipole vector to the system
            system.part[j].dipm=moment**0.5 #adding a dipole moment value to the system, dictated by the lambda given initially.
            system.part[j].rotation=[1,1,1] # activating rotation of each particle, it is turned off by default in espressomd.
        for i in ragne (npart):
            system.part[i].mol_id=int(np.floor(i/npart))
            system.part[i].rotation=[1,1,1]

    elif key=='random':
        nchains=int(npart/N)
        for i in range(nchains):
            random_list=generate_random_list(percent*N,i*N) #random list to add dipoles into completely random monomers
            for j in random_list:
                rand_vect=np.random.rand(3)  # random vector for the orientation of the dipole
                dip_hat=rand_vect/np.linalg.norm(rand_vect)
                system.part[j].mol_id=int(np.floor(i/npart))
                system.part[j].type=1
                #system.part[j].dip=dip_hat
                system.part[j].dipm=moment**0.5
        for i in range (npart):
            system.part[i].mol_id=int(np.floor(i/npart))
            system.part[i].rotation=[1,1,1]



def list_setup():
   '''Sets up lists for storing data'''
                         
    count=0
    int_cycl=50
    steps=[]
    #min_dist=[]
    #fcap=[]
    count_steps=0
    #gammas=[]
    total_en=[]
    kin_en=[]
    dipolar_en=[]
    pair_en=[]
    bond_en=[]
    ken_0=[]
    ken_1=[]

    listoflists=[steps,total_en,kin_en,dipolar_en,pair_en,bond_en,ken_0,ken_1]
    string_lists=['steps','total_en','kin_en','dipolar_en','pair_en','bond_en','ken_0','ken_1']
    return count,steps,count_steps,total_en,kin_en,dipolar_en,pair_en,bond_en,ken_0,ken_1,listoflists,string_lists



def plot_min_dist_force(key):
         ''' plots minimal distance, force cap and gamma as a function of simulation step '''
        plt.plot(steps,min_dist,'ro-',markerfacecolor=None)
        plt.ylabel("min_dist",fontsize=35)
        plt.yscale('log')
        #plt.xscale('log')
        #plt.legend(loc=(0.0,0.5),frameon=False)
        plt.xlabel("steps",fontsize=35)
        plt.savefig(tostore_plots+'/min_dist_{}.png'.format(key),dpi=300,bbox_inches='tight')
        plt.clf()
        plt.plot(steps,fcap,'ro-',markerfacecolor=None)
        plt.ylabel("force_cap",fontsize=35)
        plt.yscale('log')
        #plt.xscale('log')
        #plt.legend(loc=(0.0,0.5),frameon=False)
        plt.xlabel("steps",fontsize=35)
        plt.savefig(tostore_plots+'/force_cap__{}.png'.format(key),dpi=300,bbox_inches='tight')
        plt.clf()
        plt.plot(steps,gammas,'ro-',markerfacecolor=None)
        plt.ylabel("gamma",fontsize=35)
        plt.yscale('log')
        #plt.xscale('log')
        #plt.legend(loc=(0.0,0.5),frameon=False)
        plt.xlabel("steps",fontsize=35)
        plt.savefig(tostore_plots+'/gamma__{}.png'.format(key),dpi=300,bbox_inches='tight')
        plt.clf()

def plot_energies(key,plotlist,key_2):
        ''' plots energies as function of simulation step''' 
        plt.plot(steps,plotlist,'ro-',markerfacecolor=None)
        plt.ylabel("{} energy".format(key),fontsize=35)
       # plt.yscale('log')
        #plt.xscale('log')
        #plt.legend(loc=(0.0,0.5),frameon=False)
        plt.xlabel("steps",fontsize=35)
        plt.savefig(tostore_plots+'/energy_{}_{}.png'.format(key,key_2),dpi=300,bbox_inches='tight')
        plt.clf()



def write_to_vtk(vtk_idx):
     ''' 
     Writes vtk file containing all frames. vtk files contain vector information for dipoles, velocities as well as positions.
     Another script is needed to 'split' this all.vtk to separate frames. 
     '''
    store_vtk_file=tostore_vtk+"/part_test_{}.vtk".format(vtk_idx)
    pos_dummy=tostore_vtk+"/pos_dummy.txt"
    dip_dummy=tostore_vtk+"/dipoles_dummy.txt"
    bound_idx_dummy=tostore_vtk+"/bidx_dummy.txt"
    vel_dummy=tostore_vtk+"/veloc_dummy.txt"
    np.savetxt(pos_dummy,system.part[:].pos)
    np.savetxt(bound_idx_dummy,(system.part[:].pos/system.box_l))
    np.savetxt(vel_dummy,system.part[:].v)
    np.savetxt(dip_dummy,system.part[:].dip)
    os.system("cat ./pos_h_file.txt " +pos_dummy+ " ./bidx_h_file.txt "+ bound_idx_dummy + " ./vel_h_file.txt "+ vel_dummy +" ./dipoles_h_file.txt " + dip_dummy +" > " + store_vtk_file )
    os.system("cat " + store_vtk_file + " > ./all.vtk " )


def vtk_heads():
     ''' Writes header files for each vtk sector '''
    f=open ("pos_h_file.txt","w");
    f.write('# vtk DataFile Version 2.0')
    f.write("\n")
    f.write('particles')
    f.write('\n')
    f.write('ASCII')
    f.write('\n')
    f.write('DATASET UNSTRUCTURED_GRID')
    f.write('\n')
    f.write('POINTS 200 floats')
  #  f.write('\n')
#    f.write("POSITIONS positions float")
 #   f.write("\n")
    f.close()


    f=open ("bidx_h_file.txt","w");
    f.write('POINT_DATA 200')
    f.write("\n")

    f.write('SCALARS pbc_idx float 3')
    f.write("\n")

    f.write('LOOKUP_TABLE default')
    f.write("\n")

    f.close()

    f=open ("vel_h_file.txt","w");
    f.write('POINT_DATA 200')
    f.write("\n")
    f.write('SCALARS velocity float 3')
    f.write("\n")
    f.write('LOOKUP_TABLE default')
    f.write("\n")
    f.close()


    f=open ("dipoles_h_file.txt","w");
    f.write("\n\n")
    f.write("VECTORS dipoles float")
    f.write("\n\n")
    f.close()

def store_all_lists():
    ''' Stores all lists '''
    for i in range (len(listoflists)):
        np.savetxt(tostore_plots+"/"+string_lists[i],listoflists[i])

os.chdir(sys.argv[3])
check_folder='../mycheck_passive'
checks=glob.glob(check_folder+'/*')
checks.sort(key=os.path.getmtime)
print(checks)
checkpoint=checkpointing.Checkpoint(checkpoint_id=sys.argv[1],checkpoint_path='..') # finding checkpoints
name=checks[-1].split('/') # gettin last checkpoint
print('Name:',name)
print(name[2][0:4])
last_idx=int(name[2][0:4])
print("last idx:",last_idx)
checkpoint.load(checkpoint_index=last_idx) # loading last checkpoint
DP=200 #degree of polymerization

epsilon=1; sigma=1 #### Define WCA Params ####

##############################################################
#      Magnetics on                                          #
##############################################################


log = open("myprog_magnetics_on.log", "w") #log file
sys.stdout = log

print("Positions=\n {}".format(system.part[:].pos))
print("=================== =================== =====================")

print ("L: ", system.box_l) 
print("=================== =================== =====================")

print("number density",system.part[:].pos[:,0].size/system.box_l**3)
print("=================== =================== =====================")

print("periodicity [x y z]:",system.periodicity)
print("=================== =================== =====================")

print("Particles:",system.part[:].pos[:,0].size)
print("=================== =================== =====================")


print("#### Temp Test ######")

print("Temperatures=\n {}".format(system.part[:].temp))
print("=================== =================== =====================")

print("#### Non Bonded Test 0,0 ######")
print("Non-Bonded:\n {}".format(system.non_bonded_inter[0,0].wca.get_params()))
print("#### Non Bonded Test 1,1 ######")
print("Non-Bonded:\n {}".format(system.non_bonded_inter[1,1].wca.get_params()))
print("#### Non Bonded Test 0,1 ######")
print("Non-Bonded:\n {}".format(system.non_bonded_inter[0,1].wca.get_params()))
print("=================== =================== =====================")


print("#### Non Bonded Test 0,0 ######")
print("Non-Bonded:\n {}".format(system.non_bonded_inter[0,0].wca.get_params()))
print("#### Non Bonded Test 1,1 ######")
print("Non-Bonded:\n {}".format(system.non_bonded_inter[1,1].wca.get_params()))
print("#### Non Bonded Test 0,1 ######")
print("Non-Bonded:\n {}".format(system.non_bonded_inter[0,1].wca.get_params()))
print("=================== =================== =====================")
system.non_bonded_inter[0,0].wca.get_params()
print("\n### system.thermostat test ###")
print("system.thermostat.get_state() = {}".format(system.thermostat.get_state()))



outfile = open('polymer_magnetic.vtf', 'w') # writing vtf outfile

dip_howmany=0.5 # ratio magnetic/all particles
print("number of magnetic particles:",int(system.part[:].pos[:,0].size*dip_howmany))
dip_moment=float(sys.argv[2]) # lambda ~ dipole_moment^2
npart=system.part[:].pos[:,0].size
###### Adding dipole moments ####################
add_activity(system.part[:].pos[:,0].size,DP,'block',dip_howmany,dip_moment)
print("\n### system.dip test ###")
print("system.part[:].dip = {}".format(system.part[:].dip))

p3m = magnetostatics.DipolarP3M(prefactor=1,accuracy=1.2e-3) # p3m solver for dipolar interactions

system.actors.add(p3m)

print("DONE!")
tostore_vtk="vtk_warmup" 
vtk_idx=0
if os.path.isdir(tostore_vtk):
    pass 
else:
    os.mkdir(tostore_vtk) ## store VTK files

tostore_plots='magnetic_plots'
if os.path.isdir(tostore_plots):
    pass 
else:
    os.mkdir(tostore_plots) # store plots 
#############################################################
#      Check params                                               #
#############################################################
log = open("myprog_magnetics_warmup.log", "w") #changing log files
sys.stdout = log
print("Start Warmup after magnetics",flush=True)
print("================ ============== ===================",flush=True)
warm_steps = 10
i = 0

gamma=1e0
tstep=5e-3
system.time_step = tstep #setting time-step --> Change from initial in checkpoint

system.thermostat.set_langevin(kT=1.0, gamma=gamma) #--> Setting thermostat

# Setup lists for storing data
count,steps,count_steps,total_en,kin_en,dipolar_en,pair_en,bond_en,ken_0,ken_1,listoflists,string_lists=list_setup() ##setting up lists
vtf.writevsf(system, outfile)
vtf.writevcf(system, outfile)
energy=system.analysis.energy()
system.integrator.set_vv() 

key_2='magnetic_warmup'

active=system.part.select(type=1) #splitting active (magnetic) particles from passive ones.
passive=system.part.select(type=0)
vtk_heads()
f=open ("vectors_h_file.txt","w");
f.write("\n\n")
f.write("VECTORS vectors float")
f.write("\n\n")
f.close()
for t in range (1000):
    energy=system.analysis.energy()
    print("--------------------- Energies --------------------------",flush=True)
    print(energy)
    print("Timestep is: ",system.time_step)
    vtf.writevcf(system, outfile)
    system.part.writevtk("./dummy.vtk");
    write_to_vtk(vtk_idx);vtk_idx+=1
    system.integrator.run(warm_steps) #integrating
    count_steps+=1
    store_all_lists()

    count+=1
    steps.append(count_steps*warm_steps)
    total_en.append(energy['total']/npart)
    plot_energies('total',total_en,key_2)
    kin_en.append(energy['kinetic']/(1.5*npart))
    plot_energies('kinetic',kin_en,key_2)
    dipolar_en.append(energy['dipolar']/npart)
    plot_energies('dipolar',dipolar_en,key_2)
    pair_en.append(energy['non_bonded']/npart)
    plot_energies('non-bonded',pair_en,key_2)
    bond_en.append(energy['bonded']/npart)
    plot_energies('bonded',bond_en,key_2)
    vel=np.linalg.norm(passive.v,axis=1)
    kinetic_passive=(2/3)*(0.5*passive.mass*vel**2)
    ken_0.append(kinetic_passive.mean())
    plot_energies('T_0',ken_0,key_2)
    vel=np.linalg.norm(active.v,axis=1)
    kinetic_active=(2/3)*(0.5*active.mass*vel**2)
    ken_1.append(kinetic_active.mean())
    plot_energies('T_1',ken_1,key_2)

print("======= ======== =======")
print("WARMUP DONE!!")
store_all_lists()
#### Up to here, equilibrating system without bending potential to avoid chain breaking due to dipolar interactions ###
angle_harmonic = interactions.AngleHarmonic(bend=1.0, phi0=2 * np.pi / 3) ##adding harmonic bending otetnial
system.bonded_inter.add(angle_harmonic)
nchains=int(system.part[:].pos[:,0].size/DP)
monomers=DP
for j in range (nchains):
    for i in range(int(j*monomers),int((j+1)*monomers)-2):
        id=i
        system.part[id+1].add_bond((angle_harmonic,id,id+2))



# restore simulation temperature
system.integrator.run(warm_steps * 100) ### integrate for a few steps with angular potential on ###
print("Finished warmup")




#############################################################
#      Integration                                          #
#############################################################
log = open("myprog_magnetics_int.log", "w") # changing log
sys.stdout = log
system.integrator.set_vv() 
key_2="magnetics_equil"
print("simulating...")
check_path="mycheck_magnetic"
checkpoint=checkpointing.Checkpoint(checkpoint_id=check_path,checkpoint_path='.')
checkpoint.register("system") # storing checkpoints with dipolar interactions.
t_steps = 10000 # time steps
frames=1000 #frames to store
store_step=int(t_steps/frames)
check_idx=0
checkpoint.save(checkpoint_index=check_idx)
warm_steps=10000 # integration steps for each time step

## Setup lists for data storage as well as indices
count,steps,count_steps,total_en,kin_en,dipolar_en,pair_en,bond_en,ken_0,ken_1,listoflists,string_lists=list_setup()
vtk_heads()

############Store folders for vtk and plots##################
tostore_vtk="vtk_equil"                                     #
if os.path.isdir(tostore_vtk):                              #
    pass                                                    #
else:                                                       #
    os.mkdir(tostore_vtk)                                   #
tostore_plots='magnetic_plots_equil'                        #
if os.path.isdir(tostore_plots):                            #
    pass                                                    #
else:                                                       #
    os.mkdir(tostore_plots)                                 #
vtk_idx=0                                                   #
#############################################################

for t in range(t_steps):
    energy=system.analysis.energy()
    print("step {} of {}".format(t, t_steps))
    
    print("--------------------- Energies --------------------------")
    print(energy)
    
    system.integrator.run(warm_steps)
    vtf.writevcf(system, outfile)

    if t%store_step==0:
        system.part.writevtk("./dummy.vtk");
        write_to_vtk(vtk_idx);vtk_idx+=1

    steps.append(t*warm_steps)

    total_en.append(energy['total']/npart)
    plot_energies('total',total_en,key_2)
    kin_en.append(energy['kinetic']/(1.5*npart))
    plot_energies('kinetic',kin_en,key_2)
    dipolar_en.append(energy['dipolar']/npart)
    plot_energies('dipolar',dipolar_en,key_2)
    pair_en.append(energy['non_bonded']/npart)
    plot_energies('non-bonded',pair_en,key_2)
    bond_en.append(energy['bonded']/npart)
    plot_energies('bonded',bond_en,key_2)
    vel=np.linalg.norm(passive.v,axis=1)
    kinetic_passive=(2/3)*(0.5*passive.mass*vel**2)
    ken_0.append(kinetic_passive.mean())
    plot_energies('T_0',ken_0,key_2)
    vel=np.linalg.norm(active.v,axis=1)
    kinetic_active=(2/3)*(0.5*active.mass*vel**2)
    ken_1.append(kinetic_active.mean())
    plot_energies('T_1',ken_1,key_2)

#    check_idx+=1
 #   if t%store_step==0:
  #      checkpoint.save(checkpoint_index=check_idx)
store_all_lists()

outfile.close()
