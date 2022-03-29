
"""
This reads the checkpoint of an equilibrated dipolar system and adds a magnetic field"
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
import glob 
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
                dip_phi = np.random.random(()) *2. * np.pi
                dip_cos_theta = 2*np.random.random(()) -1
                dip_sin_theta = np.sin(np.arccos(dip_cos_theta))
                rand_vect=np.array([dip_sin_theta *np.sin(dip_phi),dip_sin_theta *np.cos(dip_phi),dip_cos_theta])
                dip_hat=rand_vect/np.linalg.norm(rand_vect)
                system.part[j].mol_id=int(np.floor(i/npart))
                system.part[j].type=1
                system.part[j].dip=dip_hat
                system.part[j].dipm=moment
        for i in range (npart):
            system.part[i].mol_id=int(np.floor(i/npart))
            system.part[i].rotation=[1,1,1]



def list_setup():
    # Setup lists for storing data
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
        plt.plot(steps,plotlist,'ro-',markerfacecolor=None)
        plt.ylabel("{} energy".format(key),fontsize=35)
       # plt.yscale('log')
        #plt.xscale('log')
        #plt.legend(loc=(0.0,0.5),frameon=False)
        plt.xlabel("steps",fontsize=35)
        plt.savefig(tostore_plots+'/energy_{}_{}.png'.format(key,key_2),dpi=300,bbox_inches='tight')
        plt.clf()


def write_to_vtk(vtk_idx):
    store_vtk_file=tostore_vtk+"/part_test_{}.vtk".format(vtk_idx)
    dip_dummy=tostore_vtk+"/dipoles_dummy.txt"
    np.savetxt(tostore_vtk+"/dipoles_dummy.txt",system.part[:].dip)
    os.system("cat ./dummy.vtk ./vectors_h_file.txt " + dip_dummy + " > " + store_vtk_file )

def store_all_lists():
    for i in range (len(listoflists)):
        np.savetxt(tostore_plots+"/"+string_lists[i],listoflists[i])

os.chdir(sys.argv[3])
check_folder='../mycheck_passive'
checks=glob.glob(check_folder+'/*')
checks.sort(key=os.path.getmtime)
print(checks)
checkpoint=checkpointing.Checkpoint(checkpoint_id=sys.argv[1],checkpoint_path='..')
name=checks[-1].split('/')
print('Name:',name)
print(name[2][0:4])
last_idx=int(name[2][0:4])
print("last idx:",last_idx)
checkpoint.load(checkpoint_index=last_idx)

epsilon=1; sigma=1 #### Define WCA Params ####
# magnetic field times dipole moment
dummy=np.random.randn(3)
H_dipm = 1
H_field = H_dipm*(dummy/np.linalg.norm(dummy))
##############################################################
#      Restart from lst frame of magnetic snapshot                                         #
##############################################################


log = open("myprog_magnetics_on.log", "w")
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


print("\n### system.thermostat test ###")
print("system.thermostat.get_state() = {}".format(system.thermostat.get_state()))

print("\n### system.dip test ###")
print("system.part[:].dip = {}".format(system.part[:].dip))


outfile = open('polymer_dipoled_field.vtf', 'w')


p3m = magnetostatics.DipolarP3M(prefactor=1,accuracy=1.2e-3)

system.actors.add(p3m)

print("DONE!")
tostore_vtk="vtk_warmup_field"
vtk_idx=0
if os.path.isdir(tostore_vtk):
    pass 
else:
    os.mkdir(tostore_vtk)

tostore_plots='magnetic_plots_field'
if os.path.isdir(tostore_plots):
    pass 
else:
    os.mkdir(tostore_plots)
#############################################################
#      Check params                                               #
#############################################################
log = open("myprog_magnetics_field_warmup.log", "w")
sys.stdout = log
print("Start Warmup after magnetics equil - Now applying field",flush=True)
print("================ ============== ===================",flush=True)
warm_steps = 10
i = 0

gamma=1e0
tstep=5e-3
system.time_step = tstep

system.thermostat.set_langevin(kT=1.0, gamma=gamma)

# Setup lists for storing data
count,steps,count_steps,total_en,kin_en,dipolar_en,pair_en,bond_en,ken_0,ken_1,listoflists,string_lists=list_setup()
vtf.writevsf(system, outfile)
vtf.writevcf(system, outfile)
energy=system.analysis.energy()
system.integrator.set_vv() 

key_2='magnetic_warmup'

active=system.part.select(type=1)
passive=system.part.select(type=0)

f=open ("vectors_h_file.txt","w");
f.write("\n\n")
f.write("VECTORS vectors float")
f.write("\n\n")
# f.close()
# for t in range (1000):
#     energy=system.analysis.energy()
#     print("--------------------- Energies --------------------------",flush=True)
#     print(energy)
#     print("Timestep is: ",system.time_step)
#     vtf.writevcf(system, outfile)
#     system.part.writevtk("./dummy.vtk");
#     write_to_vtk(vtk_idx);vtk_idx+=1
#     system.integrator.run(warm_steps)
#     count_steps+=1
#     store_all_lists()

#    #z gamma=gamma*0.999995
#     count+=1
#     steps.append(count_steps*warm_steps)
#     total_en.append(energy['total']/npart)
#     plot_energies('total',total_en,key_2)
#     kin_en.append(energy['kinetic']/(1.5*npart))
#     plot_energies('kinetic',kin_en,key_2)
#     dipolar_en.append(energy['dipolar']/npart)
#     plot_energies('dipolar',dipolar_en,key_2)
#     pair_en.append(energy['non_bonded']/npart)
#     plot_energies('non-bonded',pair_en,key_2)
#     bond_en.append(energy['bonded']/npart)
#     plot_energies('bonded',bond_en,key_2)
#     vel=np.linalg.norm(passive.v,axis=1)
#     kinetic_passive=(2/3)*(0.5*passive.mass*vel**2)
#     ken_0.append(kinetic_passive.mean())
#     plot_energies('T_0',ken_0,key_2)
#     vel=np.linalg.norm(active.v,axis=1)
#     kinetic_active=(2/3)*(0.5*active.mass*vel**2)
#     ken_1.append(kinetic_active.mean())
#     plot_energies('T_1',ken_1,key_2)

# print("======= ======== =======")
# print("WARMUP DONE!!")
# store_all_lists()

# angle_harmonic = interactions.AngleHarmonic(bend=1.0, phi0=2 * np.pi / 3)
# system.bonded_inter.add(angle_harmonic)
# nchains=int(system.part[:].pos[:,0].size/DP)
# monomers=DP
# for j in range (nchains):
#     for i in range(int(j*monomers),int((j+1)*monomers)-2):
#         id=i
#         system.part[id+1].add_bond((angle_harmonic,id,id+2))



# # restore simulation temperature
# system.thermostat.set_langevin(kT=1.0, gamma=1.0)
# system.integrator.run(warm_steps * 100)
# print("Finished warmup")


H_constraint = espressomd.constraints.HomogeneousMagneticField(H=H_field)
system.constraints.add(H_constraint)

#############################################################
#      Integration                                          #
#############################################################
log = open("myprog_magnetics_field_integrate.log", "w")
sys.stdout = log
system.integrator.set_vv() 
key_2="magnetics_equil"
print("simulating...")
check_path="mycheck_magnetic"
checkpoint=checkpointing.Checkpoint(checkpoint_id=check_path,checkpoint_path='.')
checkpoint.register("system")
t_steps = 10000 
frames=1000
store_step=int(t_steps/frames)
check_idx=0
checkpoint.save(checkpoint_index=check_idx)
warm_steps=10000

## Setup lists for data storage as well as indices
count,steps,count_steps,total_en,kin_en,dipolar_en,pair_en,bond_en,ken_0,ken_1,listoflists,string_lists=list_setup()

############################################################
tostore_vtk="vtk_equil_field_on"                            #
if os.path.isdir(tostore_vtk):                              #
    pass                                                    #
else:                                                       #
    os.mkdir(tostore_vtk)                                   #
tostore_plots='magnetic_plots_field_on_equil'               #
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

    check_idx+=1
    if t%store_step==0:
        checkpoint.save(checkpoint_index=check_idx)
store_all_lists()

outfile.close()
