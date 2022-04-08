
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

def read_traj(f,npart):
    out=[]
    for i in range (npart):
            line = f.readline()
            elems=str.split(line.strip()," ")
            #print(elems)
            for el in elems: 
                #print(el)
                out.extend([float(el)])
                if len(elems)>0: elsize=len(elems);
    out=np.reshape(np.asarray(out),(int(len(out)/elsize),elsize))
    f.readline()
    f.readline()
    return out

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
    pos_dummy=tostore_vtk+"/pos_dummy.txt"
    dip_dummy=tostore_vtk+"/dipoles_dummy.txt"
    bound_idx_dummy=tostore_vtk+"/bidx_dummy.txt"
    vel_dummy=tostore_vtk+"/veloc_dummy.txt"
    np.savetxt(pos_dummy,system.part[:].pos)
    np.savetxt(bound_idx_dummy,(system.part[:].pos/system.box_l))
    np.savetxt(vel_dummy,system.part[:].v)
    np.savetxt(dip_dummy,system.part[:].dip)
    os.system("cat ./pos_h_file.txt " +pos_dummy+ " ./bidx_h_file.txt "+ bound_idx_dummy + " ./vel_h_file.txt "+ vel_dummy +" ./dipoles_h_file.txt " + dip_dummy +" > " + store_vtk_file )
    os.system("cat ./dummy.vtk " + store_vtk_file + " > ./all.vtk " )
    #os.sytem("rm -rf " + store_vtk_file)
def vtk_heads():
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
    for i in range (len(listoflists)):
        np.savetxt(tostore_plots+"/"+string_lists[i],listoflists[i])

''' 
For reading VTK: Total 612 lines 
Positions: np.genfromtxt(filename, skip_header=5, skip_footer=407)
Velocities: np.genfromtxt(filename, skip_header=209, skip_footer=204)
Dipoles: np.genfromtxt(filename, skip_header=413)
'''
os.chdir(sys.argv[1])
track=np.arange(9900,10000)
idx=np.random.randint(9)
skipframes=track[idx]
print('Skipframes',track[idx])
check_folder='../vtk_equil'
filename=check_folder+'/part_test_%d.vtk'%999
print(sys.argv[1])
print(filename)
L=73.36369427
f=open('../polymer_magnetic.vtf', 'r')
DP=200
nch=1
npart=nch*DP
bonds=npart
angles=0
for ii in range (int(npart)+bonds+angles+2):
        f.readline()
print(f.readline())
for ii in range (skipframes):
        data=read_traj(f,npart)

pos=read_traj(f,npart)[:,1:]
pos=pos-np.floor(pos/L)
vel=np.genfromtxt(filename, skip_header=208, skip_footer=201)
dip=np.genfromtxt(filename, skip_header=412)
print('Positions:', pos.shape)
print(pos)
print('============================')
print('Velocities:', vel.shape)
print(vel)
print('============================')
print('Dipoles:', dip.shape)
print(dip)
print('============================')


''' Adding Positions, Bonds, Angles, Velocities and Dipoles to system'''


system = espressomd.System(box_l=3*[L])
system.set_random_state_PRNG()
system.time_step=0.0005
system.min_global_cut=1.0
system.cell_system.skin=0.4
system.thermostat.set_langevin(kT=1.0,gamma=1.0,seed=42)
system.cell_system.set_domain_decomposition(use_verlet_lists=True)
#### Pair Potential #######
system.non_bonded_inter[0,0].wca.set_params(epsilon=1, sigma=1)
system.non_bonded_inter[1,1].wca.set_params(epsilon=1, sigma=1)
system.non_bonded_inter[0,1].wca.set_params(epsilon=1, sigma=1)
######## Bonded Potential ####
fene = interactions.FeneBond(k=30, d_r_max=1.5)
system.bonded_inter.add(fene)
# magnetic field times dipole moment
H_dipm = int(sys.argv[2])
dummy=np.random.randn(3)
if os.path.isfile('H_%d.field'%H_dipm)==True:
    H_field=np.genfromtxt('H_%d.field'%H_dipm)
else:
    H_field = H_dipm*(dummy/np.linalg.norm(dummy))
    np.savetxt('H_%d.field'%H_dipm,H_field)
############################################################

### Add positions and bonds#######
#### Add dipoles #######
idx=0
how_many=0.5
###pos_fin=np.zeros(pos.shape)
#pos_fin[0,:]=pos[0,:]
for ii in range(pos.shape[0]):
    print('Pos:',ii)
    # if ii>0: 
    #     dist=np.abs(pos[ii,:]-pos[ii-1,:])
    #     dist_magn=np.linalg.norm(dist)
    #     print(dist_magn)
    #     t_idx=0
    #     pos_fin[ii,:]=pos[ii,:]
    #     if dist_magn>1.5:
    #         if dist[0]<=L/2:
    #             pos_fin[ii,0]=pos[ii,0]+L
    #         else: 
    #             pos_fin[ii,0]=pos[ii,0]-L
    #     if dist[1]>1.3:
    #         if dist[1]<=L/2:
    #             pos_fin[ii,1]=pos[ii,1]+L
    #         else: 
    #            pos_fin[ii,1]=pos[ii,1]-L
    #     if dist[2]>1.3:
    #         if dist[2]<=L/2:
    #            pos_fin[ii,2]=pos[ii,2]+L
    #         else: 
    #             pos_fin[ii,2]=pos[ii,2]-L
    system.part.add(id=ii,pos=pos[ii,:])
    i=ii
    if idx<int(how_many*dip.shape[0]):
        system.part[i].mol_id=int(np.floor(i/pos.shape[0]))
        system.part[i].type=1
        system.part[i].dip=dip[i]
        system.part[i].dipm=np.linalg.norm(dip[i])
        system.part[i].rotation=[1,1,1]
        idx+=1
    else:
        system.part[i].mol_id=int(np.floor(i/pos.shape[0]))
        system.part[i].type=0
        system.part[i].rotation=[1,1,1]
p3m = magnetostatics.DipolarP3M(prefactor=1,accuracy=1.2e-3)

system.actors.add(p3m)

#p3m.tune(False)


###### Add velocities #######
for ii in range(vel.shape[0]):
    system.part[ii].v=vel[ii]

### Add bonds ####
for ii in range(pos.shape[0]):
    if ii==0:
        system.part[ii].add_bond((fene,ii+pos.shape[0]-1))
    else:
        system.part[ii].add_bond((fene,ii-1))

# for i in range (dip.shape[0]):
#     if idx<int(how_many*dip.shape[0]):
#         system.part[i].mol_id=int(np.floor(i/pos.shape[0]))
#         system.part[i].type=1
#         system.part[i].dip=dip[i]
#         system.part[i].dipm=np.linalg.norm(dip[i])
#         system.part[i].rotation=[1,1,1]
#         idx+=1
#     else:
#         system.part[i].mol_id=int(np.floor(i/pos.shape[0]))
#         system.part[i].type=0
#         system.part[i].rotation=[1,1,1]


### Add angular potential ######
DP=pos.shape[0]
npart=pos.shape[0]
angle_harmonic = interactions.AngleCosine(bend=1.5, phi0=np.pi)
system.bonded_inter.add(angle_harmonic)
system.part[0].add_bond((angle_harmonic,(0+1)*(pos.shape[0])-1,0*DP+1))
system.part[(0+1)*pos.shape[0]-2].add_bond((angle_harmonic,(0+1)*(DP)-1,(0)*(DP)))
for i in range(0,DP-2):
    id=i
    system.part[id+1].add_bond((angle_harmonic,id,id+2))
#energy=system.analysis.energy();
#print("--------------------- Energies Frame 0 --------------------------",flush=True)

#print(energy,flush=True)

tostore_vtk="vtk_warmup_field"
vtk_idx=0
if os.path.isdir(tostore_vtk):
    pass 
else:
    os.mkdir(tostore_vtk)

tostore_plots='magnetic_plots_warmup_field'
if os.path.isdir(tostore_plots):
    pass 
else:
    os.mkdir(tostore_plots)



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
print("#########FIELD#######")
print(H_field)

#### Warmup ####
active=system.part.select(type=1)
passive=system.part.select(type=0)
outfile = open('first_frame.vtf', 'w')
vtf.writevsf(system, outfile)
print("VSF OK")
vtf.writevcf(system, outfile)
print("VTF OK")


# outfile = open('polymer_dipoled_field_warm.vtf', 'w')

#f=open ("vectors_h_file.txt","w");
#f.write("\n\n")
#f.write("VECTORS vectors float")
#f.write("\n\n")
#f.close()
# key_2='magnetic_warmup_field'
# for t in range (1000):
#     energy=system.analysis.energy()
#     print("--------------------- Energies --------------------------",flush=True)
#     print(energy)
#     print("Timestep is: ",system.time_step)
#     vtf.writevcf(system, outfile)
#     system.part.writevtk("./dummy.vtk");
#     write_to_vtk(vtk_idx);vtk_idx+=1
#     system.integrator.run(1)
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

# outfile.close()



# Setup lists for storing data
print('Start Warmup')

for i in range(100):
    system.integrator.run(1)
    print("progress: {:3.0f}%, dipolar energy: {:9.2f}".format(
        (i + 1) * 100. / 100,system.analysis.energy()["dipolar"]), end="\r")
    vtf.writevcf(system, outfile)
outfile.close()
##################### Done ########################################
vtk_heads()



##################################################################
#                                                                #
#    Restart from equilibrated frame of magnetic snapshot done   #
#                                                                #
##################################################################




outfile = open('polymer_dipoled_field.vtf', 'w')

vtf.writevsf(system, outfile)
vtf.writevcf(system, outfile)


print("DONE!")



H_constraint = espressomd.constraints.HomogeneousMagneticField(H=H_field)
system.constraints.add(H_constraint)
print('Constrain test')
print(system.constraints[0])
#############################################################
#      Integration                                          #
#############################################################
log = open("myprog_magnetics_field_integrate.log", "w")
sys.stdout = log
system.integrator.set_vv() 
key_2="magnetics_field_equil"
print("simulating...")
##check_path="mycheck_field_on"
#checkpoint=checkpointing.Checkpoint(checkpoint_id=check_path,checkpoint_path='.')
#checkpoint.register("system")
t_steps = 10000 
frames=10000
store_step=int(t_steps/frames)
#check_idx=0
#checkpoint.save(checkpoint_index=check_idx)
warm_steps=10000
system.time_step=5e-3
system.cell_system.tune_skin(min_skin=0.4, max_skin=2., tol=0.2, int_steps=100)
print('tuned skin = {}'.format(system.cell_system.skin))

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

vtk=open("./dummy.vtk",'w');

for t in range(t_steps):
    energy=system.analysis.energy()
    print("step {} of {}".format(t, t_steps))
    
    print("--------------------- Energies --------------------------")
    print(energy)
    print("progress: {:3.0f}%, dipolar energy: {:9.2f}".format(
        (i + 1) * 100. / t_steps, energy["dipolar"]), end="\r")


    #if t==2:
#	system.
    system.integrator.run(warm_steps)
    vtf.writevcf(system, outfile)

    if t%1==0:
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


store_all_lists()

outfile.close()
