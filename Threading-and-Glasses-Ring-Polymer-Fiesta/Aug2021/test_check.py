print("Importing....")
import espressomd
espressomd.assert_features(["WCA"])
from espressomd import thermostat
from espressomd import interactions
from espressomd import polymer
import numpy as np 
from matplotlib import pyplot as plt
from espressomd.interactions import FeneBond, AngleCosine
from espressomd.io.writer import vtf  # pylint: disable=import-error
from espressomd import checkpointing
import sys
import os
import re
print("Done Importing",flush=True)
print('Work Dir:',os.getcwd(),flush=True)
print('Work Dir Type:',type(os.getcwd()),flush=True)
path_to_use=str(sys.argv[1])
print("Type is:", type(path_to_use),flush=True)
print("Instance", isinstance(path_to_use,str),flush=True)
if not isinstance(path_to_use,int):
    print("If not works",flush=True)
print("Boolean:",bool(re.compile(r"[^a-zA-Z0-9_\-]").search(path_to_use)),flush=True)
print("Path is:",os.path.isdir("./"+path_to_use),flush=True)
checkpoint=checkpointing.Checkpoint(checkpoint_id=path_to_use,checkpoint_path=os.getcwd())
last_idx=checkpoint.get_last_checkpoint_index()
print("last idx:",last_idx)
checkpoint.load(checkpoint_index=last_idx)
outfile = open('check.vtf', 'w') 
print("Positions=\n {}".format(system.part[:].pos),flush=True)
vtf.writevsf(system, outfile)
vtf.writevcf(system,outfile)
DP=200
n_chains=150
density=0.85 
box = ((n_chains*DP)/density)**(1/3)
npart=DP*n_chains
print ("L: ", box) 
print("Particles:",npart)
print("# Chains:",int(n_chains),flush=True)
#npart=DP*nchains
hot_idx=[]
dp_count=0
for i in range (npart):
    #print(checkpoint[i,0],checkpoint[i,1], checkpoint[i,2], checkpoint[i,3])
    if dp_count<int((1/8)*DP):
            system.part[i].temp=3.0
            system.part[i].mol_id=int(np.floor(i/DP))
            system.part[i].type=2

        #hot_idx.append(i)
    else:
            system.part[i].temp=1.0
            system.part[i].mol_id=int(np.floor(i/DP))
            system.part[i].type=1
        #system.part.add(id=int(checkpoint[i,0]),mol_id=int(np.floor(i/DP)), type=1,pos=(checkpoint[i,1], checkpoint[i,2], checkpoint[i,3]),temp=1.0)
    
    dp_count+=1
    if dp_count==DP:
        dp_count=0

print("#### Temp Test ######")

print("Temperatures=\n {}".format(system.part[:].temp))
print("#### Non Bonded Test######")
print("Non-Bonded:\n {}".format(system.non_bonded_inter[0,0].wca.get_params()))
print("\n### system.thermostat test ###")
print("system.thermostat.get_state() = {}".format(system.thermostat.get_state()))

seed=np.random.randint(1e6)
index_r=0
#index_r=sys.argv[2]
filename='/binfl/lv71642/manos4/active_rings/active_%d_%d.vtf'%(index_r,seed)
import os
#while os.path.isfile(filename):
#	index_r+=1
#filename='active_%d.vtf'%(index_r)
outfile = open(filename, 'w') 
vtf.writevsf(system, outfile)
vtf.writevcf(system, outfile)
system.time_step = 0.005



#############################################################
#      Integration                                          #
#############################################################
hot=system.part.select(type=2)
cold=system.part.select(type=1)
t_steps = int(1E7)
framestep=1000
tostore=np.zeros((t_steps,6))
header="time\ttot_energy\tkin_energy_cold\tkin_energy_hot\tbox_size\tdens\t"
#seed=np.random.randint()
checkpoint=checkpointing.Checkpoint(checkpoint_id='mycheck_active_%d_%d'%(index_r,seed),checkpoint_path='/binfl/lv71642/manos4/active_rings')
checkpoint.register("system")
gamma=5
#seed=np.random.randint()
system.thermostat.set_langevin(kT=1.0,gamma=gamma,seed=seed)
check_id=0
for t in range(t_steps):
    system.integrator.run(10)
    count=0
    
    print("Step:",t)
    print("gamma:",gamma,flush=True)
    if t%framestep==0 and gamma<=1:
        print("Steps are is:",system.time/system.time_step,flush=True)
        vtf.writevcf(system, outfile)
    checkpoint.save(checkpoint_index=check_id)
    check_id+=1
        #print("Frame:", int(t/framestep))
    #en_1=system.part[[cold]].analysis.energy()
    #en_2=system.part[[hot]].analysis.energy()
    en=system.analysis.energy()
    #ken_1=en_1['kinetic']/(1.5*system.number_of_particles(1))
    #ken_2=en_1['kinetic']/(1.5*system.number_of_particles(2))
    #ken=en['kinetic']/(1.5*npart)
    box=system.box_l[0]
    dens=npart/box**3
    tostore[t,0]=system.time
    tostore[t,1]=en['total']
    tostore[t,2]=np.mean(cold.temp)
    tostore[t,3]=np.mean(hot.temp)
    tostore[t,4]=box
    tostore[t,5]=dens
    if (t%1e5)==0:
        np.savetxt("log_active_%d_%d.out"%(index_r,seed),tostore[:t],delimiter='\t',header=header)
    if (t%1000)==0 and gamma>1 :
        gamma=gamma-0.05*gamma
        gamma=int(gamma)
        system.thermostat.set_langevin(kT=1.0,gamma=gamma)

np.savetxt("log_active_%d_%d.out"%(index_r,seed),tostore[:t],delimiter='\t',header=header)

checkpoint.save()
outfile.close()
