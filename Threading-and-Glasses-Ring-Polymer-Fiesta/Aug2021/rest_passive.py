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
check_path=sys.argv[1]

checkpoint=checkpointing.Checkpoint(checkpoint_id=check_path,checkpoint_path='.')
print("Checking...")
print(checkpoint.has_checkpoints())
last_idx=checkpoint.get_last_checkpoint_index()
checkpoint.load(checkpoint_index=last_idx)
outfile = open('check.vtf', 'w') 
print("Positions=\n {}".format(system.part[:].pos))
vtf.writevsf(system, outfile)
vtf.writevcf(system,outfile)
DP=200
n_chains=150
density=0.85 
box = ((n_chains*DP)/density)**(1/3)
npart=DP*n_chains
print ("L: ", box) 
print("Particles:",npart)
print("# Chains:",int(n_chains))
#npart=DP*nchains


print("#### Temp Test ######")

print("Temperatures=\n {}".format(system.part[:].temp))
print("#### Non Bonded Test######")
print("Non-Bonded:\n {}".format(system.non_bonded_inter[0,0].wca.get_params()))
print("\n### system.thermostat test ###")
print("system.thermostat.get_state() = {}".format(system.thermostat.get_state()))

from datetime import datetime
date = datetime. now(). strftime("%Y_%m_%d-%I:%M:%S_%p")
index_chk=0
check_path='mycheck_passive_%s'%(date)

filename='passive_%s.vtf'%(date)
#import os
#while os.path.isfile(filename):
#    index_r+=1
#filename='active_%d.vtf'%(index_r)
outfile = open(filename, 'w') 
vtf.writevsf(system, outfile)
vtf.writevcf(system, outfile)
system.time_step = 0.005



#############################################################
#      Integration                                          #
#############################################################

t_steps = int(1E6)
framestep=1000
tostore=np.zeros((t_steps,5))
header="time\ttot_energy\tkin_energy\ttbox_size\tdens\t"
checkpoint=checkpointing.Checkpoint(checkpoint_id='mycheck_passive_%s'%(date),checkpoint_path='.')
checkpoint.register("system")
gamma=5
system.thermostat.set_langevin(kT=1.0,gamma=gamma)
check_id=0
for t in range(t_steps):
    system.integrator.run(10)
    count=0
    
    #print("Step:",t)
    if t%framestep==0 :
        print("Steps are is:",system.time/system.time_step)
        vtf.writevcf(system, outfile)
    checkpoint.save(checkpoint_index=check_id)
    check_id+=1
        #print("Frame:", int(t/framestep))
    en=system.analysis.energy()
    ken=en['kinetic']/(1.5*npart)
    box=system.box_l[0]
    dens=npart/box**3
    tostore[t,0]=system.time
    tostore[t,1]=en['total']
    tostore[t,2]=ken
    tostore[t,3]=box
    tostore[t,4]=dens
    if (t%1e5)==0:
        np.savetxt("log_passive_%s.out"%(date),tostore,delimiter='\t',header=header)
    if (t%1000)==0 and gamma>1 :
        gamma=gamma-0.05*gamma
        gamma=int(gamma)
        system.thermostat.set_langevin(kT=1.0,gamma=gamma)

np.savetxt("log_passive_%s.out"%(date),tostore,delimiter='\t',header=header)

checkpoint.save()
outfile.close()
