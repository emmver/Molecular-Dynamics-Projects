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

print("#### Temp Test ######")

print("Temperatures=\n {}".format(system.part[:].temp))
print("#### Non Bonded Test######")
print("Non-Bonded:\n {}".format(system.non_bonded_inter[0,0].wca.get_params()))
print("\n### system.thermostat test ###")
print("system.thermostat.get_state() = {}".format(system.thermostat.get_state()))
#gamma=system.thermostat.get_state(gamma)
#print(gamma)
print((system.thermostat.get_state()[0][1]))

#print("Gamma is:", system.thermostat.get_gamma())
