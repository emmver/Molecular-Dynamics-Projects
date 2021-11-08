
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



checkpoint=checkpointing.Checkpoint(checkpoint_id=sys.argv[1],checkpoint_path='.')
last_idx=checkpoint.get_last_checkpoint_index()
print("last idx:",last_idx)
checkpoint.load(checkpoint_index=last_idx)

energy=system.analysis.energy()
print("--------------------- Energies --------------------------",flush=True)
print(energy)




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

DP=200
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


print("Bonded Interactions")
print(system.part[:].bonds)
