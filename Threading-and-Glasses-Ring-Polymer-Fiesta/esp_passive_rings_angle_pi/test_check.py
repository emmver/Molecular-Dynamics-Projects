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


checkpoint=checkpointing.Checkpoint(checkpoint_id='mycheck',checkpoint_path='.')
checkpoint.load()
outfile = open('check.vtf', 'w') 
vtf.writevsf(system, outfile)
