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


##################### Parsing LAMMPS Data ####################
count = len(open("indata.lammps").readlines(  ))
npart=30000
skip_h=16
skip_foot=count-npart-skip_h-4
data=np.genfromtxt("indata.lammps", skip_header=16,skip_footer=skip_foot)
dens=0.85
box=(npart/dens)**(1/3)
L=box/2
data=np.genfromtxt("indata.lammps", skip_header=16,skip_footer=skip_foot)
data[:,3:6]+=box/2
system = espressomd.System(box_l=[box,box,box])
fene = FeneBond(k=30, d_r_max=1.5,r_0=1.0)
system.bonded_inter.add(fene)
system.set_random_state_PRNG()
np.random.seed(seed=system.seed)
system.time_step = 0.01
system.cell_system.skin = 0.4
system.thermostat.set_langevin(kT=1.0, gamma=1.0, seed=42)
#system.cell_system.set_n_square(use_verlet_lists=False)

n_chains=150
DP=200
##### Adding positions from LAMMPS--Shifted box origin to [0,0,0]##################

for i in range (npart):
    system.part.add(id=i,type=int(np.floor(i/DP)), pos=(data[i,3], data[i,4], data[i,5]))
outfile = open('polymer_warmup.vtf', 'w') 

##### Setting up Bonds############
system.part[0*DP].add_bond((fene,(0+1)*(DP)-1))
for i in range(DP-1):
    system.part[i].add_bond((fene,i+1))
    
for i in range (1,n_chains):
    system.part[i*DP].add_bond((fene,((i+1)*DP-1)))
    for j in range (0,DP-1):
        system.part[j+i*(DP)].add_bond((fene,(j+1)+i*(DP)))
        
##### Setting up Angles############
angle_harmonic = AngleCosine(bend=1.5, phi0=0)
system.bonded_inter.add(angle_harmonic)
#f.write("ADDING CLOSING ANGLE:%d||%d||%d\n"%(0*DP,(0+1)*(DP)-1,(0)*(DP)+1))
system.part[0*DP].add_bond((angle_harmonic,(0+1)*(DP)-1,0*DP+1))
#f.write("ADDING SECOND TO CLOSING ANGLE:%d||%d||%d\n"%((0+1)*(DP)-1,(0+1)*DP-2,(0)*(DP)))
system.part[(0+1)*DP-2].add_bond((angle_harmonic,(0+1)*(DP)-1,(0)*(DP)))
for i in range(0,DP-2):
    #f.write("ADDING ANGLE:%d||%d||%d\n"%(i+1,i,i+2))
    system.part[i+1].add_bond((angle_harmonic,i,i+2))
#f.write('#######################33First Chain Done!#####################\n')
for i in range (1,n_chains):
    #f.write("ADDING CLOSING ANGLE:%d||%d||%d\n"%((i+1)*(DP)-1,i*DP,(i)*(DP)+1))
    system.part[i*DP].add_bond((angle_harmonic,(i+1)*(DP)-1,(i)*(DP)+1))
    #f.write("ADDING SECOND TO CLOSING ANGLE:%d||%d||%d\n"%((i+1)*(DP)-1,(i+1)*DP-2,(i)*(DP)))
    system.part[(i+1)*DP-2].add_bond((angle_harmonic,(i+1)*(DP)-1,(i)*(DP)))

    for j in range (0,DP-2):
        system.part[(j+1)+i*(DP)].add_bond((angle_harmonic,j+i*(DP),(j+2)+i*(DP)))
        
    #f.write('#####################Done with Chain %d################\n'%(i))

vtf.writevsf(system, outfile)
#### WARM-UP ######
print('\n')
print("WARM-UP")
LJ_EPS = 1.0
LJ_SIG = 1.0
LJ_CUT = 2**(1/6)* LJ_SIG
LJ_CAP = 0.1
system.non_bonded_inter[0, 0].wca.set_params(
    epsilon=LJ_EPS, sigma=LJ_SIG)
system.force_cap = LJ_CAP
WARM_STEPS = 10
WARM_N_TIME = 1
MIN_DIST = 0.95
# warmup with zero temperature to remove overlaps
system.thermostat.set_langevin(kT=0.0, gamma=1.0)
i = 0
act_min_dist = system.analysis.min_dist()
while i<WARM_N_TIME:
    vtf.writevcf(system, outfile)
    system.integrator.run(WARM_STEPS)
    system.part[:].v = [0, 0, 0]
    i += 1
    LJ_CAP = LJ_CAP*1.01
    system.force_cap = LJ_CAP

system.force_cap = 0
system.integrator.run(WARM_STEPS * 10)
# restore simulation temperature
system.thermostat.set_langevin(kT=1.0, gamma=1.0)
system.integrator.run(10000)
#system.integrator.run(WARM_STEPS * 10)
outfile.close()
#############################################################
#      Integration                                          #
#############################################################
from datetime import datetime
date = datetime. now(). strftime("%Y_%m_%d")
print(date)
index_chk=0
check_path='mycheck_passive_%s'%(date)
print("Type:",type(check_path))
print(check_path)
#while os.path.isdir(check_path)==True:
#    index_chk+=1;print("check idx:",index_chk)
#    check_path='mycheck_active_%d'%(index_chk)

outfile = open('polymer_equil_%s.vtf'%(date), 'w') 
vtf.writevsf(system, outfile)        
t_steps = int(1E7)
framestep=1000
tostore=np.zeros((t_steps,5))
header="time\ttot_energy\tkin_energy\tbox_size\tdens"
checkpoint=checkpointing.Checkpoint(checkpoint_id=check_path,checkpoint_path='.')
checkpoint.register("system")
check_idx=0
for t in range(t_steps):
    system.integrator.run(WARM_STEPS)
    if t%framestep==0:
        vtf.writevcf(system, outfile)
        np.savetxt("log.out",tostore[:t],delimiter='\t',header=header)
        checkpoint.save(checkpoint_index=check_idx)
        check_idx+=1
    en=system.analysis.energy()
    ken=en['kinetic']/(1.5*npart)
    box=system.box_l[0]
    dens=npart/box**3
    tostore[t,0]=system.time
    tostore[t,1]=en['total']
    tostore[t,2]=ken
    tostore[t,3]=box
    tostore[t,4]=dens
np.savetxt("log_%s.out"%(date),tostore[:t],delimiter='\t',header=header)
#np.savetxt("log_noh.out",tostore[:t],delimiter='\t')
checkpoint.save()
outfile.close()
