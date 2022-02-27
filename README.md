# Molecular Dynamics Projects

Here are scripts and function libraries I'm using to run and analyze Coarse-grained molecular dynamics simulations using packages such as LAMMPS,ESSPResSo and/or HOOMD-Blue.
The analysis is done in a post-processing manner using scripts that read and then analyze the trajectories dumped by the two packages. 

There are three different projects: 

1) Colloidal Gels - Monodisperse attractive colloidal gels, this project is mainly run on the LAMMPS package, with some scripts provided. The analysis is done fully in python, using several different scripts. 

2) Denpols - Coarse-graining MD simulations of dendronized polymers, some scripts for running with HOOMD-Blue and analysis packages are here. 

3) Threading and Glasses Ring Polymers Fiesta - This is a big project containing two sub-projects. One has to do with ring polymers that carry magnetic dipoles and studying their conformation and dynamics under the influence or absence of an external magnetic field. The second project has to do with a binary blends (i.e. two different molecular weights) of ring polymers and how this blending will eventually affect the flow properties (i.e. viscosity) of the blend.  The project contains python scripts for analysis and running with ESSPResSo (for the magnetic polymers project) as well as scripts for running with LAMMPS (for the binary blends project).

This repository is for personal use, cataloguing, syncing and in general ease of access to my code, however, it is made public in the hopes, and the slight chance, it might somehow prove useful to someone. 
