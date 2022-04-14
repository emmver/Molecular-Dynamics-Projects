
"""
This reads files produced from analysis scripts and plots the pair potential energy per particle and the number of bonds as a function of cycles of deformation
"""

import numpy as np
import logging
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


strains=[0.5,1,3,5,10,30,50,70,100,200,400,700,1000]
# strains_Pe01=[0.5,1,3,5,10,30,50,70,100,200,400,700,1000]
# strains_Pe1=[0.5,1,3,5,10,30,50,70,100,200,400,700,1000]
# strains_Pe10=[0.5,1,3,5,10,30,50,70,100,200,400,700,1000]
