%matplotlib inline
import numpy as np # Library for multi-dimensional arrays and math operations on them 
import pandas as pd # Library for data manipulation and analysis
import os # Interfacing python with the operating system
import matplotlib.pyplot as plt #Plotting library
import seaborn as sns #Plotting library with some statistical tools 
import nltk # Natural Language Processing toolkit, similar as textblob below, both have useful functions
from textblob import TextBlob
import string
import time
import re # Basically a string processing and manipulation library
from datetime import datetime ## Manipulating datetime inputs
from matplotlib.ticker import FuncFormatter
import glob
import matplotlib.pylab as pl
### My plotting style is inputted here #####
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 150
plt.rcParams["font.family"] = "Ubuntu"
plt.style.use('C:/Users/Toumba/Documents/plotstyle.mplstyle')
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
#############################################