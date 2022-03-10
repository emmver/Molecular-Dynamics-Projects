import numpy as np 
#import MDAnalysis as mda
from matplotlib import pyplot as plt 
from numba import njit
import sys
from timeit import default_timer as timer
from scipy.stats import gaussian_kde,iqr
import glob
import ctypes
import matplotlib as mpl
import os
from numba import njit
import tarfile
import glob 
import os 
import sys
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import KFold
from sklearn.neighbors import KernelDensity
import joblib
import scipy.stats as st
from PIL import Image,ImageDraw,ImageFont


import seaborn as sns 
os.chdir(sys.argv[1])

lin_path=sys.argv[2]
ring_path=sys.argv[3]
dictio=['active','passive','total']
lambdas=np.arange(4)+1
fig, axarr = plt.subplots(4,5,figsize=(16,5),dpi=300)

xidx=1
yidx=1
wspace=0.05
hspace=0.0
#### Setting array elements with text ######
for i in range(4):
    for j in range(5):
        axarr[i,j].axis('off')
        print(i,j)

ft_size=12
axarr[1,0].text(0.5,0.5,'Active',fontsize=ft_size)
axarr[2,0].text(0.5,0.5,'Passive', fontsize=ft_size)
axarr[3,0].text(0.5,0.5,'Total', fontsize=ft_size)
axarr[0,1].text(0.5,0.5,r'$\lambda = 1$', fontsize=ft_size)
axarr[0,2].text(0.5,0.5,r'$\lambda = 2$', fontsize=ft_size)
axarr[0,3].text(0.5,0.5,r'$\lambda = 3$', fontsize=ft_size)
axarr[0,4].text(0.5,0.5,r'$\lambda = 4$', fontsize=ft_size)


for key in dictio:
    files_path='/gyration_results/'+key+'/2d_distributions'
    xidx=1
    for l in lambdas: 
        image1=Image.open(ring_path+files_path+'/lambda_%d_prolat_rg_kde.jpg'%(l) )
        draw = ImageDraw.Draw(image1)
        font = ImageFont.truetype("/usr/share/fonts/truetype/ubuntu/UbuntuMono-B.ttf", 146)
        draw.text((400, 100),"Ring",'white',font=font)
        image2=Image.open(lin_path+files_path+'/lambda_%d_prolat_rg_kde.jpg'%l)
        draw = ImageDraw.Draw(image2)
        font = ImageFont.truetype("/usr/share/fonts/truetype/ubuntu/UbuntuMono-B.ttf", 146)
        draw.text((400, 100),"Linear",'white',font=font)
        image1_size=image1.size
        image2_size=image2.size
        new_image=Image.new('RGB',(2*image1_size[0], int(image1_size[1])), (250,250,250))
        new_image.paste(image1,(0,0))
        new_image.paste(image2,(image1_size[0],0))
        #axarr[yidx,xidx].text(1,1,'Ring',fontsize=ft_size/2)
        axarr[yidx,xidx].imshow(new_image)
        #axarr[yidx,xidx].text(1.5,1.05,'Linear')
        xidx+=1
        #new_image.show()
    yidx+=1

plt.subplots_adjust(wspace=wspace, hspace=hspace)

#plt.show()
fig.savefig('./all_2d_dist.jpg',dpi=300,bbox_inches='tight')