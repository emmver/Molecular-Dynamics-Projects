import warnings
import os
import sys
import freud
import matplotlib.cm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib as mpl
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import normalize
from sklearn.svm import SVC
from umap import UMAP
from unitcell import  UnitCell as unitcell
warnings.filterwarnings("ignore")

def get_features(box, positions, structure):
    voro = freud.locality.Voronoi()
    voro.compute(system=(box, positions))
    nlist = voro.nlist.copy()
    nlist.filter(nlist.weights > 0.1)
    features = {}
    for l in [4, 6, 8, 10, 12]:
        ql = freud.order.Steinhardt(l=l, weighted=True)
        ql.compute(system=(box, positions), neighbors=nlist)
        features[f"q{l}"] = ql.particle_order

    return features



def write_to_lammps(filename,box,n,positions):
    f=open(filename,'w')
    f.write('ITEM: TIMESTEP\n')
    f.write('0\n')
    f.write('ITEM: NUMBER OF ATOMS\n')
    f.write('%d\n'%n)
    f.write('ITEM: BOX BOUNDS xx yy zz pp pp pp\n')
    f.write('%.2f %.2f\n'%(-box,box))
    f.write('%.2f %.2f\n'%(-box,box))
    f.write('%.2f %.2f\n'%(-box,box))
    f.write('ITEM: ATOMS id x y z\n')
    for i in range(positions[:,0].size):
        f.write('%d %.3f %.3f %.3f\n'%(i,positions[i,0],positions[i,1],positions[i,2]))

    f.close()
    
############################# Plotting style #############################################################
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

########################################################################################################
################### Generate structure and calculate local order parameters ############################

os.chdir(sys.argv[1])
plot_dir='training_plots'
if os.path.isdir(plot_dir)==False:
    os.mkdir(plot_dir)



noises = np.linspace(0,0.1,10)
scores=np.zeros(noises.shape)
for idx,noise in enumerate(noises):
    print('For Noise:',noise)
    N = 32000
    structures = {}
    n = round((N / 2) ** (1 / 3))
    structures["bcc"] = unitcell.bcc().generate_system(n, sigma_noise=noise)
    n = round((N / 4) ** (1 / 3))
    structures["fcc"] = unitcell.fcc().generate_system(n, sigma_noise=noise)
    n = round((N / 1) ** (1 / 3))
    structures["sc"] = unitcell.sc().generate_system(n, sigma_noise=noise)
    n=round((N/1)**(1/3))
    a1=[1,-3**0.5  , 0]
    a2=[0.5,+3**0.5  , 0]
    a3=[0,0, (8/3)**0.5]
    n = round((N / 4) ** (1 / 3))
    structures["hcp"]= unitcell.hcp().generate_system(n, sigma_noise=noise)
    for name, (box, positions) in structures.items():
        print(name, "has", len(positions), "particles.")
       # print(positions[1,0])
        if noise<0.0001:
            print('Print Noise',noise)
            write_to_lammps(name,n/2,len(positions),positions)
   # sys.exit()


    structure_features = {}
    for name, (box, positions) in structures.items():
        structure_features[name] = get_features(box, positions, name)


    for l in [4, 6]:
        plt.figure(figsize=(16, 9), dpi=300)
        for name in structures.keys():
            plt.hist(
                structure_features[name][f"q{l}"],
                range=(0, 1),
                bins=100,
                label=name,
                alpha=0.7,
            )
        plt.title(fr"$q_{{{l}}}$")
        plt.legend()
        for lh in plt.legend().legendHandles:
            lh.set_alpha(1)
        plt.savefig(plot_dir+'/q_%d_noise_%.3f.jpg'%(l,noise),dpi=300,bbox_inches="tight")
        #plt.show()
        plt.clf()

    ############## Training of Support Vector Machine (SVM) ##############

    structure_dfs = {}
    for i, structure in enumerate(structure_features):
        df = pd.DataFrame.from_dict(structure_features[structure])
        df["class"] = i
        structure_dfs[structure] = df    

    df = pd.concat(structure_dfs.values()).reset_index(drop=True)

    X = df.drop("class", axis=1).values
    X = normalize(X)
    y = df["class"].values
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.33, random_state=42
    )

    svm = SVC()
    svm.fit(X_train, y_train)
    print("Score:", svm.score(X_test, y_test))

    from umap import UMAP

    umap = UMAP(random_state=42)

    X_reduced = umap.fit_transform(X)

    plt.figure(figsize=(10, 10), dpi=300)
    for i in range(max(y) + 1):
        indices = np.where(y == i)[0]
        plt.scatter(
            X_reduced[indices, 0],
            X_reduced[indices, 1],
            color=matplotlib.cm.tab10(i),
            s=8,
            alpha=0.2,
            label=list(structure_features.keys())[i],
        )
    plt.legend()
    for lh in plt.legend().legendHandles:
        lh.set_alpha(1)
    plt.savefig(plot_dir+'/struc_pred_noise_%.3f.jpg'%(noise),dpi=300,bbox_inches="tight")
    plt.clf()
    #plt.show()
    scores[idx]=svm.score(X_test, y_test)

plt.plot(noises,scores,'ko-')
plt.xlabel('Noise')
plt.ylabel('SVM Score')
plt.savefig(plot_dir+'/struc_pred_noise_%.3f.jpg'%(noise),dpi=300,bbox_inches="tight")

