import numpy as np 
from matplotlib import pyplot as plt

mylist=[10,16,25,36,63,100,160,1000]
gen=int(input("Gen?"))
for i in mylist:
    print(i)
    filename='avg_bond_dist_g{}{}.txt'.format(int(gen),int(i))
    data=np.genfromtxt(filename,skip_header=1)


    plt.plot(data[1:,0],data[1:,1],linewidth="3",label=r"G=%.1d  N$_{bb}$=%.1d"%(gen,i))
    plt.xticks(fontsize='20')
    plt.yticks(fontsize='20')
    plt.locator_params(axis='x',nbins=7)
    plt.xlabel("cosθ",fontsize='35')
    plt.ylabel("P(cosθ)",fontsize='35')
    plt.yscale('log')
    plt.xlim(0.5,1)
    plt.ylim(1e-4,100)
    #plt.yscale('log')
    #plt.xscale('log')
    plt.legend(loc=1,bbox_to_anchor=(1.65,1.0),ncol=2,fontsize="medium",frameon=False)
plt.savefig("G{}_bonddist.png".format(gen),bbox_inches="tight",dpi=300)
plt.clf()
for i in mylist:
    print(i)
    filename='avg_bond_dist_g{}{}.txt'.format(int(gen),int(i))
    data=np.genfromtxt(filename,skip_header=1)


    plt.plot(data[1:,0],data[1:,1],linewidth="3",label=r"G=%.1d  N$_{bb}$=%.1d"%(gen,i))
    plt.xticks(fontsize='20')
    plt.yticks(fontsize='20')
    plt.locator_params(axis='x',nbins=7)
    plt.xlabel("cosθ",fontsize='35')
    plt.ylabel("P(cosθ)",fontsize='35')
    #plt.yscale('log')
    plt.xlim(0.5,1)
    plt.ylim(1e-4,35)
    #plt.yscale('log')
    #plt.xscale('log')
    plt.legend(loc=1,bbox_to_anchor=(1.65,1.0),ncol=2,fontsize="medium",frameon=False)
plt.savefig("G{}_bonddist_lin.png".format(gen),bbox_inches="tight",dpi=300)
