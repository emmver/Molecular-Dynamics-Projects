import numpy as np
import sys
import string

from utils import read_lammpstrj, mybroadcast_to, C_pd_L, rotation_matrix, center_of_mass

def print_lammpstrj(myinput, thisL, myname):

        f = open(myname,"w")
        howmany = myinput.shape[0]

        f.write("ITEM: TIMESTEP\n%ld\n" % 0);
        f.write("ITEM: NUMBER OF ATOMS\n%ld\n" % howmany);
        f.write("ITEM: BOX BOUNDS pp pp pp\n%.0f %.1f\n%.0f %.1f\n%.0f %.1f\n" % (0.,thisL[0],0.,thisL[1],0.,thisL[2]));
        f.write("ITEM: ATOMS id type x y z\n");
        for i in range(howmany):
                f.write("%d %d %.5f %.5f %.5f 0 0 0\n" % (i+1,0, myinput[i,0], myinput[i,1], myinput[i,2]));

        f.close()


def make_dendron(gen):

	outmap = [[0,1],[0,2]]
	idx = [0, [1,2]]
    
	ww=3
	thetas = [20,20,20,20,20,20,20]
    
	for i in range (2,gen+1):
		blobs=2**(gen+1)-1
        	dendron=np.zeros((blobs,3))
        	a=np.zeros((1,3))
        	theta0 = 60
        	b=np.array([1*np.cos(theta0*np.pi/180),1*np.sin(theta0*np.pi/180),0])+np.array([0,0.1,0.1])
        	c=np.array([1*np.cos((180-theta0)*np.pi/180),1*np.sin((180-theta0)*np.pi/180),0])+np.array([0,0.1,-0.1])
       		dendron_1 = [a, [b,c]]
        	dendron_2 = [a, [b,c]]
        
        	outmap = [[0,1],[0,2]]
        	idx = [0, [1,2]]
    
        	ww=3
        	thetas = [20,20,20,20,20,20,20]
        	for i in range (2,gen+1):
			dendron_1.append([])
			dendron_2.append([])
			idx.append([])
			curr_g = dendron_2[i-1]
			curridx = idx[i-1]

			rot=np.array([np.cos(thetas[i]*np.pi/180),np.sin(thetas[i]*np.pi/180),0,-np.sin(thetas[i]*np.pi/180),np.cos(thetas[i]*np.pi/180),0,0,0,1])
			rot=rot.reshape(3,3)
			rot1=np.copy(rot) 
			rot1[0,1] = -rot1[0,1]; rot1[1,0] = -rot1[1,0] #Changing first and second row?
        
        		for j in range(2**(i-1)):
                		v = dendron_1[i-1][j]
                		ref_idx = curridx[j]
               			for k in [0,1]:
                        
                    			if k == 0:
                        			dendron_1[i].extend([np.dot(rot,v)])
                        	    		dendron_2[i].extend([np.dot(rot,v)+curr_g[j]])
                            			idx[i].extend([ww]); outmap.append([ref_idx,ww]); ww+=1
                       			else:	
						dendron_1[i].extend([np.dot(rot1,v)])
                   				dendron_2[i].extend([np.dot(rot1,v)+curr_g[j]])
						idx[i].extend([ww]); outmap.append([ref_idx,ww]); ww+=1
                            
		kk=1
        	for i in range(1,gen+1):
        		kstart=kk
            		for j in dendron_2[i]:
            			dendron[kk,:] = j
                		kk+=1

	return dendron, outmap

def assemble_unit(monomer, side_mon, dendron): #assembles monomer linker and dendron
    
	out = np.empty((dendron.shape[0]+1+1,3)) #dendron.shape[0] gives number of rows 
    
	out[0,:] = monomer 
	out[1,:] = side_mon
	_dir = side_mon - monomer; _dir /= np.linalg.norm(_dir)  
	CdM = center_of_mass(dendron[1:],(1,1,1), 'nopbc')
	_dir2 = CdM - dendron[0]; _dir2 /= np.linalg.norm(_dir2)
	nn, nn, mat = rotation_matrix(_dir2, _dir)
	oldd = np.copy(dendron); 
	for i in range(oldd.shape[0]):
		dendron[i] = np.dot(mat,oldd[i])
    	out[2:,:] = dendron + mybroadcast_to(monomer + 2*_dir, (out.shape[0]-2,3)) 
    
	return out


if len(sys.argv) != 5:
	print "usage is python"+sys.argv[0]+" infile oldNbackbone newgen outfile"
	exit(1)

myfile = sys.argv[1]
Nmon = int(sys.argv[2])
newgen = int(sys.argv[3])
outfold = sys.argv[4]

newN = Nmon*(2+(2**(newgen+1)-1))
ldendron = 2**(newgen+1)-1 #length of dendron array
nbonds_dendr = np.sum(np.asarray([2**i for i in range(1,newgen+1)]))

fz = open(myfile,"r")
_time, mybox, nparticles, data = read_lammpstrj(fz)
fz.close()

CoM = center_of_mass(data[:,2:5], mybox, 'pbc')
dd = mybox[:]/2. - CoM
temp = C_pd_L(data[:,2:5] + mybroadcast_to(dd,(data.shape[0],3)), mybox) 
data[:,2:5] = np.copy(temp)
##mybox = [mybox[0], mybox[0],mybox[0]]

backbone = data[np.where(data[:,1] == 1)]
sides = data[np.where(data[:,1] == 2)]

allparticles = np.empty((newN,3),dtype = float)

f = open(outfold+"/data_dendron.dat", "w+") #opens file and writes the following f.write

f.write("LAMMPS data file\n\n")

f.write("%ld atoms\n" % (Nmon*(2+(2**(newgen+1)-1)))) #number of atoms
f.write("%ld bonds\n" % (Nmon-1+Nmon*2+Nmon*nbonds_dendr) ) #number of bonds
f.write("%ld angles\n" % (Nmon + Nmon-2) )#number of angles
f.write("0 dihedrals\n0 impropers\n\n")

f.write("3 atom types\n3 bond types\n2 angle types\n0 dihedral types\n0 improper types\n\n")

f.write("0.0 %.1f xlo xhi\n0.0 %.1f ylo yhi\n0.0 %.1f zlo zhi\n\n" % (mybox[0], mybox[1], mybox[2]))

f.write("Masses\n\n1 1\n2 1\n3 1\n\n") #masses of particle types

f.write("Atoms\n\n")

dendron, bondmap= make_dendron(newgen) #makes dendron coordinates and mapping

kk = 1
radius0 = 1.1
for i in range(Nmon):
	monomer = backbone[i,2:5]; side_m = sides[i,2:5]
        dendronized = assemble_unit(monomer, side_m, dendron) #calls function with monomer coordinates and dendron coordinates 
        #uses allparticles array's corrseponding row to copy the dendronized array 
        allparticles[i*(2+(2**(newgen+1)-1)):(i+1)*(2+(2**(newgen+1)-1)),:] = np.copy(dendronized) 
        #writes coordinates of atoms of type 1 changing line with loop, changing number with kk
        f.write("%ld 1 1 %.5f %.5f %.5f\n" % (kk, dendronized[0,0],dendronized[0,1],dendronized[0,2])); kk+=1
        f.write("%ld 1 2 %.5f %.5f %.5f\n" % (kk, dendronized[1,0],dendronized[1,1],dendronized[1,2])); kk+=1
        for j in range((2**(newgen+1)-1)):
                f.write("%ld 1 3 %.5f %.5f %.5f\n" % (kk, dendronized[j+2,0],dendronized[j+2,1],dendronized[j+2,2]))
                kk+=1

f.write("\n")

print_lammpstrj(allparticles, mybox, outfold+"/initial_condition.lammpstrj")

f.write("Velocities\n\n")

for i in range(Nmon*(2+(2**(newgen+1)-1))):
    temp = np.random.randn(3)
    f.write("%ld %.3f %.3f %.3f\n" % (i+1, temp[0], temp[1], temp[2]))
    
    
f.write("\n")

f.write("Bonds\n\n")

kk = 1
for k in range(Nmon-1):    
    f.write("%ld 1 %ld %ld\n" % (kk, k*(2+(2**(newgen+1)-1))+1, (k+1)*(2+(2**(newgen+1)-1))+1))
    kk+=1
    
for k in range(Nmon):    
    f.write("%ld 2 %ld %ld\n" % (kk, k*(2+(2**(newgen+1)-1))+1, k*(2+(2**(newgen+1)-1))+2))
    kk+=1
    f.write("%ld 2 %ld %ld\n" % (kk, k*(2+(2**(newgen+1)-1))+2, k*(2+(2**(newgen+1)-1))+3))
    kk+=1

##print bondmap
for k in range(Nmon):
    for w, l in bondmap:
        ##print w,l 
        f.write("%ld 3 %ld %ld\n" % (kk, k*(2+(2**(newgen+1)-1))+3+w, k*(2+(2**(newgen+1)-1))+3+l))
        
f.write("\n")

f.write("Angles\n\n")

kk=1

for k in range(Nmon-2):
    f.write("%ld 1 %ld %ld %ld\n" % ( kk, k*(2+(2**(newgen+1)-1))+1, (k+1)*(2+(2**(newgen+1)-1))+1, (k+2)*(2+(2**(newgen+1)-1))+1))
    kk+=1        
                       
for k in range(Nmon):
    f.write("%ld 2 %ld %ld %ld\n" % ( kk, k*(2+(2**(newgen+1)-1))+1, k*(2+(2**(newgen+1)-1))+2, k*(2+(2**(newgen+1)-1))+3))
    kk+=1

f.close()

