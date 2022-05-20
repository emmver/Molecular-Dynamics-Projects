import numpy as np
import string

def C_pd_L(input, elle):

        out=input-np.floor(input/elle)*elle
        return out

def read_lammpstrj(f):

    out = list()
    box = np.empty(3)

    trash = f.readline(); trash = f.readline();
    elems = string.split(trash," ");  time = float(elems[0]);
    trash = f.readline(); trash = f.readline();
    elems = string.split(trash," ");  npart = long(elems[0]);
    trash = f.readline();
    trash = f.readline();
    elems = string.split(trash," ");  box[0] = float(elems[1])-float(elems[0])
    trash = f.readline();
    elems = string.split(trash," ");  box[1] = float(elems[1])-float(elems[0]);
    trash = f.readline();
    elems = string.split(trash," ");  box[2] = float(elems[1])-float(elems[0]);
    trash = f.readline();

    i = 0
    for i in range(npart):

        line  = f.readline()
        elems = string.split(line.strip(" \n")," ")
        for el in elems:
                out.extend([float(el)])
                if len(elems) > 0 : elsize = len(elems)



    #f.close()
    out = np.reshape(np.asarray(out), (len(out)/elsize,elsize))
    nparticles = out.shape[0]

    return time, box, nparticles, out[out[:,0].argsort()]


def mybroadcast_to(_input, _shape):

        out_flatshape = np.prod(_shape)
        out = np.empty(out_flatshape,dtype = _input.dtype)
        inshape = _input.size

        for i in range(out_flatshape/inshape):
                out[i*inshape:(i+1)*inshape] = _input

        return np.reshape(out, _shape, order='C')



def rotation_matrix(input_vec, input_ref):

        out = np.empty((3,3),dtype = float)

        n1 = np.linalg.norm(input_vec)
        if n1 != 1.:
                v1 = input_vec/n1
        else:
                v1 = input_vec

        n2 = np.linalg.norm(input_ref)
        if n2 != 1.:
                v2 = input_ref/n2
        else:
                v2 = input_ref

        c = np.dot(v1,v2);
        if c == -1:
                print "vectors are opposite: rotation not possible!"
                exit(1)

        if np.array_equal(v1,v2) == True:
                out = np.identity(3)
        else:
                v = np.cross(v1, v2)
                #s = np.linalg.norm(v)

                vmat = np.reshape(np.asarray([0,-v[2], v[1], v[2], 0, -v[0], -v[1], v[0], 0]), (3,3))

                out = np.identity(3) + vmat + np.dot(vmat,vmat)*(1./(1.+c))

        return v1, v2, out


def center_of_mass(myinput, box, mytype):

        outx = 0; outy = 0; outz = 0

        if mytype != 'pbc' and mytype != 'nopbc':
                print "only possible types: pbc and nopbc"
                exit(1)

        if mytype == 'pbc':
                xi_x = 0; xi_y = 0; xi_z = 0;
                zeta_x = 0; zeta_y = 0; zeta_z = 0;
                thisN = myinput.shape[0]

                for i in range(thisN):
                        #outx = outx + input[i,0]; outy = outy + input[i,1]; outz = outz + input[i,2]; 
                        tx = myinput[i,0]*2.*np.pi/box[0]; xi_x = xi_x + np.cos(tx); zeta_x = zeta_x + np.sin(tx);
                        ty = myinput[i,1]*2.*np.pi/box[1]; xi_y = xi_y + np.cos(ty); zeta_y = zeta_y + np.sin(ty);
                        tz = myinput[i,2]*2.*np.pi/box[2]; xi_z = xi_z + np.cos(tz); zeta_z = zeta_z + np.sin(tz);

                xi_x = xi_x/float(thisN); xi_y = xi_y/float(thisN); xi_z = xi_z/float(thisN);
                zeta_x = zeta_x/float(thisN); zeta_y = zeta_y/float(thisN); zeta_z = zeta_z/float(thisN);

                outx = box[0]*(np.arctan2(-zeta_x,-xi_x)+np.pi)*0.5/np.pi;
                outy = box[1]*(np.arctan2(-zeta_y,-xi_y)+np.pi)*0.5/np.pi;
                outz = box[2]*(np.arctan2(-zeta_z,-xi_z)+np.pi)*0.5/np.pi;

        if mytype == 'nopbc':
                outx = np.sum(myinput[:,0]); outy = np.sum(myinput[:,1]); outz = np.sum(myinput[:,2]);

                outx = outx/myinput.shape[0]; outy = outy/myinput.shape[0]; outz = outz/myinput.shape[0];

        return np.reshape(np.asarray([outx, outy, outz]),(3))




