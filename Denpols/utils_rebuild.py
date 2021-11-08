import numpy as np

def C_pd_L(_input, elle):

    out = _input - np.floor(_input / elle) * elle
    return out


def read_lammpstrj(f):

    out = list()
    # f = open(inputfile, 'r')
    box = np.empty(3)

    trash = f.readline()
    trash = f.readline(); elems = trash.split(" "); time = float(elems[0])
    trash = f.readline()
    trash = f.readline(); elems = trash.split(" "); npart = int(elems[0])
    trash = f.readline()
    trash = f.readline(); elems = trash.split(" "); box[0] = float(elems[1]) - float(elems[0])
    trash = f.readline(); elems = trash.split(" "); box[1] = float(elems[1]) - float(elems[0])
    trash = f.readline(); elems = trash.split(" "); box[2] = float(elems[1]) - float(elems[0])
    trash = f.readline()

    if box[0] == box[1] and box[1] == box[2]:
        box = box[0]

    i = 0
    for i in range(npart):
        line = f.readline()
        elems = line.strip(" \r\n").split(" ")
        for el in elems:
            out.extend([float(el)])
            if len(elems) > 0:
                elsize = len(elems)

    out = np.asarray(out)
    out = out.reshape( ((out.size // elsize), int(elsize)) )
    nparticles = out.shape[0]

    return time, box, nparticles, out[out[:, 0].argsort()]

def mybroadcast_to(_input, _shape):

    _input = np.asarray(_input)
    _input = np.ravel(_input)
    out_flatshape = np.prod(_shape)
    out = np.empty(out_flatshape, dtype=_input.dtype)
    inshape = _input.size

    for i in range(int(out_flatshape / inshape)):
        out[i * inshape : (i + 1) * inshape] = _input

    return np.reshape(out, _shape, order="C")
    
def rotation_matrix(input_vec, input_ref):

    out = np.empty((3, 3), dtype=float)

    n1 = np.linalg.norm(input_vec)
    if n1 != 1.0:
        v1 = input_vec / n1
    else:
        v1 = input_vec

    n2 = np.linalg.norm(input_ref)
    if n2 != 1.0:
        v2 = input_ref / n2
    else:
        v2 = input_ref

    c = np.dot(v1, v2)
    if c == -1:
        out[0,0] = 1-2*v1[0]**2; out[0,1] = -2*v1[0]*v1[1]; out[0,2] = -2*v1[0]*v1[2];
        out[1,0] = out[0,1]; out[1,1] = 1-2*v1[1]**2; out[1,2] = -2*v1[1]*v1[2];
        out[2,0] = out[0,2]; out[2,1] = out[1,2]; out[2,2] = 1-2*v1[2]**2
    else:

        if np.array_equal(v1, v2) == True:
            out = np.identity(3)
        else:
            v = np.cross(v1, v2)
            vmat = np.reshape( np.asarray([0, -v[2], v[1], v[2], 0, -v[0], -v[1], v[0], 0]), (3, 3) )
            out = np.identity(3) + vmat + np.dot(vmat, vmat) * (1.0 / (1.0 + c))

    return v1, v2, out


def center_of_mass(_input, box, mytype):

    box1 = np.asarray(box)

    if mytype != "pbc" and mytype != "nopbc":
        print("only possible types: pbc and nopbc")
        exit(1)

    if mytype == "pbc":

        if box.size == 1:
            box1 = np.asarray([box, box, box])
        else:
            box1 = np.copy(box)

        t = _input*2.*np.pi/box1

        xi = np.average(np.cos(t), axis=0)
        zeta = np.average(np.sin(t), axis=0)

        out = box1*(np.arctan2(-zeta,-xi) + np.pi)*0.5/np.pi

    if mytype == "nopbc":
        out = np.average(_input, axis=0)
        
    return out



