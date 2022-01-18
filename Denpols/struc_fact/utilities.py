import numpy as np


def read_lammpstrj(f):

    out = list()
    #f = open(inputfile, 'r')
    box = np.empty(3)

    trash = f.readline(); trash = f.readline();
    elems = trash.strip("\n").split(" ");  time = float(elems[0]);
    trash = f.readline(); trash = f.readline();
    elems = trash.strip("\n").split(" "); npart = int(elems[0]);
    trash = f.readline(); 
    trash = f.readline(); 
    elems = trash.strip("\n").split(" ");  box[0] = float(elems[1])-float(elems[0])
    trash = f.readline();
    elems = trash.strip("\n").split(" ");  box[1] = float(elems[1])-float(elems[0])
    trash = f.readline();
    elems = trash.strip("\n").split(" ");  box[2] = float(elems[1])-float(elems[0]);
    trash = f.readline();

    i = 0
    for i in range(npart):

        line  = f.readline();
        elems = line.strip(" \r\n").split(" ");
        for el in elems:
                out.extend([float(el)])
                if len(elems) > 0 : elsize = len(elems)

    #f.close()
    out = np.reshape(np.asarray(out), (int(len(out)/elsize),elsize))
    nparticles = out.shape[0]

    return time, box, nparticles, out[out[:,0].argsort()]


def print_lorenzo(myinput, thisl, myname):

    f = open(myname, "w")
    howmany = myinput.shape[0]

    f.write(".Box: %.5f , %.5f , %.5f\n" % (thisl, thisl, thisl))
    color = "magenta"

    for i in range(howmany):
        f.write(
            "%.5f %.5f %.5f @ 0.500 c[%s]\n"
            % (myinput[i, 0], myinput[i, 1], myinput[i, 2], color)
        )

    f.close()


def random_on_sphere():

    ransq = 100

    while ransq >= 1.0:
        ran_1 = 1.0 - 2.0 * np.random.rand()
        ran_2 = 1.0 - 2.0 * np.random.rand()
        ransq = np.power(ran_1, 2) + np.power(ran_2, 2)

    ranh = 2.0 * np.sqrt(1.0 - ransq)

    res = np.empty((3), dtype=float)

    res[0] = ran_1 * ranh
    res[1] = ran_2 * ranh
    res[2] = 1.0 - 2.0 * ransq

    return res


def CdL(_input, elle):  ##mindist

        out=_input-np.rint(_input/elle)*elle
        return out

def C_pd_L(_input, elle): ##periodic coordinates

    out=_input-np.floor(_input/elle)*elle
    return out


def center_of_mass(myinput, box, mytype):

    outx = 0; outy = 0; outz = 0

    if mytype != 'pbc' and mytype != 'nopbc':
        print("only possible types: pbc and nopbc")
        exit(1)

    if mytype == 'pbc':
        xi_x = 0; xi_y = 0; xi_z = 0;
        zeta_x = 0; zeta_y = 0; zeta_z = 0;
        thisN = myinput.shape[0]

        for i in range(thisN):
            #outx = outx + input[i,0]; outy = outy + input[i,1]; outz = outz + input[i,2]; 
            tx = myinput[i,0]*2.*np.pi/box; xi_x = xi_x + np.cos(tx); zeta_x = zeta_x + np.sin(tx);
            ty = myinput[i,1]*2.*np.pi/box; xi_y = xi_y + np.cos(ty); zeta_y = zeta_y + np.sin(ty);
            tz = myinput[i,2]*2.*np.pi/box; xi_z = xi_z + np.cos(tz); zeta_z = zeta_z + np.sin(tz);

        xi_x = xi_x/float(thisN); xi_y = xi_y/float(thisN); xi_z = xi_z/float(thisN);
        zeta_x = zeta_x/float(thisN); zeta_y = zeta_y/float(thisN); zeta_z = zeta_z/float(thisN);

        outx = box*(np.arctan2(-zeta_x,-xi_x)+np.pi)*0.5/np.pi;
        outy = box*(np.arctan2(-zeta_y,-xi_y)+np.pi)*0.5/np.pi;
        outz = box*(np.arctan2(-zeta_z,-xi_z)+np.pi)*0.5/np.pi;

    if mytype == 'nopbc':
        outx = np.sum(myinput[:,0]); outy = np.sum(myinput[:,1]); outz = np.sum(myinput[:,2]);

        outx = outx/myinput.shape[0]; outy = outy/myinput.shape[0]; outz = outz/myinput.shape[0];

    return np.reshape(np.asarray([outx, outy, outz]),(3))


def gyration(_input, CoM, box):

    out = np.sum((_input[:]-CoM)**2)

    return np.sqrt(out/float(_input.shape[0]))


def Ree_correlation(_input):

    N = _input.shape[0]
    _out = np.zeros((N))
    count = np.zeros((N))

    #norm = np.sum(_input**2,axis=1)
    norm = np.linalg.norm(_input, axis=1)

    for i in range(N):
        _out[:N-i] += np.sum(_input[i:]*_input[i],axis=1)/(norm[i:]*norm[i])
        count[:N-i] += 1

    return _out/count


def persistence_length(_input):

        N = _input.shape[0]

        corr = np.zeros((N),dtype = float)
        count = np.zeros((N),dtype = float)

        norm = np.sum(_input**2,axis=1)

        for i in range(N):
            corr[:N-i] += np.sum(_input[i:]*_input[i],axis=1)/norm[i]
            count[:N-i] += 1

        corr /= count

        lp = np.average(np.sqrt(norm))*np.sum(corr[:int(N/2)])

        return corr, lp


def persistence_length_Ree(_input, mybox,Ree):

    out = np.empty((_input.shape[0]))

    norm = np.sum(_input**2,axis=1)

    out = np.sum(_input*Ree, axis=1)

    return out/norm



def bond_correlations(mydata, box):

    _N = mydata.shape[0];
    out = np.zeros((_N)); count = np.zeros((_N))
    bonds = np.empty((_N-1,4))

    bonds[:,:3] = CdL(mydata[1:,:] - mydata[:-1,:], box) 
    bonds[:,3] = np.sum(bonds[:,:3]*bonds[:,:3], axis=1)

    for i in range(_N-2):
                       
        for j in range(i+1,_N-1):
            out[j-i-1] += np.dot(bonds[i,:3],bonds[j,:3])/bonds[i,3]
            count[j-i-1] += 1

    return out, count


def avg_bending_angle(bonds, blength, out):

    _N = bonds.shape[0];print(_N)
    M = out.size; dth = 2./M    

    for i in range(_N-1):
        temp = np.sum(bonds[i+1]*bonds[i])/np.sqrt(blength[i+1]*blength[i]) ##np.cumsum(blength[i:])
        print((temp+1.)/dth)
        print(np.floor((temp+1.)/dth)-1)
        binn = int(np.floor((temp+1.)/dth)-0.001)
        out[binn] += 1.
	

##############################################################

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
        print("vectors are opposite: rotation not possible!")
        exit(1)

    if np.array_equal(v1,v2) == True:
        out = np.identity(3)
    else:
        v = np.cross(v1, v2)
        vmat = np.reshape(np.asarray([0,-v[2], v[1], v[2], 0, -v[0], -v[1], v[0], 0]), (3,3))
        out = np.identity(3) + vmat + np.dot(vmat,vmat)*(1./(1.+c))

    return v1, v2, out

def rotate_all(input_data, rot_matrix, myaxis):

    mydim = input_data.shape[0]
    out_data = np.empty_like(input_data)

    for i in range(mydim):
        temp = input_data[i] - myaxis
        outtemp = np.dot(rot_matrix,temp)
        out_data[i] = outtemp

    return out_data


def radial_density_cylinder(_input0, _inputbb, lp, mygen, Nbins):

    NN = _inputbb.shape[0]
    zaxis = np.asarray([0,0,1])

    allr = []; allth = []

    _step = int(2.**(mygen+1.)-1.)

    i = 0
    while(i != NN):
        i1 = i+lp+1
        
        if i1 > NN-1:
            i1 = NN-1
        
        vec = _inputbb[i1] - _inputbb[i]
       
        _input = _input0[i*_step:i1*_step]

        v1, v2, rm = rotation_matrix(vec, zaxis)
        newdata = rotate_all(_input, rm, _inputbb[i])
        #if i == 0:
        #    print_lorenzo(newdata, 100., "test.dat")
        #    print_lorenzo(_input, 100., "test0.dat")

        #temp0 = newdata[np.where(newdata[:,2] < i1)]
        #temp1 = temp0[np.where(temp0[:,2] > 0)]
        
        r = np.linalg.norm(newdata[:,:2], axis=1); allr.extend(r)
        theta = np.arctan2(newdata[:,1],newdata[:,0]); allth.extend(theta)
        
        i = i1+1
     
    rhist, bin_edges = np.histogram(allr, bins = Nbins, range=(0,10), density=False)
    rbins = 0.5*(bin_edges[1:]+bin_edges[:-1])

    thhist, bin_edges = np.histogram(allth, bins = Nbins, range=(-np.pi,np.pi), density=False)
    thbins = 0.5*(bin_edges[1:]+bin_edges[:-1])

    return [rbins, rhist], [thbins, thhist]
    
    
    
    
    
def print_data(_input, _params, _type):

    if _type == 'rdf':
        ## _params = [ nbins, dr, 4.*pi*dr*Npart, Nrep, filename]        
        if len(_params) != 5:
            print ('too few arguments for printing')
            return
            
        nbins = int(_params[0]); dr = float(_params[1]); w = float(_params[2]); Nrep = int(_params[3]);    filename = str(_params[4])
        mybox = 2.*nbins*dr
        #print "nbins %d, dr %.5f, box %.5f, else = %.5f" % (nbins, dr, mybox, w)


        toprint = np.zeros((nbins, 2), dtype = float)
        toprint[:,0] = np.linspace(dr/2.,mybox/2.+dr/2., nbins) #+ (dr/2.)*np.ones((nbins),dtype = float)
        weight = w*(np.power(toprint[:,0] + dr,3) - np.power(toprint[:,0],3))
        toprint[:,1] =  (1./(Nrep)) * _input / weight
        toprint[:,0] = toprint[:,0] + dr/2.
                #scale = np.average(toprint[-30:,1]); toprint[:,1] = toprint[:,1]/scale
        np.savetxt(filename, toprint, fmt = '%.8e')


    if _type == 'mixgrfourier':
        ##_params = [qmin, qmax, qbins, 4.*np.pi*dr*rho, rho, box2, filename]
                if len(_params) != 7:
                        print ('too few arguments for printing')
                        return

                qmin = float(_params[0]); qmax = float(_params[1]); qbins = int(_params[2]); weight = float(_params[3]); rho = float(_params[4]); box2 = float(_params[5]); filename = str(_params[6]) 
                dq = (qmax-qmin)/float(qbins)

                toprint = np.zeros((qbins, 2), dtype = float)
                toprint[:,0] = np.linspace(qmin,qmax, qbins) #+ (dq/2.)*np.ones((qbins),dtype = float)
                toprint[0,1] = _input[0]*weight #- (4./3.*np.pi)*rho*box2**3
                toprint[1:,1] = _input[1:]*weight #- (4./3.*np.pi)*rho*box2**3*Sksphere(box2*toprint[1:,0])
                np.savetxt(filename, toprint, fmt = '%.8e')

    
    if _type == 'grfourier':
        ##_params = [qmin, qmax, qbins, 4.*np.pi*dr*rho1, rho, box2, filename]
        if len(_params) != 7:
                        print ('too few/many arguments for printing')
                        return
    
        qmin = float(_params[0]); qmax = float(_params[1]); qbins = int(_params[2]); weight = float(_params[3]); rho = float(_params[4]); box2 = float(_params[5]); filename = str(_params[6])
        dq = (qmax-qmin)/float(qbins) 

        toprint = np.zeros((qbins, 2), dtype = float);print(toprint)
        toprint[:,0] = np.linspace(qmin, qmax, qbins) #np.linspace(qmin + dq/2.,qmax +dq/2., qbins) 
        toprint[0,1] = _input[0]*weight + 1 #- (4./3.*np.pi)*rho*box2**3
        toprint[1:,1] =  _input[1:]*weight + 1 #- (4./3.*np.pi)*rho*box2**3*Sksphere(box2*toprint[1:,0])
        np.savetxt(filename, toprint, fmt = '%.8e')

    if _type == 'formfactor':
        ##_params = [qmin, qmax, qbins, 4.*np.pi*dr*rho1, rho, box2, filename]
        if len(_params) != 7:
            print ('too few/many arguments for printing')
            return

        qmin = float(_params[0]); qmax = float(_params[1]); qbins = int(_params[2]); weight = float(_params[3]); rho = float(_params[4]); box2 = float(_params[5]); filename = str(_params[6])
        dq = (qmax-qmin)/float(qbins)

        toprint = np.zeros((qbins, 2), dtype = float);print(toprint)
        toprint[:,0] = np.linspace(qmin, qmax, qbins) #np.linspace(qmin + dq/2.,qmax +dq/2., qbins) 
        toprint[:,1] =  (_input[:]/_input[0])  ##*weight #+ 1 #- (4./3.*np.pi)*rho*box2**3*Sksphere(box2*toprint[1:,0])
        np.savetxt(filename, toprint, fmt = '%.8e')


    if _type == 'vanhove':
        
        if len(_params) != 2:
            print ('too few/many arguments for printing')
            return

        myx = np.copy(_params[0]); filename = str(_params[1])
        toprint = np.zeros((_input.shape[1], _input.shape[0]+1), dtype=float)
        toprint[:,0] = np.copy(myx)
        toprint[:,1:] = np.copy(_input.T)
        np.savetxt(filename, toprint, fmt = '%.8e')

    return toprint

