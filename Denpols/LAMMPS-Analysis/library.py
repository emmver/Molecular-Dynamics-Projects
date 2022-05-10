import numpy as np
import string
from numba import njit 
from scipy.stats import gaussian_kde
def read_lammpstrj(f):

    out = list()
    #f = open(inputfile, 'r')
    box = np.empty(3)

    trash = f.readline(); trash = f.readline();
    elems = str.split(trash.strip()," ");#print(elems)
    if str.isdigit(elems[0])==True:
        time = float(elems[0]);
        trash = f.readline(); trash = f.readline();
        elems = str.split(trash.strip()," ");  npart = float(elems[0]);npart=int(npart);#print("Particle Numbers",npart)
        trash = f.readline();
        trash = f.readline();
        elems = str.split(trash.strip()," ");  box[0] = float(elems[1])-float(elems[0])
        trash = f.readline();
        elems = str.split(trash.strip()," ");  box[1] = float(elems[1])-float(elems[0])
        trash = f.readline();
        elems = str.split(trash.strip()," ");  box[2] = float(elems[1])-float(elems[0]);
        trash = f.readline();
        i = 0
        for i in range(npart):
    
            line  = f.readline();
            elems = str.split(line.strip(" \r\n")," ");
            for el in elems:
                    #print ('el is',el)
                    out.extend([float(el)])
                    if len(elems) > 0 : elsize = len(elems)
    
        #f.close()
        out = np.reshape(np.asarray(out), (int(len(out)/elsize),elsize))
        nparticles = out.shape[0]
    
        return time, box, nparticles, out[out[:,0].argsort()]

        
        
        
       
    else:
        time=0
        npart = 0;npart=0;#print("Particle Numbers",npart)
        box[0] = 0
        box[1] = 0
        box[2] = 0
        elsize = len(elems)
        out = np.reshape(np.asarray(out), (int(len(out)/elsize),elsize))
        nparticles = np
        
        return time, box, nparticles, out
        

    

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

###########################################################################################

def CdL(_input, elle):

    out=_input-np.rint(_input/elle)*elle
    return out

@njit
def calc_rdf(data, output, dr, myL, nbins, fact):

        for i in range(data.shape[0]-1):

                mycenter = data[i,:]

                for j in range(i+1,data.shape[0],1):

                        ##d = np.linalg.norm(CdL(data[j,:]-mycenter[:], myL))#np.sqrt(dx*dx+dy*dy+dz*dz);
                        d = np.linalg.norm(data[j,:]-mycenter[:])

                        binn = int(np.floor(d/dr));
                        if binn < nbins:
                                output[binn] = output[binn] + 2*fact
#@njit
def calc_sq(data, qbins, dr, nbins,myL):

    out = np.zeros(qbins,dtype=float)
    qmin=2*np.pi/myL
    qmax=10
    dq=(qmax-qmin)/qbins
    print(dq,qbins*dq)
    
    #qmin=mylog()
    #q=np.logspace()

    #r = (np.arange(nbins)+0.5)*dr

    for i in range(qbins):
        q = (i+0.00001)*dq;

        for j in range(nbins):
            r = (j+0.5)*dr
            out[i] += r*(data[j])*np.sin(q*r)
            #out[i] = np.trapz(r*(data-1.0)*np.sin(q*r)/q,r)
        out[i] /= q
    print (q)

    return out,qmin,qmax

#@njit
def center_of_mass(myinput, box, mytype):

        outx = 0; outy = 0; outz = 0

        if mytype != 'pbc' and mytype != 'nopbc':
                print ("only possible types: pbc and nopbc")
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

#@njit
def gyration(input, CoM, box):

        out = 0

        for i in range(input.shape[0]):
                dx = input[i,0] - CoM[0]; dy = input[i,1] - CoM[1]; dz = input[i,2] - CoM[2];
                dx = CdL(dx,box[0]); dy = CdL(dy,box[1]); dz = CdL(dz,box[2]);
                out = out + (dx*dx+dy*dy+dz*dz)

        return np.sqrt(out/float(input.shape[0]))

def kde_scipy(x, x_grid, bandwidth=0.2, **kwargs):
    """Kernel Density Estimation with Scipy"""
    # Note that scipy weights its bandwidth by the covariance of the
    # input data.  To make the results comparable to the other methods,
    # we divide the bandwidth by the sample standard deviation here.
    kde = gaussian_kde(x, bw_method=bandwidth / x.std(ddof=1), **kwargs)
    return kde.evaluate(x_grid)