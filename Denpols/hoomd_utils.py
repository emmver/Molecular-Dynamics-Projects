
######### Calculate Gyration Tensor #######################3
@njit(fastmath=True)
def rg_tens (r):
    
    N = r.shape[0]
    r2=np.zeros((3,3))
    for i in range (N):
        for j in range (N):
            r2[0,0]+=(r[i,0]-r[j,0])*(r[i,0]-r[j,0])#diag
            r2[1,0]+=(r[i,1]-r[j,0])*(r[i,1]-r[j,0])#ndiag
            r2[2,0]+=(r[i,2]-r[j,0])*(r[i,2]-r[j,0])#ndiag
            r2[0,1]+=(r[i,0]-r[j,1])*(r[i,0]-r[j,1])#ndiag
            r2[0,2]+=(r[i,0]-r[j,2])*(r[i,0]-r[j,2])#ndiag
            r2[1,1]+=(r[i,1]-r[j,1])*(r[i,1]-r[j,1])#diag
            r2[2,1]+=(r[i,2]-r[j,1])*(r[i,2]-r[j,1])#ndiag
            r2[1,2]+=(r[i,1]-r[j,2])*(r[i,1]-r[j,2])#ndiag
            r2[2,2]+=(r[i,2]-r[j,2])*(r[i,2]-r[j,2])# diag
    r2=r2/(2*N**2)
    #rg2=r2[0,0]**2+r2[1,1]**2+r2[2,2]**2
    return r2

def plot_shape(num_frames,gyr_tens):
    frames=np.arange(num_frames)+1
    plt.plot(frames,gyr_tens[0]**0.5,color='#cc0000')
    plt.xlabel('Frames')
    plt.ylabel('$R_g [\sigma]$')

    plt.show()
    plt.clf()
    ### Eigenvalues - Normalized ####
    colors = pl.cm.jet(np.linspace(0,0.8,3))
    plt.plot(frames,gyr_tens[1]/gyr_tens[0],color=colors[0],label=r'$\lambda_1$')
    plt.plot(frames,gyr_tens[2]/gyr_tens[0],color=colors[1],label=r'$\lambda_2$')
    plt.plot(frames,gyr_tens[3]/gyr_tens[0],color=colors[2],label=r'$\lambda_3$')
    plt.xlabel('Frames')
    plt.ylabel(r'$\frac{\lambda_i}{R_g^2}$')
    plt.legend(frameon=False,loc=(0.2,1.02),ncol=3)
    plt.show()
    plt.clf()
    ### Ratio of Eigenvalues ###
    colors = pl.cm.jet(np.linspace(0,0.8,3))
    plt.plot(frames,gyr_tens[1]/gyr_tens[2],color=colors[0],label=r'$\frac{\lambda_1}{\lambda_2}$')
    plt.plot(frames,gyr_tens[2]/gyr_tens[3],color=colors[1],label=r'$\frac{\lambda_2}{\lambda_3}$')
    plt.plot(frames,gyr_tens[3]/gyr_tens[1],color=colors[2],label=r'$\frac{\lambda_3}{\lambda_1}$')
    plt.xlabel('Frames')
    plt.ylabel(r'$\frac{\lambda_i}{\lambda_j}$')
    plt.legend(frameon=False,loc=(0.8,0.53))
    plt.show()
    plt.clf()

    #### Shape Parameters - Normalized #####
    b=1.5*gyr_tens[2] - gyr_tens[0]/2
    c=gyr_tens[2]-gyr_tens[3]
    k=(b**2 + 0.75 * c**2)/gyr_tens[0]**2
    plt.plot(frames,b/gyr_tens[0],color=colors[0],label=r'b')
    plt.plot(frames,c/gyr_tens[0],color=colors[1],label=r'c')
    plt.plot(frames,k,color=colors[2],label=r'$\kappa^2 $')
    plt.xlabel('Frames')
    plt.ylabel(r'Shape Parameters')
    plt.legend(frameon=False,loc=(0.05,0.05),ncol=1)
    plt.show()
    plt.clf()
################################################################


###### k-grid ##### 

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

def normsk(inputx, inputy, d):

        nk = inputx.size
        _out = np.zeros((nk,2), dtype = np.double)

        k = 0; count = 0; j = 0

        for i in range(0,inputx.size-1):
                count += 1
                _out[k,0] += inputx[i]
                _out[k,1] += inputy[i]

                if inputx[i+1] - inputx[j] > d or i+1 == inputx.size-1:
                        _out[k,0] = _out[k,0]/count
                        _out[k,1] = _out[k,1]/count
                        #print "for k = %d count was %d, outx is %.5f" % (k, count, _out[k,0])
                        count = 0; k = k+1; d *= 1.05
        return _out[:k,:]

def generate_kgrid(dqmult, mybox, qmax):

        nqbound = int(np.floor((qmax*mybox)/(2.*np.pi*dqmult*1.732050808))); print(nqbound)
        Nq = nqbound**3; print("how many k", Nq)

        output = np.empty((Nq,4), dtype = np.float)

        dkx = 2.0*np.pi*dqmult/mybox; dky = 2.0*np.pi*dqmult/mybox; dkz = 2.0*np.pi*dqmult/mybox;

        kk = int(0)
        for ix in tqdm(range(nqbound)):
                for iy in range(nqbound):
                        for iz in range(nqbound):
                                output[kk,0] = ix*dkx; output[kk,1] = iy*dky; output[kk,2] = iz*dkz;
                                output[kk,3] = np.sqrt( (ix*dkx)**2 + (iy*dky)**2 + (iz*dkz)**2 )
                                kk += 1

        output = output[output[:,3].argsort()]
        thisqmax = output[-1,3]
        print ("qmax = %.5f !!!" % thisqmax)

        return np.ravel(output[:,:3]), output[:,3], Nq


def normsk_rand(inputx, inputy):

        nk = inputy.size
        _out = np.zeros((nk,2), dtype = np.double)

        k = 0; count = 0
        j = 0
        for i in range(0,nk-1):
                count += 1
                _out[k,0] += inputx[i]
                _out[k,1] += inputy[i]

                if (i+1)%100 == 0:
                        _out[k,0] = _out[k,0]/count
                        _out[k,1] = _out[k,1]/count
                        count = 0; k = k+1; 
                        

        return _out[:(k-1),:]


def generate_kgrid_rand (NqL, q0, _ff):

        qval = np.empty((NqL),dtype=float)
        d = q0
        for i in tqdm(range(NqL),desc="Loop 1"):
                qval[i] = d
                d *= _ff

        output = np.empty((NqL*100,4), dtype = np.float)
        
        for i in tqdm(range(NqL),desc="Loop 2 (long)"):
                for j in range(100):
                        output[100*i+j,:3] = qval[i]*random_on_sphere(); output[100*i+j,3] = qval[i]

        #print "qmax = %.5f !!!" % qval[-1]

        return np.ravel(output[:,:3]), output[:,3],  NqL*100
################################################################################################


### Data printing 
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

#################################################################3
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


def Ree_correlation(_input):

    N = _input.shape[0]
    _out = np.zeros((N))
    count = np.zeros((N))

    #norm = np.sum(_input**2,axis=1)
    norm = np.linalg.norm(_input, axis=1)

    for i in tqdm(range(N),desc="End-to-End Correlation"):
        _out[:N-i] += np.sum(_input[i:]*_input[i],axis=1)/(norm[i:]*norm[i])
        count[:N-i] += 1

    return _out/count

def avg_bending_angle(bonds, blength, out):

    _N = bonds.shape[0];print(_N)
    M = out.size; dth = 2./M    

    for i in range(_N-1):
        temp = np.sum(bonds[i+1]*bonds[i])/np.sqrt(blength[i+1]*blength[i]) ##np.cumsum(blength[i:])
        print((temp+1.)/dth)
        print(np.floor((temp+1.)/dth)-1)
        binn = int(np.floor((temp+1.)/dth)-0.001)
        out[binn] += 1.


def analyze_sk(data):

  
    unwrap = data.particles.positions[:] 
    CoM = center_of_mass(unwrap, mybox, 'nopbc')
    unwrap2 = np.copy(unwrap); unwrap2[:] -= CoM
    _sk.calc_sk (ctypes.c_void_p(kgrid.ctypes.data), ctypes.c_void_p(unwrap.ctypes.data), ctypes.c_void_p(Sskarr.ctypes.data), ctypes.c_int(nparticles), ctypes.c_int(nkbins))
    types=data.particles.particle_types[:]
    temp=unwrap[np.where(types==0)]
    backbone=temp
    sk.calc_sk (ctypes.c_void_p(kgrid.ctypes.data), ctypes.c_void_p(backbone.ctypes.data), ctypes.c_void_p(Sskbb.ctypes.data), ctypes.c_int(Lbackbone), ctypes.c_int(nkbins))
    bonds = backbone[1:] - backbone[:-1]
    blength = np.sum(bonds[:]**2, axis=1)    
    avg_bending_angle(bonds, blength, btheta)


def normsk_rand(inputx, inputy):

        nk = inputy.size
        _out = np.zeros((nk,2), dtype = np.double)

        k = 0; count = 0
        j = 0
        for i in range(0,nk-1):
                count += 1
                _out[k,0] += inputx[i]
                _out[k,1] += inputy[i]

                if (i+1)%100 == 0:
                        _out[k,0] = _out[k,0]/count
                        _out[k,1] = _out[k,1]/count
                        count = 0; k = k+1; 
                        

        return _out[:(k-1),:]


def norm_sk(totrep):
    for i in tqdm(range(nkbins),desc='Norm Full Sk'):
        Ssktrans[i] = (1./(1.*totrep)) * Sskarr[i] / float(nparticles)**2
    toprint = normsk_rand(knorm, Ssktrans)
    filename=_tostore+"/G{}{}".format(mygen,Lbackbone)
    print("To Store:",filename)
    np.savetxt(filename, toprint, fmt = '%.8e')
    for i in tqdm(range(nkbins),desc='Norm Backbone Sk'):
        Ssktrans[i] = (1./(1.*totrep)) * Sskbb[i] / float(Lbackbone)**2
    toprint = normsk_rand(knorm, Ssktrans)
    filename=_tostore+"/G{}{}_bb".format(mygen,Lbackbone)
    print("To Store:",filename)


