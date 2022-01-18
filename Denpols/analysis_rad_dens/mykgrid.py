import numpy as np

from utilities import random_on_sphere

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
        for ix in range(nqbound):
                for iy in range(nqbound):
                        for iz in range(nqbound):
                                output[kk,0] = ix*dkx; output[kk,1] = iy*dky; output[kk,2] = iz*dkz;
                                output[kk,3] = np.sqrt( (ix*dkx)**2 + (iy*dky)**2 + (iz*dkz)**2 )
                                kk += 1

        output = output[output[:,3].argsort()]
        thisqmax = output[-1,3]
        print ("qmax = %.5f !!!" % thisqmax)

        #if thisqmax <= qmax:
        #       dkmult = qmax/thisqmax
        #       output = output*dkmult

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


def generate_kgrid_rand(NqL, q0, _ff):

        qval = np.empty((NqL),dtype=float)
        d = q0
        for i in range(NqL):
                qval[i] = d
                d *= _ff

        output = np.empty((NqL*100,4), dtype = np.float)
        
        for i in range(NqL):
                for j in range(100):
                        output[100*i+j,:3] = qval[i]*random_on_sphere(); output[100*i+j,3] = qval[i]

        #print "qmax = %.5f !!!" % qval[-1]

        return np.ravel(output[:,:3]), output[:,3],  NqL*100



