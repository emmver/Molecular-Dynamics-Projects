import glob 
import numpy 
import sys
import os

os.chdir(sys.argv[1])

file_name='all.vtk'
lambdas_w=[2,6,12]
DP=200
to_read=4*DP+5+4+3+3

for i in range (1,9):
    for idx in lambdas_w: 
        print('RUN:%d/lambda:%d'%(i,idx))
        file='./RUN_%d/lambda_%d/vtk_equil/'%(i,idx)+file_name
#        keep_track=0
        g=open(file,'r')
        to_store=os.path.dirname(file)
        for frame in range(1000):
            file_to_write='part_test_%d.vtk'%frame
            f=open(to_store+'/'+file_to_write,'w')
            for ii in range(to_read):
                line=g.readline()
                f.write(line)
#            keep_track+=1
            f.close()
        g.close()
