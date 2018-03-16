
# coding: utf-8

# In[3]:

import os
from configobj import ConfigObj
import numpy as np
import sys
folder=sys.argv[1]
file_root=sys.argv[2]

os.system('cp '+file_root+'_root.ini '+file_root+'.ini')

inifile=ConfigObj(file_root+'.ini')

l_max=int(inifile['l_max_scalar'])


#for each subset

nbins=int(inifile['num_redshiftwindows'])
inifile['output_root']=folder+file_root+'_ns'
inifile['dn_s']='T'
inifile.write()
os.system('./camb '+file_root+'.ini')

dcl=np.zeros((l_max-1,nbins**2+1))
dclraw=np.loadtxt(folder+file_root+'_ns_scalCovCls.dat')[:,1:]
for l in range(len(dcl[:,0])):
    dcl[l,0]=2+l
    dcl[l,1:]=(dclraw[l,:].reshape((nbins+3,nbins+3)))[3:,3:].reshape((nbins)**2)*2*np.pi/((2+l)*(l+3))
np.savetxt(folder+file_root+'_ns_dCl.dat',dcl)

#delete ini and scalcls
os.system('rm '+folder+'*.ini')
os.system('rm '+folder+'*_scalCls.dat')
