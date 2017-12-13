
# coding: utf-8

# In[3]:

import os
from configobj import ConfigObj
import numpy as np

#general stuff
folder='mt_ska_spherex/'
file_root='isw_MT_hi_ha_dz02_z02_30'

inifile=ConfigObj(file_root+'.ini')

l_max=int(inifile['l_max_scalar'])
nw=int(inifile['num_redshiftwindows'])
terms=[]
for i in range(nw):
    if i<int(nw/2):
        terms=terms+['_bhi'+str(i+1)]
    else:
        terms=terms+['_bha'+str(i-int(nw/2)+1)]

inifile['DoTerms_per_Window']='T'
inifile['num_redshiftwindows']=str(nw+1)

for j in range(nw):
    inifile['counts_density('+str(j+1)+')']= 'T'
    inifile['counts_redshift('+str(j+1)+')']= 'T'
    inifile['counts_ISW('+str(j+1)+')']= 'T'
    inifile['counts_velocity('+str(j+1)+')']= 'T'
    inifile['counts_potential('+str(j+1)+')']= 'T'
    inifile['counts_evolve('+str(j+1)+')']= 'T' 

inifile['Dobder('+str(nw+1)+')']='T'
inifile.write()

for j in range(nw):
    inifile['output_root']=folder+file_root+terms[j]
    inifile['redshift('+str(nw+1)+')'] = inifile['redshift('+str(j+1)+')'] 
    inifile['redshift_sigma('+str(nw+1)+')'] = inifile['redshift_sigma('+str(j+1)+')']
    inifile['redshift_kind('+str(nw+1)+')'] = inifile['redshift_kind('+str(j+1)+')']
    inifile['redshift_wintype('+str(nw+1)+')'] = inifile['redshift_wintype('+str(j+1)+')']
    inifile['redshift_smooth('+str(nw+1)+')'] = inifile['redshift_smooth('+str(j+1)+')']
    inifile['redshift_bias('+str(nw+1)+')'] = inifile['redshift_bias('+str(j+1)+')']
    inifile['redshift_dlog10Ndm('+str(nw+1)+')'] = inifile['redshift_dlog10Ndm('+str(j+1)+')']
    inifile['redshift_dNdz('+str(nw+1)+')'] = inifile['redshift_dNdz('+str(j+1)+')']
            
    inifile.write()
    os.system('./camb '+file_root+'.ini')
    dcl=np.zeros((l_max-1,nw**2+1))
    ders_matrix=np.zeros((nw,nw))
    dclraw=np.loadtxt(folder+file_root+terms[j]+'_scalCovCls.dat')[:,1:]
    for l in range(len(dcl[:,0])):
        ders_matrix=np.zeros((nw,nw))
        dcl[l,0]=2+l
        ders_matrix[j,:]=dclraw[l,:].reshape((nw+4,nw+4))[-1,3:-1]*2*np.pi/((2+l)*(l+3))
        ders_matrix=ders_matrix+ders_matrix.T
        dcl[l,1:]=ders_matrix.reshape(nw**2)
    np.savetxt(folder+file_root+terms[j]+'_dCl.dat',dcl)
    inifile['output_root']=''
    inifile.write()


inifile['DoTerms_per_Window']='F'
inifile['num_redshiftwindows']=str(nw)
inifile['output_root']=''
inifile['Dobder('+str(nw+1)+')']='F'
inifile.write()

#delete ini and scalcls
os.system('rm '+folder+'*.ini')
os.system('rm '+folder+'*_scalCls.dat')



