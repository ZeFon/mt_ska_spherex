
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

terms=['_eISW','_edoppler','_eSW']


#for each subset

nw_raw=int(inifile['num_redshiftwindows'])
nw=2*nw_raw
inifile['DoTerms_per_Window']='T'
inifile['num_redshiftwindows']=str(nw)

for j in range(nw_raw):
    inifile['counts_density('+str(j+1)+')']= 'T'
    inifile['counts_redshift('+str(j+1)+')']= 'F'
    inifile['counts_ISW('+str(j+1)+')']= 'F'
    inifile['counts_velocity('+str(j+1)+')']= 'F'
    inifile['counts_potential('+str(j+1)+')']= 'F'
    inifile['counts_evolve('+str(j+1)+')']= 'T' 
inifile.write()
        

for j in range(nw_raw):
    inifile['redshift('+str(j+nw_raw+1)+')'] = inifile['redshift('+str(j+1)+')']
    inifile['redshift_sigma('+str(j+nw_raw+1)+')'] = inifile['redshift_sigma('+str(j+1)+')']
    inifile['redshift_kind('+str(j+nw_raw+1)+')'] = inifile['redshift_kind('+str(j+1)+')']
    inifile['redshift_wintype('+str(j+nw_raw+1)+')'] = inifile['redshift_wintype('+str(j+1)+')']  
    inifile['redshift_smooth('+str(j+nw_raw+1)+')'] = inifile['redshift_smooth('+str(j+1)+')']
    inifile['redshift_bias('+str(j+nw_raw+1)+')'] = inifile['redshift_bias('+str(j+1)+')']
    inifile['redshift_dlog10Ndm('+str(j+nw_raw+1)+')'] = inifile['redshift_dlog10Ndm('+str(j+1)+')']
    inifile['redshift_dNdz('+str(j+nw_raw+1)+')'] = inifile['redshift_dNdz('+str(j+1)+')']
inifile.write()

for i in range(len(terms)):
    inifile['output_root']=folder+file_root+'_den'+terms[i]
    
    for j in range(nw_raw):
        inifile['counts_evolve('+str(j+nw_raw+1)+')']= 'T'
        if terms[i]=='_eISW':
            inifile['counts_ISW('+str(j+nw_raw+1)+')']= 'T'
            inifile['counts_ISW('+str(j+1)+')']= 'T'
        elif terms[i]=='_edoppler':
            inifile['counts_velocity('+str(j+nw_raw+1)+')']= 'T'
            inifile['counts_velocity('+str(j+1)+')']= 'T'
        elif terms[i]=='_eSW':
            inifile['counts_potential('+str(j+nw_raw+1)+')']= 'T'
            inifile['counts_potential('+str(j+1)+')']= 'T'
            
    inifile.write()
    os.system('./camb '+file_root+'.ini')
    cl=np.zeros((l_max-1,nw_raw**2+1))
    dcl=np.zeros((l_max-1,nw_raw**2+1))
    dclraw=np.loadtxt(folder+file_root+'_den'+terms[i]+'_scalCovCls.dat')[:,1:]
    ders_matrix=np.zeros((nw_raw,nw_raw))
    for l in range(len(dcl[:,0])):
        dcl[l,0]=2+l
        fid_matrix=(dclraw[l,:].reshape((nw+3,nw+3))[3:(3+nw_raw),3:(3+nw_raw)])*2*np.pi/((2+l)*(l+3))
        ders_matrix=(dclraw[l,:].reshape((nw+3,nw+3))[(3+nw_raw):,3:(nw_raw+3)])*2*np.pi/((2+l)*(l+3))
        ders_matrix=ders_matrix+ders_matrix.T
        dcl[l,1:]=ders_matrix.reshape(nw_raw**2)
        cl[l,1:]=fid_matrix.reshape(nw_raw**2)
    np.savetxt(folder+file_root+'_den'+terms[i]+'_dCl.dat',dcl)
    np.savetxt(folder+file_root+'_den'+terms[i]+'_fid_Cl.dat',dcl)
    inifile['output_root']=''
    for j in range(nw_raw):
        if terms[i]=='_eISW':
            inifile['counts_ISW('+str(j+nw_raw+1)+')']= 'F'
            inifile['counts_ISW('+str(j+1)+')']= 'F'
        elif terms[i]=='_edoppler':
            inifile['counts_velocity('+str(j+nw_raw+1)+')']= 'F'
            inifile['counts_velocity('+str(j+1)+')']= 'F'
        elif terms[i]=='_eSW':
            inifile['counts_potential('+str(j+nw_raw+1)+')']= 'F'
            inifile['counts_potential('+str(j+1)+')']= 'F'
    inifile.write()

inifile['DoTerms_per_Window']='F'
inifile['num_redshiftwindows']=str(nw_raw)
inifile.write()

#delete ini and scalcls
os.system('rm '+folder+'*.ini')
os.system('rm '+folder+'*_scalCls.dat')
