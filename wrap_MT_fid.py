
# coding: utf-8

# In[3]:

import os
from configobj import ConfigObj
from numpy import *

file_root='isw_MT_hi_ha_dz02_z02_30'
folder='mt_ska_spherex/'
root='isw_MT_hi_ha_dz02_z02_30'

parametername=['Ocdm_','ns_','Ob_','w_','H0_']
h=array([7.8e-03,0.01,0.01,0.05,0.7])
fid=array([0.26,0.9667,0.05,-1.0,67.74])


inifile=ConfigObj(file_root+'.ini')
inifile['scalar_amp(1)']=str(2.142e-09)
inifile['omega_cdm']=str(fid[0])
inifile['omega_lambda']=str(1.0-fid[2]-fid[0])
inifile['scalar_spectral_index(1)']=str(fid[1])
inifile['omega_baryon']=str(fid[2])
inifile['w']=str(fid[3])
inifile['hubble']=str(fid[4])

l_max=600
inifile['l_max_scalar']=str(l_max)

inifile['output_root']=folder+root+'_fid'
#for each
nT=2
nw=14 
nbins=nT*nw

inifile['num_redshiftwindows'] = str(nbins)
inifile['DoTerms_per_Window']='F'
inifile['counts_density']= 'T'
inifile['counts_redshift']= 'T'
inifile['counts_ISW']= 'T'
inifile['counts_velocity']= 'T'
inifile['counts_potential']= 'T'
inifile['counts_evolve']= 'T' 
inifile.write()
os.system('./camb '+file_root+'.ini')


def derivcl(interval,m2,m1,p1,p2):
    return (m2-8*m1+8*p1-p2)/(12.0*interval)


cl=zeros((l_max-1,(nbins)**2+1))
clraw=loadtxt(folder+root+'_fid_scalCovCls.dat')[:,1:]
for l in range(len(cl[:,0])):
    cl[l,0]=2+l
    cl[l,1:]=(clraw[l,:].reshape((nbins+3,nbins+3)))[3:,3:].reshape((nbins)**2)*2*pi/((2+l)*(l+3))
savetxt(folder+root+'_fid_Cl.dat',cl)

dcl=zeros((l_max-1,(nbins)**2+1))
ext=['2m','1m','1p','2p']
for j in range(len(parametername)):
    deriv_dic={}
    val=fid[j]+h[j]*array([-2,-1,1,2])
    for i in range(4):
        inifile['output_root']=folder+root+'_'+parametername[j]+ext[i]
        if j==0:
            inifile['omega_cdm']=str(val[i])
            inifile['omega_lambda']=str(1.0-fid[2]-val[i])
        if j==1:
            inifile['scalar_spectral_index(1)']=str(val[i])
        if j==2:
            inifile['omega_baryon']=str(val[i])
            inifile['omega_cdm']=str(fid[2]+fid[0]-val[i])
        if j==3:
            inifile['w']=str(val[i])
        if j==4:
            inifile['hubble']=str(val[i])
        inifile.write()
        os.system('./camb '+file_root+'.ini')
        deriv_dic[str(i)]=loadtxt(folder+file_root+'_'+parametername[j]+ext[i]+'_scalCovCls.dat')[:,1:]
    dclraw=derivcl(h[j],deriv_dic['0'],deriv_dic['1'],deriv_dic['2'],deriv_dic['3'])
    for l in range(len(dcl[:,0])):
        dcl[l,0]=2+l
        dcl[l,1:]=(dclraw[l,:].reshape((nbins+3,nbins+3)))[3:,3:].reshape((nbins)**2)*2*pi/((2+l)*(l+3))
    savetxt(folder+root+'_'+parametername[j]+'dCl.dat',dcl)
    inifile['scalar_amp(1)']=str(2.142e-09)
    inifile['omega_cdm']=str(fid[0])
    inifile['omega_lambda']=str(1.0-fid[2]-fid[0])
    inifile['scalar_spectral_index(1)']=str(fid[1])
    inifile['omega_baryon']=str(fid[2])
    inifile['w']=str(fid[3])
    inifile['hubble']=str(fid[4])
    inifile.write()
    
inifile['output_root']=''
inifile.write()

#delete ini and scalcls
os.system('rm '+folder+'*.ini')
os.system('rm '+folder+'*_scalCls.dat')



