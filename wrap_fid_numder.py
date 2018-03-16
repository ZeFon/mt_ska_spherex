
# coding: utf-8

# In[3]:

import os
from configobj import ConfigObj
from numpy import *
import sys
folder=sys.argv[1]
file_root=sys.argv[2]

par_val_h={
    'H0':['hubble',67.74,0.7], #interval untested;
    'w':['w',-1,0.05],
    'Ob': ['omega_baryon',0.05,0.01],
    'Ocdm': ['omega_cdm',0.26,7.8e-03],
    'OL': ['omega_lambda',0.69,0.01],
    'On': ['omega_neutrino',0.0,0.001] #untested
    }

os.system('cp '+file_root+'_root.ini '+file_root+'.ini')

parametername=['Ocdm','Ob','w','H0']
fid=[]
h=[]
for i in range(len(parametername)):
    fid.append(par_val_h[parametername[i]][1])
    h.append(par_val_h[parametername[i]][2])


inifile=ConfigObj(file_root+'.ini')
l_max=int(inifile['l_max_scalar'])

inifile['output_root']=folder+file_root+'_fid'
inifile.write()
nbins=int(inifile['num_redshiftwindows'])

os.system('./camb '+file_root+'.ini')


def derivcl(interval,m2,m1,p1,p2):
    return (m2-8*m1+8*p1-p2)/(12.0*interval)


cl=zeros((l_max-1,nbins**2+1))
clraw=loadtxt(folder+file_root+'_fid_scalCovCls.dat')[:,1:]
for l in range(len(cl[:,0])):
    cl[l,0]=2+l
    cl[l,1:]=(clraw[l,:].reshape((nbins+3,nbins+3)))[3:,3:].reshape((nbins)**2)*2*pi/((2+l)*(l+3))
savetxt(folder+file_root+'_fid_Cl.dat',cl)

dcl=zeros((l_max-1,(nbins)**2+1))
ext=['_2m','_1m','_1p','_2p']
for j in range(len(parametername)):
    deriv_dic={}
    val=fid[j]+h[j]*array([-2,-1,1,2])
    for i in range(4):
        inifile['output_root']=folder+file_root+'_'+parametername[j]+ext[i]
        inifile[par_val_h[parametername[j]][0]]=str(val[i])
        #ensure flat universe. CDM implies a change on OL. Ob implies a change in Ocdm (for fixed O_m)
        if 'Ocdm'==parametername[j]:
            inifile['omega_lambda']=str(1.0-par_val_h['Ob'][1]-val[i])
        if 'OL'==parametername[j]:
            inifile['omega_cdm']=str(1.0-par_val_h['Ob'][1]-val[i])
        if 'Ob'==parametername[j]:
            inifile['omega_cdm']=str(par_val_h['Ocdm'][1]+par_val_h['Ob'][1]-val[i])
        inifile.write()
        os.system('./camb '+file_root+'.ini')
        deriv_dic[str(i)]=loadtxt(folder+file_root+'_'+parametername[j]+ext[i]+'_scalCovCls.dat')[:,1:]
    dclraw=derivcl(h[j],deriv_dic['0'],deriv_dic['1'],deriv_dic['2'],deriv_dic['3'])
    for l in range(len(dcl[:,0])):
        dcl[l,0]=2+l
        dcl[l,1:]=(dclraw[l,:].reshape((nbins+3,nbins+3)))[3:,3:].reshape((nbins)**2)*2*pi/((2+l)*(l+3))
    savetxt(folder+file_root+'_'+parametername[j]+'_dCl.dat',dcl)
    
    inifile[par_val_h[parametername[j]][0]]=str(fid[j])
    if 'Ocdm'==parametername[j]:
        inifile['omega_lambda']=str(1.0-par_val_h['Ob'][1]-par_val_h['Ocdm'][1])
    if 'OL'==parametername[j]:
        inifile['omega_cdm']=str(1.0-par_val_h['Ob'][1]-par_val_h['OL'][1])
    if 'Ob'==parametername[j]:
        inifile['omega_cdm']=str(par_val_h['Ocdm'][1])
    inifile.write()
    
inifile['output_root']=''
inifile.write()

#delete ini and scalcls
os.system('rm '+folder+'*.ini')
os.system('rm '+folder+'*_scalCls.dat')



