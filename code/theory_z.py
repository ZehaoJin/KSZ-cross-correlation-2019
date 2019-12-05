# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 00:08:07 2019

@author: zehaojin
"""



import numpy as np
import matplotlib.pyplot as plt
#import scipy
#from astropy.cosmology import FlatLambdaCDM
from astropy.cosmology import Planck15
##H0=67.7 km / (Mpc s), Om0=0.307, Tcmb0=2.725 K, Neff=3.05, m_nu=[0.   0.   0.06] eV, Ob0=0.0486)
import camb
from camb import model, initialpower
from scipy.special import eval_legendre
#from scipy.interpolate import Rbf
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d
import healpy
from astropy.coordinates import SkyCoord
from astropy.io import fits

def Cl(l,pk_interp,dn_dz,z_bins,zmin,zmax):
    ##
    delta_h=3
    ##
    Omega_m=0.307
    Omega_Lambda=1.0-Omega_m
    Omega_b=0.0486
    sigma_8=0.84
    T_cmb=2.725e6 #uK
    H0=67.7 #km/s/Mpc
    h=H0/100.0
    H_i=h*100.0
    H_e=H_i-delta_h
    sigma_T=6.6524e-29 #m^2
    c=3.0e5 #km/s
    b_g=1.18 #galaxy bias/linear bias
    X_e=1.0/1.6726219e-27   #1/kg  1proton/1 electron per unit mass,ignoring helium
    z_size=zmax-zmin
    cosmo=Planck15

    G=6.674e-11 #m^3 kg^-1 s^-2
    #1 megaParsec =3.08567758 × 10^19 kilometers
    Hs=H0/3.08567758e19  #H in s^-1
    rho_c=3*(Hs**2)/(8*np.pi*G)  #kg m^-3

    z_num=1000.0
    z=np.linspace(zmin,zmax,z_num+1)
    dz=z_size/(z_num)
    z=z[:-1]+dz/2.0

    D_A_co=cosmo.comoving_distance(z).value #Mpc
    v_H=(H_i-H_e)*D_A_co/(1+z)
    k=l/D_A_co  #Mpc^-1
    kh=k/h

    pk=np.zeros(int(z_num))
    for i in range(int(z_num)):
        pk[i]=pk_interp(kh[i],z[i])
    #pk=pk_interp(kh,z)
    #print kh.shape
    #print z.shape
    #print pk.shape
    P_e=pk*(h**3)

    Delta_e2=k**3/(2.0*np.pi**2)*P_e
    #print Delta_e2
    #Delta_e=np.sqrt(Delta_e2)
    #print Delta_e

    #http://www.iac.es/congreso/isapp2012/media/Longair-lectures/Longair2.pdf
    #dz_dt=-H0*(1.0+z)*np.sqrt((1.0+z)**2*(Omega_m*z+1.0)-Omega_Lambda*z*(z+2.0))
    dz_dt=-Hs*(1.0+z)*np.sqrt((1.0+z)**2*(Omega_m*z+1.0)-Omega_Lambda*z*(z+2.0))
    dt_dz=1.0/dz_dt
    n_e=Omega_b*rho_c*X_e #m^-3
    dtau_dz=sigma_T*n_e*c*1.0e3*dt_dz
    #print dtau_dz
    #integrand_T=(v_H)*(dtau_dz)*(np.sqrt(np.pi/l)*Delta_e2)*np.sqrt(D_A_co/(c/H_i))*(T_cmb/c)

    #integrand_g=b_g*(dn_dz)*(np.sqrt(np.pi/l))*np.sqrt(D_A_co/(c/H_i))*(T_cmb)

    #integrand_T=(v_H)*(dtau_dz)*(np.sqrt(np.pi/l)*Delta_e2)*np.sqrt(D_A_co/(c/H_i))*(1.0/c)

    #integrand_g=b_g*(dn_dz)*(np.sqrt(np.pi/l))*np.sqrt(D_A_co/(c/H_i))*(1.0)

    integrand=v_H/c*dtau_dz*(np.pi/l)*Delta_e2*b_g*dn_dz*(D_A_co*H_i/c)

    #Txg=np.sum(integrand_T*integrand_g*dz)
    Txg=np.sum(integrand*dz)

    C_l=2*np.pi*Txg/(l*(l+1.0))
    return C_l,Txg


#theta_max=10  #deg
#theta_bins=100
#theta_deg=np.linspace(0,theta_max,theta_bins)
l_max=4096


theta_deg=np.load('datalog/errorbar_r_SuperCosmos.npy')
theta_cos=np.cos(np.deg2rad(theta_deg))
theta_bins=int(theta_deg.size)

w_theta=np.zeros([l_max,theta_bins])
legendre=np.zeros([l_max,theta_bins])

zz=np.load('zz.npy')
kk=np.load('kk.npy')
pk=np.load('pk.npy')
pk_interp=interp2d(kk,zz,pk,kind='linear')

beamfile=fits.open('COM_CMB_IQU-smica_2048_R3.00_full.fits')
beam=beamfile[2].data['INT_BEAM'] #0~4096  shape(4097,)
beamfile.close()



for z_bins in range(4):

    catalog=np.load('datalog/wisecatalog.npy')

    scosmask=healpy.read_map('WISExSCOSmask.fits')
    coord=SkyCoord(catalog[0],catalog[1],frame='icrs',unit='deg').galactic
    l,b=coord.l.deg,coord.b.deg
    catalog=catalog[:,scosmask[healpy.ang2pix(256,l,b,nest=False,lonlat=True)]==1]
    #print catalog.shape
    catalog=catalog[:,catalog[2].argsort()]

    if z_bins==3:
        catalog=catalog[:,catalog[0].size/4*z_bins:]
    else:
        catalog=catalog[:,catalog[0].size/4*z_bins:catalog[0].size/4*(z_bins+1)]

    if z_bins==0:
        catalog=catalog[:,catalog[2]>=0.01]

    print('bin',z_bins,'size',catalog[2].shape,'z range',catalog[2,0],'~',catalog[2,-1])

    catalog=catalog[2]
    total_galaxy_number=float(catalog.size)
    #print catalog.shape
    z_bin=np.linspace(catalog[0],catalog[-1],1000.0+1)
    dn_dz,binedge=np.histogram(catalog,z_bin)

    print total_galaxy_number

    dn_dz=dn_dz/(total_galaxy_number*((catalog[-1]-catalog[0])/1000.0))
    #print dn_dz

    l_list=np.linspace(1,l_max,l_max)
    cl_list=np.zeros(l_max)

    grid_scale=1
    l_grid=np.linspace(1,l_max,(l_max-1)/grid_scale+1)
    cl_grid=np.zeros((l_max-1)/grid_scale+1)
    Txg_grid=np.zeros((l_max-1)/grid_scale+1)
    print 'generating Cl grid for z',z_bins
    for l in range(1,(l_max-1)/grid_scale+2):
        cl,Txg=Cl(l_grid[int(l-1)],pk_interp,dn_dz,z_bins,catalog[0],catalog[-1])
        cl_grid[int(l-1)]=cl
        Txg_grid[int(l-1)]=Txg
    cl_interp=interp1d(l_grid,cl_grid,kind='linear')


    print 'interp Cl...'
    for l in l_list:
        if l%100==0:
            print 'running calculation for l=',l
        cl=cl_interp(l)
        cl_list[int(l-1)]=cl
        bl=beam[int(l)]
        #w_theta+=(2*l+1)/(4*np.pi)*cl*eval_legendre(int(l),theta_cos)
        w_theta[int(l-1)]=(2*l+1)/(4*np.pi)*cl*bl*eval_legendre(int(l),theta_cos)

    cl_interp=None

    w_theta_sum=np.sum(w_theta,0)
    #w_theta_sum=w_theta_sum*(2.725*total_galaxy_number/(4*np.pi))
    w_theta_sum=w_theta_sum*(2.725)
    np.save('w_theta_z_{}.npy'.format(str(z_bins)),w_theta_sum)
    np.save('theta_deg_z_{}.npy'.format(str(z_bins)),theta_deg)
    #print Txg_grid,w_theta_sum
    #print theta_deg,w_theta_sum

    #'''
    plt.plot(l_list,cl_list,ls='--',label='linear interp',c='blue')
    plt.scatter(l_grid,cl_grid,label='grid',c='red')
    plt.xlabel('l')
    plt.ylabel('Cl')
    plt.legend()
    plt.show()


    plt.scatter(l_grid,Txg_grid,c='red')
    plt.xlabel('l')
    plt.ylabel('l(l+1)Cl/2pi')
    plt.ylim(-1e-7,0)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.show()




    plt.plot(theta_deg,w_theta_sum)
    plt.xscale('log')
    plt.xlabel(r'$\theta$ (degrees)')
    plt.ylabel(r'$w(\theta) (K)$')
    #plt.ylim(-1e-6,1e-6)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.title('theory prediction z {}'.format(str(z_bins)))

    plt.show()
    #'''
