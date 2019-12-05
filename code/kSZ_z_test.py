# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 00:08:07 2019

@author: zehaojin
"""


'''
http://ssa.roe.ac.uk/WISExSCOS.html
'''

#import matplotlib
#matplotlib.use('Qt5Agg')

from astropy.io import fits
import treecorr
import healpy
from astropy.coordinates import SkyCoord



import numpy as np
import matplotlib.pyplot as plt

'''
catalog=np.load('datalog/wisecatalog_z.npy')
print('z_min:',np.min(catalog[2]),'z_max:',np.max(catalog[2]))
plt.hist(catalog[2],100)
plt.xlabel('z')
plt.ylabel('# of galaxies')
plt.title('SuperCosmos z histogram')
plt.show()

plt.scatter(catalog[0],catalog[1],c=catalog[2],s=0.005,cmap='rainbow',edgecolors='none')
plt.xlabel('RA(deg)')
plt.ylabel('DEC(deg)')
plt.title('SuperCosmos redshift scatter plot')
plt.colorbar()
plt.show()

raw_input=()
'''

def written_as_a_function_to_save_memory(z_bins,randoms,result):
    print('loading SuperCosmos catalog and mask...')

    catalog=np.load('datalog/wisecatalog_z.npy') #ra,deg,z

    #print(catalog.shape)

    scosmask=healpy.read_map('WISExSCOSmask.fits')


    num=30000000
    #coord=SkyCoord(catalog[0],catalog[1],frame='galactic',unit='deg').icrs

    #catalog[0],catalog[1]=coord.ra.deg,coord.dec.deg
    coord=SkyCoord(catalog[0],catalog[1],frame='icrs',unit='deg').galactic
    l,b=coord.l.deg,coord.b.deg




    catalog=catalog[:,scosmask[healpy.ang2pix(256,l,b,nest=False,lonlat=True)]==1]

    '''
    print(catalog.shape)
    print('z_min:',np.min(catalog[2]),'z_max:',np.max(catalog[2]))
    plt.hist(catalog[2],100)
    plt.xlabel('z')
    plt.ylabel('# of galaxies')
    plt.title('SuperCosmos z histogram(masked)')
    plt.show()

    plt.scatter(catalog[0],catalog[1],c=catalog[2],s=0.005,cmap='rainbow',edgecolors='none')
    plt.xlabel('RA(deg)')
    plt.ylabel('DEC(deg)')
    plt.title('SuperCosmos redshift scatter plot(masked)')
    plt.colorbar()
    plt.show()

    raw_input=()
    '''
    catalog=catalog[:,catalog[2].argsort()]

    if z_bins==3:
        catalog=catalog[:,catalog[0].size/4*z_bins:]
    else:
        catalog=catalog[:,catalog[0].size/4*z_bins:catalog[0].size/4*(z_bins+1)]

    if z_bins==0:
        catalog=catalog[:,catalog[2]>=0.01]

    #print('bin',z_bins,'size',catalog[2].shape,'z range',catalog[2,0],'~',catalog[2,-1])
    #catalog=catalog[:,catalog[0].size]

    #cat_galaxy=treecorr.Catalog(ra=catalog[0],dec=catalog[1],ra_units='deg',dec_units='deg',k=np.ones(catalog[0].size))
    cat_galaxy=treecorr.Catalog(ra=catalog[0],dec=catalog[1],ra_units='deg',dec_units='deg')

    print('Done!\n')



    print 'generating random galaxy catalog'
    #plt.scatter(catalog[0],catalog[1],s=0.01)
    #plt.xlabel('RA(deg)')
    #plt.ylabel('DEC(deg)')
    #plt.show()
    ra_min = np.min(cat_galaxy.ra)
    ra_max = np.max(cat_galaxy.ra)
    dec_min = np.min(cat_galaxy.dec)
    dec_max = np.max(cat_galaxy.dec)
    print('ra range = %f .. %f' % (ra_min, ra_max))
    print('dec range = %f .. %f' % (dec_min, dec_max))

    rand_ra = np.random.uniform(ra_min, ra_max, num)
    rand_sindec = np.random.uniform(np.sin(dec_min), np.sin(dec_max), num)
    rand_dec = np.arcsin(rand_sindec)

    coord=SkyCoord(rand_ra,rand_dec,frame='icrs',unit='rad').galactic
    l,b=coord.l.deg,coord.b.deg



    rand_ra=rand_ra[scosmask[healpy.ang2pix(256,l,b,nest=False,lonlat=True)]==1]
    rand_dec=rand_dec[scosmask[healpy.ang2pix(256,l,b,nest=False,lonlat=True)]==1]

    #plt.scatter(np.rad2deg(rand_ra),np.rad2deg(rand_dec),s=0.01)
    #plt.xlabel('RA(deg)')
    #plt.ylabel('DEC(deg)')
    #plt.show()
    print('Done!\n')

    cat_rand = treecorr.Catalog(ra=rand_ra, dec=rand_dec, ra_units='radians', dec_units='radians')


    #load planck data
    print('loading Planck catalog and mask...')

    planckdata=fits.open('COM_CMB_IQU-smica-nosz_2048_R3.00_full.fits')

    planckmask=fits.open('HFI_Mask_GalPlane-apo0_2048_R2.00.fits')


    planckImap=planckdata[1].data['I_STOKES']

    planckImask=planckmask[1].data['GAL080']

    planckdata.close()
    planckmask.close()

    planckpix=np.arange(0,planckImap.size)

    planckImap=planckImap[planckImask==1]
    planckpix=planckpix[planckImask==1]

    planck_ra,planck_dec=healpy.pix2ang(nside=2048,ipix=planckpix,nest=True,lonlat=True)


    coord=SkyCoord(planck_ra,planck_dec,frame='galactic',unit='deg').icrs
    planck_ra,planck_dec=coord.ra.deg,coord.dec.deg

    cat_planck=treecorr.Catalog(ra=planck_ra,dec=planck_dec,ra_units='deg',dec_units='deg',k=planckImap)

    print('Done!\n')


    print('calculating cross-relation...')
    nk = treecorr.NKCorrelation(min_sep=0.01, max_sep=10, nbins=35, sep_units='deg')
    rk = treecorr.NKCorrelation(min_sep=0.01, max_sep=10, nbins=35, sep_units='deg')
    nk.process(cat_galaxy,cat_planck)
    rk.process(cat_rand,cat_planck)

    xi, varxi = nk.calculateXi(rk)
    sig = np.sqrt(varxi)
    r = np.exp(nk.meanlogr)

    #print xi


    print('Done!\n')


    #print('Plotting')


    #plt.plot(r, xi, color='blue')
    #plt.errorbar(r[xi>0], xi[xi>0], yerr=sig[xi>0], lw=1, ls='',ecolor='g')
    #leg = plt.errorbar(-r, xi, yerr=sig, color='blue')

    #plt.xscale('log')
    #plt.xlabel(r'$\theta$ (degrees)')
    #plt.ylabel(r'$w(\theta)$')
    #plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    #plt.legend([leg], [r'$w(\theta)$'], loc='lower left')
    #plt.title('SuperCosmos x Planck at {} z bin'.format(str(z_bins)))


    #plt.show()

    #np.save('datalog/SuperCosmos_z_{}.npy'.format(str(z_bins)),np.array([xi,r,sig]))
    result[z_bins,randoms]=np.array([xi,r,sig])
    print randoms,' runs'
    print('{} bins datalog saved!'.format(str(z_bins)))

    nk.clear()
    rk.clear()
    cat_galaxy.clear_cache()
    cat_rand.clear_cache()
    cat_planck.clear_cache()

    catalog=None
    scosmask=None
    coord=None
    l=None
    b=None
    ra_min,ra_max,dec_min,dec_max=None,None,None,None
    rand_ra,rand_sindec,rand_dec=None,None,None
    planck_ra,planck_dec=None,None
    planckImap,planckpix=None,None
    xi,r,sig,varxi=None,None,None,None


    return result


#result=np.zeros((4,50,3,35))
result=np.load('datalog/randomgalaxytest.npy')
for randoms in range(15,21):
    for z_bins in range(4):
        written_as_a_function_to_save_memory(z_bins,randoms,result)
        np.save('datalog/randomgalaxytest.npy',result)


np.save('datalog/SuperCosmos_z_0_reso.npy',np.mean(result[0],axis=0))
np.save('datalog/SuperCosmos_z_1_reso.npy',np.mean(result[1],axis=0))
np.save('datalog/SuperCosmos_z_2_reso.npy',np.mean(result[2],axis=0))
np.save('datalog/SuperCosmos_z_3_reso.npy',np.mean(result[3],axis=0))

#written_as_a_function_to_save_memory(0)
