# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 20:53:04 2019

@author: zehaojin
"""

#import pyfits
import treecorr
import healpy
from astropy.coordinates import SkyCoord
from astropy.io import fits
import astropy.utils.data

import numpy as np
import matplotlib.pyplot as plt
import time

#import wget
#import os



def loadGalaxy():
    print 'loading SuperCosmos catalog and mask...'

    catalog=np.load('datalog/wisecatalog.npy')
    scosmask=healpy.read_map('WISExSCOSmask.fits',verbose=False)

    num=30000000

    coord=SkyCoord(catalog[0],catalog[1],frame='icrs',unit='deg').galactic
    l,b=coord.l.deg,coord.b.deg

    catalog=catalog[:,scosmask[healpy.ang2pix(256,l,b,nest=False,lonlat=True)]==1]

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
    #print('ra range = %f .. %f' % (ra_min, ra_max))
    #print('dec range = %f .. %f' % (dec_min, dec_max))

    rand_ra = np.random.uniform(ra_min, ra_max, num)
    rand_sindec = np.random.uniform(np.sin(dec_min), np.sin(dec_max), num)
    rand_dec = np.arcsin(rand_sindec)

    coord=SkyCoord(rand_ra,rand_dec,frame='icrs',unit='rad').galactic
    l,b=coord.l.deg,coord.b.deg



    rand_ra=rand_ra[scosmask[healpy.ang2pix(256,l,b,nest=False,lonlat=True)]==1]
    rand_dec=rand_dec[scosmask[healpy.ang2pix(256,l,b,nest=False,lonlat=True)]==1]

    cat_rand = treecorr.Catalog(ra=rand_ra, dec=rand_dec, ra_units='radians', dec_units='radians')

    #plt.scatter(np.rad2deg(rand_ra),np.rad2deg(rand_dec),s=0.01)
    #plt.xlabel('RA(deg)')
    #plt.ylabel('DEC(deg)')
    #plt.show()
    print('Done!\n')

    return cat_galaxy,cat_rand

def loadCluster():
    print('loading MCXC catalog and mask...')

    num=1000000

    MCXCfile=fits.open('MCXC.fits')
    MCXCdata=MCXCfile[1].data
    MCXCfile.close()

    MCXCmaskfile=fits.open('HFI_Mask_GalPlane-apo0_2048_R2.00.fits')
    MCXCmask=MCXCmaskfile[1].data['GAL080']
    MCXCmaskfile.close()

    coord=SkyCoord(MCXCdata['RA'],MCXCdata['DEC'],frame='icrs',unit='deg').galactic
    l,b=coord.l.deg,coord.b.deg

    MCXCdata=MCXCdata[MCXCmask[healpy.ang2pix(2048,l,b,nest=True,lonlat=True)]==1]

    #ra260-280 dec60-70 suspicious area
    MCXCmask2=np.ones(MCXCdata.size)
    for i in range(MCXCdata.size):
        if (260<MCXCdata[i]['RA']<280 and 60<MCXCdata[i]['DEC']<70):
            MCXCmask2[i]=0
    MCXCdata=MCXCdata[MCXCmask2==1]

    cat_galaxy=treecorr.Catalog(ra=MCXCdata['RA'],dec=MCXCdata['DEC'],ra_units='deg',dec_units='deg')
    #cat_galaxy=treecorr.Catalog(ra=MCXCdata['RA'],dec=MCXCdata['DEC'],ra_units='deg',dec_units='deg',k=np.ones(MCXCdata['RA'].size))
    print('Done!\n')

    #plt.scatter(MCXCdata['RA'],MCXCdata['DEC'],s=0.5)
    #plt.xlabel('RA(deg)')
    #plt.ylabel('DEC(deg)')
    #plt.show()
    print 'generating random galaxy catalog'

    ra_min = np.min(cat_galaxy.ra)
    ra_max = np.max(cat_galaxy.ra)
    dec_min = np.min(cat_galaxy.dec)
    dec_max = np.max(cat_galaxy.dec)
    print('ra range = %f .. %f' % (ra_min, ra_max))
    print('dec range = %f .. %f' % (dec_min, dec_max))

    rand_ra = np.random.uniform(ra_min, ra_max, num)
    rand_sindec = np.random.uniform(np.sin(dec_min), np.sin(dec_max), num)
    rand_dec = np.arcsin(rand_sindec)

    #plt.hist(rand_ra)
    #plt.show()
    #plt.hist(rand_sindec)
    #plt.show()

    planckmask=fits.open('HFI_Mask_GalPlane-apo0_2048_R2.00.fits')
    planckImask=planckmask[1].data['GAL080']
    planckmask.close()

    coord=SkyCoord(rand_ra,rand_dec,frame='icrs',unit='rad').galactic
    l,b=coord.l.deg,coord.b.deg

    rand_ra=rand_ra[planckImask[healpy.ang2pix(2048,l,b,nest=True,lonlat=True)]==1]
    rand_dec=rand_dec[planckImask[healpy.ang2pix(2048,l,b,nest=True,lonlat=True)]==1]

    #ra260-280 dec60-70 suspicious area
    MCXCrandmask2=np.ones(rand_ra.size)
    for i in range(rand_ra.size):
        if (260<np.rad2deg(rand_ra[i])<280 and 60<np.rad2deg(rand_dec[i])<70):
            MCXCrandmask2[i]=0
    rand_ra=rand_ra[MCXCrandmask2==1]
    rand_dec=rand_dec[MCXCrandmask2==1]

    cat_rand = treecorr.Catalog(ra=rand_ra, dec=rand_dec, ra_units='radians', dec_units='radians')
    #plt.scatter(np.rad2deg(rand_ra),np.rad2deg(rand_dec),s=0.5)
    #plt.scatter(MCXCdata['RA'],MCXCdata['DEC'],s=0.5)
    #plt.xlabel('RA(deg)')
    #plt.ylabel('DEC(deg)')
    #plt.show()
    #input('stop here')
    print('Done!\n')

    return cat_galaxy,cat_rand

def Calculate(cat_galaxy,cat_rand,result,num_of_runs,i,planck_ra,planck_dec):
    start_time=time.time()
    print 'run ',i,'/',num_of_runs
    print 'loading simulation and mask...'
    #download
    #wget.download('http://pla.esac.esa.int/pla/aio/product-action?SIMULATED_MAP.FILE_ID=dx12_v3_smica_nosz_cmb_mc_{}_raw.fits'.format(str(i).zfill(5)),'data/')
    #wget.download('http://pla.esac.esa.int/pla/aio/product-action?SIMULATED_MAP.FILE_ID=dx12_v3_smica_nosz_noise_hm1_mc_{}_raw.fits'.format(str(i).zfill(5)),'data/')

    #plancksim_n=healpy.read_map('data/dx12_v3_smica_nosz_cmb_mc_{}_raw.fits'.format(str(i).zfill(5)))
    #plancknoise=healpy.read_map('data/dx12_v3_smica_nosz_noise_hm1_mc_{}_raw.fits'.format(str(i).zfill(5)))

    plancksim_n_file=fits.open('http://pla.esac.esa.int/pla/aio/product-action?SIMULATED_MAP.FILE_ID=dx12_v3_smica_nosz_cmb_mc_{}_raw.fits'.format(str(i).zfill(5)),cache=True,ignore_missing_end=True,show_progress=False)
    plancknoise_file=fits.open('http://pla.esac.esa.int/pla/aio/product-action?SIMULATED_MAP.FILE_ID=dx12_v3_smica_nosz_noise_hm1_mc_{}_raw.fits'.format(str(i).zfill(5)),cache=True,ignore_missing_end=True,show_progress=False)

    plancksim_n=healpy.read_map(plancksim_n_file,verbose=False)
    plancknoise=healpy.read_map(plancknoise_file,verbose=False)

    plancksim_n_file.close()
    plancknoise_file.close()

    plancksim=plancksim_n+plancknoise

    #print plancksim_n.shape, plancknoise.shape,plancksim.shape

    plancksim_n=None
    plancknoise=None
    plancksim=healpy.ud_grade(plancksim,2048, order_in='RING', order_out='NEST')

    plancksim=plancksim[planckImask==1]

    cat_sim=treecorr.Catalog(ra=planck_ra,dec=planck_dec,ra_units='deg',dec_units='deg',k=plancksim)


    print 'Calculating correlation...'
    print 'run ',i,'/',num_of_runs

    nk = treecorr.NKCorrelation(min_sep=0.01, max_sep=10, nbins=35, sep_units='deg')
    rk = treecorr.NKCorrelation(min_sep=0.01, max_sep=10, nbins=35, sep_units='deg')
    nk.process(cat_galaxy,cat_sim)
    rk.process(cat_rand,cat_sim)

    xi, varxi = nk.calculateXi(rk)


    print 'Done'

    result[i]=xi
    r = np.exp(nk.meanlogr)

    if catalogname=='1':
    	np.save('datalog/errorbar_SuperCosmos.npy',result)
    	np.save('datalog/errorbar_r_SuperCosmos.npy',r)
    if catalogname=='2':
    	np.save('datalog/errorbar_MCXC.npy',result)
    	np.save('datalog/errorbar_r_MCXC.npy',r)

    #print 'xi:',xi
    #print 'r:',r

    #print xi.size
    #print r.size

    nk.clear()
    rk.clear()
    cat_sim.clear_cache()

    #os.remove('data/dx12_v3_smica_nosz_cmb_mc_{}_raw.fits'.format(str(i).zfill(5)))
    #os.remove('data/dx12_v3_smica_nosz_noise_hm1_mc_{}_raw.fits'.format(str(i).zfill(5)))
    astropy.utils.data.clear_download_cache()

    end_time=time.time()
    print 'time used: ',end_time-start_time,'s'

    return result



catalogname=input('which catalog to use? 1.SuperCosmos 2.MCXC\n')
catalogname=str(catalogname)
num_of_runs=raw_input('Number of simulation runs? default=300\n')
num_of_runs=str(num_of_runs)
if num_of_runs=='':
    num_of_runs='300'
num_of_runs=int(num_of_runs)


result=np.ndarray((num_of_runs,35))

if catalogname=='1':
    cat_galaxy,cat_rand=loadGalaxy()
if catalogname=='2':
    cat_galaxy,cat_rand=loadCluster()

planckmask=fits.open('HFI_Mask_GalPlane-apo0_2048_R2.00.fits')
planckImask=planckmask[1].data['GAL080']
planckmask.close()
planckpix=np.arange(0,planckImask.size)
planckpix=planckpix[planckImask==1]
planck_ra,planck_dec=healpy.pix2ang(nside=2048,ipix=planckpix,nest=True,lonlat=True)
planckpix=None

coord=SkyCoord(planck_ra,planck_dec,frame='galactic',unit='deg').icrs
planck_ra,planck_dec=coord.ra.deg,coord.dec.deg

for i in range(num_of_runs):
    result=Calculate(cat_galaxy,cat_rand,result,num_of_runs,i,planck_ra,planck_dec)
























'''
catalog=np.load('datalog/wisecatalog.npy')
scosmask=healpy.read_map('WISExSCOSmask.fits')

num=30000000

coord=SkyCoord(catalog[0],catalog[1],frame='icrs',unit='deg').galactic
l,b=coord.l.deg,coord.b.deg

catalog=catalog[:,scosmask[healpy.ang2pix(256,l,b,nest=False,lonlat=True)]==1]

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

cat_rand = treecorr.Catalog(ra=rand_ra, dec=rand_dec, ra_units='radians', dec_units='radians')

#plt.scatter(np.rad2deg(rand_ra),np.rad2deg(rand_dec),s=0.01)
#plt.xlabel('RA(deg)')
#plt.ylabel('DEC(deg)')
#plt.show()
print('Done!\n')



num_of_runs=10
result=np.ndarray((10,35))

planckmask=pyfits.open('HFI_Mask_GalPlane-apo0_2048_R2.00.fits')
planckImask=planckmask[1].data['GAL080']
planckmask.close()
planckpix=np.arange(0,planckImask.size)
planckpix=planckpix[planckImask==1]
planck_ra,planck_dec=healpy.pix2ang(nside=2048,ipix=planckpix,nest=True,lonlat=True)
planckpix=None

coord=SkyCoord(planck_ra,planck_dec,frame='galactic',unit='deg').icrs
planck_ra,planck_dec=coord.ra.deg,coord.dec.deg

for i in range(num_of_runs):
    start_time=time.time()
    print 'run ',i,'/',num_of_runs
    print 'loading simulation and mask...'
    plancksim_n=healpy.read_map('dx12_v3_smica_nosz_cmb_mc_{}_raw.fits'.format(str(i).zfill(5)))
    plancknoise=healpy.read_map('dx12_v3_smica_nosz_noise_hm1_mc_{}_raw.fits'.format(str(i).zfill(5)))


    plancksim=plancksim_n+plancknoise
    plancksim_n=None
    plancknoise=None
    plancksim=healpy.ud_grade(plancksim,2048, order_in='RING', order_out='NEST')

    plancksim=plancksim[planckImask==1]

    cat_sim=treecorr.Catalog(ra=planck_ra,dec=planck_dec,ra_units='deg',dec_units='deg',k=plancksim)


    print 'Calculating correlation...'

    nk = treecorr.NKCorrelation(min_sep=0.01, max_sep=10, nbins=35, sep_units='deg')
    rk = treecorr.NKCorrelation(min_sep=0.01, max_sep=10, nbins=35, sep_units='deg')
    nk.process(cat_galaxy,cat_sim)
    rk.process(cat_rand,cat_sim)

    xi, varxi = nk.calculateXi(rk)


    print 'Done'

    result[i]=xi
    r = np.exp(nk.meanlogr)

    np.save('datalog/errorbar.npy',result)
    np.save('datalog/errorbar_r.npy',r)

    #print 'xi:',xi
    #print 'r:',r

    #print xi.size
    #print r.size

    nk.clear()
    rk.clear()
    cat_sim.clear_cache()
    plancksim=None
    xi=None

    end_time=time.time()
    print 'time used: ',end_time-start_time

#np.save('datalog/errorbar.npy',result)
#np.save('datalog/errorbar_r.npy',r)

'''
