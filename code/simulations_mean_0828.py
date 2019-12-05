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



def Calculate(planckImask,num_of_runs,i):
    start_time=time.time()
    print 'run ',i,'/',num_of_runs
    print 'downloading planck sim map...'
    #download

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

    mean=np.mean(plancksim)

    print 'mean: ',mean

    astropy.utils.data.clear_download_cache()

    end_time=time.time()
    print 'time used: ',end_time-start_time,'s'

    return mean





num_of_runs=raw_input('Number of simulation runs? default=300\n')
num_of_runs=str(num_of_runs)
if num_of_runs=='':
    num_of_runs='300'
num_of_runs=int(num_of_runs)


result=np.zeros(num_of_runs)

scosmask=healpy.read_map('WISExSCOSmask.fits')
planckImask=healpy.pixelfunc.ud_grade(scosmask,2048,order_in='RING',order_out='NEST')
scosmask=None


for i in range(num_of_runs):
    mean=Calculate(planckImask,num_of_runs,i)
    result[i]=mean

plt.hist(result,50)
print 'overall mean: ',np.mean(result)
























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
