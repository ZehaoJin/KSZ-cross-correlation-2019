# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 22:16:48 2019

@author: zehaojin
"""

import pyfits
import treecorr
import healpy

import numpy as np
import matplotlib.pyplot as plt
from collections import deque


from astropy_healpix import HEALPix
from astropy.coordinates import ICRS
from astropy.coordinates import SkyCoord
from astropy import units as u

import time


'''
MCXCfile=pyfits.open('MCXC.fits')

MCXCdata=MCXCfile[1].data
'''


planckdata=pyfits.open('COM_CMB_IQU-smica-nosz_2048_R3.00_full.fits')
planckmask=pyfits.open('HFI_Mask_GalPlane-apo0_2048_R2.00.fits')


healpy.mollview(planckmask[1].data['GAL020'],nest=True,title='GAL20')
healpy.mollview(planckmask[1].data['GAL080'],nest=True,title='GAL80')


healpy.mollview(planckdata[1].data['I_STOKES'],nest=True,title='I_STOKES')


'''

planckImap=planckdata[1].data['I_STOKES']

planckImask=planckmask[1].data['GAL080']

planckdata.close()
planckmask.close()

planckImap=np.array(planckImap)
planckImask=np.array(planckImask)


hp = HEALPix(nside=2048, order='nested', frame=ICRS())


##laggy to run!!!
#cat_supercosmos=treecorr.Catalog('wiseScosPhotoz160708.csv',delimiter=',',ra_col=7,dec_col=8,ra_units='deg',dec_units='deg',k_col=16)


planck_ra=deque()
planck_dec=deque()
planck_k=deque()

start_time = time.time()

for i in range(10000):
#for i in range(200000):
    
    if planckImask[i]==1:
        coord=hp.healpix_to_skycoord(i)
        
        planck_ra.append(coord.ra.deg)
        planck_dec.append(coord.dec.deg)
        planck_k.append(planckImap[i])
        


planck_ra=np.array(planck_ra)
planck_dec=np.array(planck_dec)       
planck_k=np.array(planck_k) 

print(time.time()-start_time)
'''
