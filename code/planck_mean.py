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

import pyfits
import treecorr
import healpy
from astropy.coordinates import SkyCoord



import numpy as np
import matplotlib.pyplot as plt


single_fre=raw_input('1.Full(default) 2.44 GHz Single Frequency test 3.143 GHz 4.217 GHz 5.545 GHz\n')
single_fre=str(single_fre)
if single_fre=='':
    single_fre='1'


catalogname=input('which catalog to use? 1.SuperCosmos 2.MCXC\n')

catalogname=str(catalogname)

if catalogname=='1':
    print('loading SuperCosmos catalog and mask...')

    #load supercosmos catalog
    #catalog=np.loadtxt('wiseScosPhotoz160708.csv',skiprows=1,delimiter=",",usecols=(7,8,11,16,18),max_rows=3000)   #RA, DEC,Ebv(extinction),z,mask from "all sky survey"
    #catalog=np.loadtxt('wiseScosPhotoz160708.csv',skiprows=1,delimiter=",",usecols=(7,8,16))   #ra,dec,z (in degrees)
    #print(catalog)

    #catalog=catalog.transpose()

    catalog=np.load('datalog/wisecatalog.npy')

    scosmask=healpy.read_map('WISExSCOSmask.fits')



    #num=catalog[0].size
    num=30000000
    #coord=SkyCoord(catalog[0],catalog[1],frame='galactic',unit='deg').icrs

    #catalog[0],catalog[1]=coord.ra.deg,coord.dec.deg
    coord=SkyCoord(catalog[0],catalog[1],frame='icrs',unit='deg').galactic
    l,b=coord.l.deg,coord.b.deg




    catalog=catalog[:,scosmask[healpy.ang2pix(256,l,b,nest=False,lonlat=True)]==1]

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



if catalogname=='2':
    print('loading MCXC catalog and mask...')

    num=1000000

    MCXCfile=pyfits.open('MCXC.fits')
    MCXCdata=MCXCfile[1].data
    MCXCfile.close()

    MCXCmaskfile=pyfits.open('HFI_Mask_GalPlane-apo0_2048_R2.00.fits')
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

    #scosmask=healpy.read_map('WISExSCOSmask.fits')
    planckmask=pyfits.open('HFI_Mask_GalPlane-apo0_2048_R2.00.fits')
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

    #plt.scatter(np.rad2deg(rand_ra),np.rad2deg(rand_dec),s=0.5)
    #plt.scatter(MCXCdata['RA'],MCXCdata['DEC'],s=0.5)
    #plt.xlabel('RA(deg)')
    #plt.ylabel('DEC(deg)')
    #plt.show()
    #input('stop here')
    print('Done!\n')

cat_rand = treecorr.Catalog(ra=rand_ra, dec=rand_dec, ra_units='radians', dec_units='radians')


#load planck data
print('loading Planck catalog and mask...')
if single_fre=='1':
    planckdata=pyfits.open('COM_CMB_IQU-smica-nosz_2048_R3.00_full.fits')
    #planckdata=pyfits.open('COM_CMB_IQU-smica_2048_R3.00_full.fits')
if single_fre=='2':
    planckdata=pyfits.open('LFI_SkyMap_044-field-IQU_1024_R3.00_full.fits')
if single_fre=='3':
    planckdata=pyfits.open('HFI_SkyMap_143-field-IQU_2048_R3.00_full.fits')
if single_fre=='4':
    planckdata=pyfits.open('HFI_SkyMap_217-field-IQU_2048_R3.00_full.fits')
if single_fre=='5':
    planckdata=pyfits.open('HFI_SkyMap_545-field-Int_2048_R3.00_full.fits')

planckmask=pyfits.open('HFI_Mask_GalPlane-apo0_2048_R2.00.fits')


planckImap=planckdata[1].data['I_STOKES']

planckImask=planckmask[1].data['GAL080']

planckdata.close()
planckmask.close()

if single_fre=='2':
    planckImask=healpy.pixelfunc.ud_grade(planckImask,1024,order_in='NEST',order_out='NEST')

#hp = HEALPix(nside=2048, order='nested', frame=ICRS()

planckpix=np.arange(0,planckImap.size)



planckImask=healpy.pixelfunc.ud_grade(scosmask,2048,order_in='RING',order_out='NEST')



planckImap=planckImap[planckImask==1]

planckpix=planckpix[planckImask==1]

planck_ra,planck_dec=healpy.pix2ang(nside=2048,ipix=planckpix,nest=True,lonlat=True)




coord=SkyCoord(planck_ra,planck_dec,frame='galactic',unit='deg').icrs
planck_ra,planck_dec=coord.ra.deg,coord.dec.deg

'''
plt.scatter(planck_ra,planck_dec,s=0.01,c=planckImap,cmap='rainbow',edgecolors='none')
plt.xlabel('RA(deg)')
plt.ylabel('DEC(deg)')
plt.title('Planck map after mask')
plt.colorbar()
plt.show()
input('stop here')
'''
print 'planck mean:',np.mean(planckImap)

raw_input('stop')

cat_planck=treecorr.Catalog(ra=planck_ra,dec=planck_dec,ra_units='deg',dec_units='deg',k=planckImap)

print('Done!\n')





#print planckpix.size,planckImap.size,plancknoise.size

print('calculating cross-relation...')
#'''cross correlation
nk = treecorr.NKCorrelation(min_sep=0.01, max_sep=10, nbins=35, sep_units='deg')
rk = treecorr.NKCorrelation(min_sep=0.01, max_sep=10, nbins=35, sep_units='deg')
#nk = treecorr.KKCorrelation(min_sep=0.01, max_sep=3.0, nbins=50, sep_units='radians')
nk.process(cat_galaxy,cat_planck)
rk.process(cat_rand,cat_planck)

xi, varxi = nk.calculateXi(rk)
sig = np.sqrt(varxi)
r = np.exp(nk.meanlogr)
#'''
'''auto correlation
nn = treecorr.NNCorrelation(min_sep=0.01, max_sep=10, nbins=35, sep_units='deg')
rr = treecorr.NNCorrelation(min_sep=0.01, max_sep=10, nbins=35, sep_units='deg')
nr = treecorr.NNCorrelation(min_sep=0.01, max_sep=10, nbins=35, sep_units='deg')
nn.process(cat_galaxy)
rr.process(cat_rand)
nr.process(cat_galaxy, cat_rand)
xi, varxi = nn.calculateXi(rr, nr)
sig = np.sqrt(varxi)
r = np.exp(nn.meanlogr)
'''
print xi


print('Done!\n')


print('Plotting')

#not sure what r is

'''auto correlation plot
plt.plot(r, xi, color='blue')
plt.plot(r,np.zeros(xi.size),color='red')

plt.errorbar(r, xi, yerr=sig, lw=1, ls='',ecolor='g')

#plt.xlabel(r'$\theta$ (rad)')
plt.xlabel(r'$\theta$ (degrees)')
plt.xscale('log')
if catalogname=='1':
    plt.title('SuperCosmos x SuperCosmos')
if catalogname=='2':
    plt.title('MCXC x MCXC')
plt.show()
'''

#'''cross correlation plot
plt.plot(r, xi, color='blue')
#plt.plot(r, -xi, color='blue', ls=':')
plt.errorbar(r[xi>0], xi[xi>0], yerr=sig[xi>0], lw=1, ls='',ecolor='g')
#plt.errorbar(r[xi<0], -xi[xi<0], yerr=sig[xi<0], lw=1, ls='',ecolor='g')
leg = plt.errorbar(-r, xi, yerr=sig, color='blue')

plt.xscale('log')
#plt.yscale('log', nonposy='clip')
plt.xlabel(r'$\theta$ (degrees)')
#plt.xlabel(r'$\theta$ (rad)')
plt.ylabel(r'$w(\theta)$')
plt.legend([leg], [r'$w(\theta)$'], loc='lower left')
if catalogname=='1':
    if single_fre=='1':
        plt.title('SuperCosmos x Planck')
    if single_fre=='2':
        plt.title('SuperCosmos x Planck 44GHz')
    if single_fre=='3':
        plt.title('SuperCosmos x Planck 143GHz')
    if single_fre=='4':
        plt.title('SuperCosmos x Planck 217GHz')
    if single_fre=='5':
        plt.title('SuperCosmos x Planck 545GHz')
if catalogname=='2':
    if single_fre=='1':
        plt.title('MCXC x Planck')
    if single_fre=='2':
        plt.title('MCXC x Planck 44GHz')
    if single_fre=='3':
        plt.title('MCXC x Planck 143GHz')
    if single_fre=='4':
        plt.title('MCXC x Planck 217GHz')
    if single_fre=='5':
        plt.title('MCXC x Planck 545GHz')
plt.show()
#'''
'''save data
if single_fre=='1':
    if catalogname=='1':
        np.save('datalog/SuperCosmos.npy',np.array([xi,r,sig]))
        print('datalog saved!')
    if catalogname=='2':
        np.save('datalog/MCXC.npy',np.array([xi,r,sig]))
        print('datalog saved!')

'''
