# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 00:08:07 2019

@author: zehaojin
"""


import numpy as np
import matplotlib.pyplot as plt

'''
ISWfile=np.load('datalog/ISW.npy')
CIBfile=np.load('datalog/CIB.npy')

ISWdata = ISWfile[0]
ISWr    = ISWfile[1]
ISWr0   = np.append([0.01],ISWr)



CIBdata = CIBfile[0]*69.0*1.0e-6
CIBr    = CIBfile[1]
CIBr0   = np.append([0.01],CIBr)



datafile=np.load('datalog/SuperCosmos.npy')
data=datafile[0]
datar=datafile[1]

print ISWr
print CIBr
print datar

final=data-ISWdata-CIBdata

final=np.average(final[:15],weights=(CIBr[:15]-CIBr0[:15]))-np.average(final[15:24],weights=(CIBr[15:24]-CIBr0[15:24]))

print final
'''
#'''
for z in range(4):
    ISWfile=np.load('datalog/ISW_z_{}.npy'.format(str(z)))
    CIBfile=np.load('datalog/CIB_z_{}.npy'.format(str(z)))

    ISWdata = ISWfile[0]
    ISWr    = ISWfile[1]
    ISWr0   = np.append([0.01],ISWr)

    #print ISWdata
    #print ISWr[:15]-ISWr0[:15]

    CIBdata = CIBfile[0]
    CIBr    = CIBfile[1]
    CIBr0   = np.append([0.01],CIBr)

    theory_x=np.load('theta_deg_z_{}.npy'.format(str(z)))
    theory_y=np.load('w_theta_z_{}.npy'.format(str(z)))

    data=np.load('datalog/SuperCosmos_z_{}.npy'.format(str(z)))

    plt.plot(data[1],data[0], color='blue',label='galaxy X CMB')
    plt.plot(ISWr,ISWdata, color='black',label='galaxy X ISW')
    plt.plot(CIBr,CIBdata, color='green',label='galaxy X CIB')
    plt.plot(theory_x,theory_y,color='red',label='theory galaxy X CMB')
    plt.plot(data[1],data[0]-ISWdata-CIBdata,color='grey',label='ISW and CIB subtracted')
    plt.plot([0,1,10],[0,0,0],color='red',ls=':')
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

    z_value=np.array([(0.14293+0.0100002)/2,(0.202978+0.14293)/2,(0.25900501+0.202978)/2,(0.58010101+0.25900501)/2])
    plt.title('z around {0:.5}'.format(str(z_value[z])))

    plt.xscale('log')
    plt.xlabel(r'$\theta$ (degrees)')
    plt.ylabel(r'$w(\theta)$(K)')
    plt.legend()
    plt.show()
#'''
