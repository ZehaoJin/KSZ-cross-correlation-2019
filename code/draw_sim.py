# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 01:41:54 2019

@author: zehaojin
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm


catalogname=input('Plot for 1.SuperCosmos 2.MCXC\n')
catalogname=str(catalogname)

for z in range(4):
    data=np.load('datalog/SuperCosmos_z_{}.npy'.format(str(z)))
    simulation=np.load('datalog/sim_SuperCosmos_{}.npy'.format(str(z)))
    #simulation_neg=np.load('datalog/errorbar_neg.npy')
    r=np.load('datalog/sim_r_SuperCosmos_{}.npy'.format(str(z)))

    num_of_runs,num_of_bins=simulation.shape

    for i in range(num_of_runs):
        plt.plot(r,simulation[i],label='sim'+str(i))

    plt.plot(data[1],data[0], color='red',label='real data',lw=3)
    plt.xscale('log')
    #plt.legend(loc=0,fancybox=True)
    plt.title('Planck_smica_nosz simulation 0~{},z_bin={}'.format(str(num_of_runs-1),str(z)))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.xlabel(r'$\theta$ (degrees)')
    plt.ylabel(r'$w(\theta)$')
    if z!=0:
        plt.ylim(-0.0000015,0.0000015)
    plt.show()
