# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 00:08:07 2019

@author: zehaojin
"""



import numpy as np
import matplotlib.pyplot as plt


data=np.load('datalog/SuperCosmos.npy')
data_x=data[1]
data_y=data[0]



theory_y=np.load('w_theta.npy')
theory_x=np.load('datalog/errorbar_r_SuperCosmos.npy')


print 'theory_x:'
print theory_x
print 'theory_y'
print theory_y

#'''
plt.plot(data_x,data_y,label='data')
plt.plot(theory_x,-theory_y,label='-theory')
plt.xscale('log')
#plt.yscale('log')
plt.xlabel(r'$\theta$ (degrees)')
plt.ylabel(r'$w(\theta) (K)$')
#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.title('CMB x Galaxy')
plt.legend()

plt.show()
#'''
