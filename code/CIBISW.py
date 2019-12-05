# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 00:08:07 2019

@author: zehaojin
"""


import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

'''
ISWfile=np.load('datalog/ISW.npy')
CIBfile=np.load('datalog/CIB.npy')

ISWdata = ISWfile[0]
ISWr    = ISWfile[1]
ISWr0   = np.append([0.01],ISWr)

print ISWdata
print ISWr[:15]-ISWr0[:15]

CIBdata = CIBfile[0]
CIBr    = CIBfile[1]
CIBr0   = np.append([0.01],CIBr)

ISW_signal=np.average(ISWdata[:15],weights=(ISWr[:15]-ISWr0[:15]))-np.average(ISWdata[15:24],weights=(ISWr[15:24]-ISWr0[15:24]))

CIB_signal=np.average(CIBdata[:15],weights=(CIBr[:15]-CIBr0[:15]))-np.average(CIBdata[15:24],weights=(CIBr[15:24]-CIBr0[15:24]))

print 'ISW:',ISW_signal
print 'CIB:',CIB_signal
'''


import numpy as np
import matplotlib.pyplot as plt


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

    ISW_signal=np.average(ISWdata[:15],weights=(ISWr[:15]-ISWr0[:15]))-np.average(ISWdata[15:24],weights=(ISWr[15:24]-ISWr0[15:24]))

    CIB_signal=np.average(CIBdata[:15],weights=(CIBr[:15]-CIBr0[:15]))-np.average(CIBdata[15:24],weights=(CIBr[15:24]-CIBr0[15:24]))


    #print 'ISW:',ISW_signal
    #print 'CIB:',CIB_signal

    data=np.load('datalog/SuperCosmos_z_{}.npy'.format(str(z)))
    simulation=np.load('datalog/sim_SuperCosmos_{}.npy'.format(str(z)))
    #simulation_neg=np.load('datalog/errorbar_neg.npy')
    r=np.load('datalog/sim_r_SuperCosmos_{}.npy'.format(str(z)))

    num_of_runs,num_of_bins=simulation.shape

    theory_y=np.load('w_theta_z_{}.npy'.format(str(z)))

    r0=np.append([0.01],r)
    datar0=np.append([0.01],data[1])
    theoryr0=np.append([0],r)

    #print data[1],r,(r[:20]-r0[:20]),(r[20:24]-r0[20:24]),(data[1,:20]-datar0[:20]),(data[1,20:24]-datar0[20:24])

    sim_signal=np.zeros(num_of_runs)
    xzeros=np.zeros(num_of_runs)

    '''
    for i in range(num_of_runs):
        sim_signal[i]=np.average(simulation[i,:20],weights=(r[:20]-r0[:20]))-np.average(simulation[i,20:24],weights=(r[20:24]-r0[20:24]))

    data_signal=np.average(data[0,:20],weights=(data[1,:20]-datar0[:20]))-np.average(data[0,20:24],weights=(data[1,20:24]-datar0[20:24]))
    '''
    for i in range(num_of_runs):
        sim_signal[i]=np.average(simulation[i,:15],weights=(r[:15]-r0[:15]))-np.average(simulation[i,15:24],weights=(r[15:24]-r0[15:24]))

    data_signal=np.average(data[0,:15],weights=(data[1,:15]-datar0[:15]))-np.average(data[0,15:24],weights=(data[1,15:24]-datar0[15:24]))
    theory_signal=np.average(theory_y[:15],weights=(r[:15]-theoryr0[:15]))-np.average(theory_y[15:24],weights=(r[15:24]-theoryr0[15:24]))


    convert_to_deltaH=3.0/theory_signal

    theory_signal*=convert_to_deltaH
    sim_signal*=convert_to_deltaH
    data_signal*=convert_to_deltaH
    CIB_signal*=convert_to_deltaH
    ISW_signal*=convert_to_deltaH


    mean=np.mean(sim_signal)
    std=np.std(sim_signal)



    z_value=np.array([(0.14293+0.0100002)/2,(0.202978+0.14293)/2,(0.25900501+0.202978)/2,(0.58010101+0.25900501)/2])
    #print sim_signal
    plt.scatter(xzeros+z_value[z],sim_signal,c='blue',s=3,label=str(num_of_runs)+'simulations')
    plt.scatter(z_value[z],data_signal,color='black',marker='*',s=64,label='real data')
    #####plt.scatter(z_value[z],theory_signal,color='black',marker='*',s=64,label='theory')

    #plt.errorbar(z_value[z],mean,yerr=std*2,capsize=6,ecolor='green',label='2sigma')
    #plt.errorbar(z_value[z],mean,yerr=std,capsize=6,ecolor='red',label='1sigma')

    #####plt.errorbar(z_value[z],data_signal,yerr=std*2,capsize=6,ecolor='green',label='2sigma')
    #####plt.errorbar(z_value[z],data_signal,yerr=std,capsize=6,ecolor='red',label='1sigma')
    plt.scatter(z_value[z],CIB_signal,color='blue',marker='*',s=128,label='CIB')
    plt.scatter(z_value[z],ISW_signal,color='grey',marker='*',s=128,label='ISW')

    final_signal=data_signal-CIB_signal-ISW_signal
    plt.scatter(z_value[z],final_signal,color='red',marker='*',s=128,label='data-CIB-ISW')
    plt.errorbar(z_value[z],final_signal,yerr=std*2,capsize=6,ecolor='green',label='2sigma')
    plt.errorbar(z_value[z],final_signal,yerr=std,capsize=6,ecolor='red',label='1sigma')


    if z==0:
        plt.legend(loc=0)


plt.plot([0,0.2,0.4,0.6],[0,0,0,0],c='green',ls=':')
plt.title('mean(0~10 arcmin/0.167deg)-mean(0.167 deg-1 deg) for SuperCosmos')
plt.xlabel('z')
plt.ylabel(r'$\Delta H$')
#plt.ylim(-2e-7,2e-7)

#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

plt.show()
