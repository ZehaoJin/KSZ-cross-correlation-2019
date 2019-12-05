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



if catalogname=='1':
    data=np.load('datalog/SuperCosmos.npy')
    simulation=np.load('datalog/errorbar_SuperCosmos.npy')
    #simulation_neg=np.load('datalog/errorbar_neg.npy')
    r=np.load('datalog/errorbar_r_SuperCosmos.npy')
if catalogname=='2':
    data=np.load('datalog/MCXC.npy')
    simulation=np.load('datalog/errorbar_MCXC.npy')
    #simulation_neg=np.load('datalog/errorbar_neg.npy')
    r=np.load('datalog/errorbar_r_MCXC.npy')

num_of_runs,num_of_bins=simulation.shape

###
#num_of_runs=20
###
theory_y=np.load('w_theta.npy')


r0=np.append([0],r)
datar0=np.append([0],data[1])


#print data[1],r,(r[:20]-r0[:20]),(r[20:24]-r0[20:24]),(data[1,:20]-datar0[:20]),(data[1,20:24]-datar0[20:24])
print r

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
theory_signal=np.average(theory_y[:15],weights=(r[:15]-r0[:15]))-np.average(theory_y[15:24],weights=(r[15:24]-r0[15:24]))


#data_signal*=(18672608.0/(4*np.pi))
#sim_signal*=(18672608.0/(4*np.pi))


mean=np.mean(sim_signal)
std=np.std(sim_signal)
print 'sim mean:',np.mean(sim_signal)
print 'sim std:',np.std(sim_signal)

'''
con_coef=0.99
alpha=1.-con_coef
q=1-alpha/2
z_critical=norm.ppf(q)
standard_error=std/np.sqrt(num_of_runs)
CI_lower=mean-z_critical*standard_error
CI_upper=mean+z_critical*standard_error
print q,CI_lower,CI_upper
plt.plot([-1,0,1],[CI_lower,CI_lower,CI_lower],c='red',ls=':')
plt.plot([-1,0,1],[CI_upper,CI_upper,CI_upper],c='red',ls=':')
#68.27% in 1 std
#95.45% in 2 std
'''

sorted=np.sort(sim_signal)
#print sorted[int((1-0.68)/2*300)]  #?1 std
#print sorted[int((1-0.95)/2*300)]  #?2 std



print 'data: ',data_signal
print 'theory: ',theory_signal
plt.scatter(xzeros,sim_signal,c='blue',s=3,label=str(num_of_runs)+'simulations')
plt.scatter(0,data_signal,color='red',marker='*',s=64,label='real data')
plt.scatter(0,theory_signal,color='black',marker='*',s=64,label='theory')
plt.plot([-1,0,1],[0,0,0],c='green',ls=':')
plt.errorbar(0,mean,yerr=std*2,capsize=6,ecolor='green',label='2sigma')
plt.errorbar(0,mean,yerr=std,capsize=6,ecolor='red',label='1sigma')
#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))


if catalogname=='1':
    plt.title('mean(0~10 arcmin/0.167deg)-mean(0.167 deg-1 deg) for SuperCosmos')
    plt.ylim(-2e-7,2e-7)
    plt.ylabel('(K)')
if catalogname=='2':
    plt.title('mean(0~10 arcmin/0.167deg)-mean(0.167 deg-1 deg) for MCXC')
    plt.ylim(-1e-5,1e-5)
plt.legend(loc=0)
plt.show()


'''
for i in range(num_of_runs):
    plt.plot(r,simulation[i],label='sim'+str(i))
    #plt.plot(r,simulation_neg[i],ls=':',label='-sim'+str(i))
#for i in range(num_of_runs):
    #plt.plot(r,simulation[i],label='sim'+str(i))
    #plt.plot(r,simulation_neg[i],ls=':',label='-sim'+str(i))

plt.plot(data[1],data[0], color='red',label='real data',lw=3)
plt.xscale('log')
#plt.legend(loc=0,fancybox=True)
plt.title('Planck_smica_nosz simulation 0~{}'.format(str(num_of_runs-1)))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xlabel(r'$\theta$ (degrees)')
plt.ylabel(r'$w(\theta)$')
plt.show()
'''
'''
plt.hist(sim_signal,25)
plt.plot(data_signal,10,c='r',marker='*')
plt.show()
'''







##################
'''
total=simulation+simulation_neg
for i in range(num_of_runs):
    plt.plot(r,total[i],label='sim'+str(i))
plt.show()
'''




'''
for i in range(num_of_runs):
    plt.plot(r,simulation[i],label='sim'+str(i))

plt.plot(data[1],data[0], color='red',label='real data',lw=3)
plt.xscale('log')
#plt.yscale('log', nonposy='clip')
#plt.legend(loc=0,fancybox=True)
plt.title('Planck_smica_nosz simulation 0~{}'.format(str(num_of_runs-1)))
plt.xlabel(r'$\theta$ (degrees)')
plt.show()
'''



'''
simulation_abs=np.abs(simulation)

errorbar1=np.mean(simulation,axis=0)

errorbar2=np.mean(simulation_abs,axis=0)

plt.plot(data[1],errorbar1,label='direct avg')
plt.plot(data[1],errorbar2,label='abs avg')
plt.title('size of errorbar')
plt.legend(loc=0,fancybox=True,framealpha=0.5)

plt.show()

sig=np.abs(errorbar1)

plt.plot(data[1],data[0], color='blue')
plt.plot(data[1],-data[0], color='blue', ls=':')
plt.errorbar(data[1][data[0]>0], data[0][data[0]>0], yerr=sig[data[0]>0], lw=1, ls='',ecolor='g')
plt.errorbar(data[1][data[0]<0], -data[1][data[0]<0], yerr=sig[data[0]<0], lw=1, ls='',ecolor='g')
leg = plt.errorbar(-data[1], data[0], yerr=sig, color='blue')

plt.xscale('log')
plt.yscale('log', nonposy='clip')
plt.xlabel(r'$\theta$ (degrees)')

plt.legend([leg], [r'$w(\theta)$'], loc='lower left')
plt.xlim([0.01,10])
plt.title('w(theta) with error bar')
plt.show()
'''
