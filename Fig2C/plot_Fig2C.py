import matplotlib.pyplot as plt
import numpy as np
import scipy.stats

import matplotlib
font = {'family' : 'normal',
#         'weight' : 'bold',
        'size'   : 18}

matplotlib.rc('font', **font)

def get_post_spikes(v):
    nspikes = 0
    last_v = v[0]
    for i in range(len(v)):
        current_v = v[i]
        if current_v > 0.0 and last_v <=0.0:
            nspikes+= 1
        last_v = current_v
    return nspikes

miaves1=[]
mierrs1=[]
miaves2=[]
mierrs2=[]
for ni in range(16):
    mili1 = []
    mili2 = []
    for trial in range(1,10):
        dpath = 'datag6/%s_100.0_4.0_1.0_1.0_%s'%(ni,trial)
        midat = np.loadtxt('%s/mi.txt'%dpath,delimiter=',')

        mili1.append(midat[0,0])
        mili2.append(midat[-1,0])
    miaves1.append(np.mean(mili1))
    mierrs1.append(scipy.stats.sem(mili1))
    miaves2.append(np.mean(mili2))
    mierrs2.append(scipy.stats.sem(mili2))

    
fig,ax = plt.subplots(nrows =1, ncols = 1, figsize = (5,3))
aa = [i for i in range(16)]
ax.scatter(aa,miaves1,c='steelblue', marker='o',label="$x_1$; soma")
ax.errorbar(aa,miaves1, yerr = mierrs1, color = 'steelblue', capsize = 3.5)

ax.scatter(aa,miaves2,c='steelblue', marker='^',label="$x_1$; dendrite")
ax.errorbar(aa,miaves2, yerr = mierrs2, marker='^', color = 'steelblue', capsize = 3.5)

ax.set_xticks([0,5,10,10])
ax.set_xlabel('$n$')
# ax.set_yticks([0,0.005,0.01])
ax.set_ylabel('$I(x_1; y)$')
ax.tick_params(direction='in', which='both')
plt.subplots_adjust(left=0.3, right=0.9, top=0.9, bottom=0.2)
 
plt.legend(prop={'size': 12})
plt.savefig('Fig2Cmi.eps')


faves1=[]
ferrs1=[]
faves2=[]
ferrs2=[]
for ni in range(16):
    fli1 = []
    fli2 = []
    for trial in range(1,10):
        dpath = 'datag6/%s_100.0_4.0_1.0_1.0_%s'%(ni,trial)
        vdat = np.loadtxt('%s/vars.txt'%dpath)
        dvdat = np.loadtxt('%s/dv.txt'%dpath)

        fli1.append(get_post_spikes(vdat[:,1]))
        fli2.append(get_post_spikes(dvdat[:,-2]))
    faves1.append(np.mean(fli1))
    ferrs1.append(scipy.stats.sem(fli1))
    faves2.append(np.mean(fli2))
    ferrs2.append(scipy.stats.sem(fli2))

fig,ax = plt.subplots(nrows =1, ncols = 1, figsize = (5,3))
faves1 = np.array(faves1)
faves2 = np.array(faves2)
ferrs1 = np.array(ferrs1)
ferrs2 = np.array(ferrs2)
aa = [i for i in range(16)]
ax.scatter(aa,faves1/30,c='steelblue', marker='o',label="soma")
ax.errorbar(aa,faves1/30, yerr = ferrs1/30, color = 'steelblue', capsize = 3.5)

ax.scatter(aa,faves2/30,c='steelblue', marker='^',label="dendrite")
ax.errorbar(aa,faves2/30, yerr = ferrs2/30, marker='^', color = 'steelblue', capsize = 3.5)

ax.set_xticks([0,5,10,15])
ax.set_xlabel('$n$')
ax.set_yticks([0,5,10])
ax.set_ylabel('$f$ (Hz)')

ax.tick_params(direction='in', which='both')
plt.subplots_adjust(left=0.3, right=0.9, top=0.9, bottom=0.2)
 
plt.legend(prop={'size': 12})
plt.savefig('Fig2Cf.eps')