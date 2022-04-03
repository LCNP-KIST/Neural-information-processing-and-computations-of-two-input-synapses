import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib
import scipy.optimize

font = {'family' : 'normal',
#         'weight' : 'bold',
        'size'   : 16}

matplotlib.rc('font', **font)


def ad_sigmoid(x,a1,a2,b1,b2,c,d):
    return d*1.0/(1.0+np.exp(-a1*(x-a2))) *  (1.0/(1.0+np.exp(b1*(x-b2)))*(1.0-c/d) + c/d )

def sigmoid(x,a1,a2,d):
    return d*1.0/(1.0+np.exp(-a1*(x-a2)))


alpha = 0.1

nli = [i for i in range(0,16)]

npo = 16

mi1ar = np.zeros(npo*npo).reshape(npo,npo)
mi2ar = np.zeros(npo*npo).reshape(npo,npo)

for q1i, q1 in enumerate(nli):
    for q2i, q2 in enumerate(nli):
        mis1 = []
        mis2 = []
        for trial in range(5):
            midat = np.loadtxt("data_Fig3C/%s_%s_%s_%s/mi.txt"%(q1,q2,alpha,trial), delimiter=',')
            mis1.append(midat[2,0])
            mis2.append(midat[2,1])
    
        mi1ar[q2i,q1i] = np.mean(mis1)
        mi2ar[q2i,q1i] = np.mean(mis2)

mitot = mi1ar + mi2ar



# Fit type 1 to sigmoid function
acti1 = []
for ij in range(npo):
    acti1.append(mitot[0,ij])
    
xx = np.arange(npo)
res = scipy.optimize.curve_fit(sigmoid,xx,acti1,p0=[1,0.8,0.15])

# Plot type 1
a1,a2,d =  res[0]
xar = np.linspace(-1,15,200)

fig, ax = plt.subplots(1,1,figsize=(6,4))

ax.plot(xar,sigmoid(xar,a1,a2,d),c='black')
ax.scatter(xx,acti1,c=(146/255,206/255,80/255))
ax.set_xlim(-1,15)

ax.set_xlabel('$n_1$')
ax.set_ylabel('$I(x_1;y) + I(x_2;y)$')

plt.subplots_adjust(left=0.2,right=0.9,top=0.9,bottom=0.2)

plt.savefig("Fig3C1.eps")


# Fit type 2
acti2 = []
for ij in range(npo):
    acti2.append(mitot[ij,ij])
    
res = scipy.optimize.curve_fit(ad_sigmoid,xx,acti2,p0=[1,0.8,0.5,10,0.05,0.15])

# Plot type 2
a1,a2,b1,b2,c,d =  res[0]
ax.scatter(xx,acti2,c=acti2,cmap=cm.coolwarm)
xar = np.linspace(-2,15,200)

fig, ax = plt.subplots(1,1,figsize=(6,4))

# ax.plot(xar,sigmoid(xar,a1,a2,d),c='black')
ax.plot(xar,ad_sigmoid(xar,a1,a2,b1,b2,c,d),c=(234/255,79/255,240/255))
ax.plot(xar,ad_sigmoid(xar,a1,a2,b1,b2,c,d),c='black')
# ax.scatter(xx,acti,c=acti,cmap=cm.coolwarm)
ax.scatter(xx,acti2,c=(234/255,79/255,240/255))
ax.set_xticks([0,5,10,15])
# ax.set_xticklabels(['0','10','20','30'])

ax.set_xlabel('$n_1 (=n_2)$ ')
ax.set_ylabel('$I(x_1;y) + I(x_2;y)$')

plt.subplots_adjust(left=0.2,right=0.9,top=0.9,bottom=0.2)

plt.savefig("Fig3C2.eps")
