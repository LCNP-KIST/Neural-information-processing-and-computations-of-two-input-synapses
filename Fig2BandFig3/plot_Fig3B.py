import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from matplotlib import cm
import matplotlib

font = {'family' : 'normal',
#         'weight' : 'bold',
        'size'   : 16}

matplotlib.rc('font', **font)

alpha = 0.1

nli = [i for i in range(0,16)]

npo = 16

#retrieve data as arrays
mi1ar = np.zeros(npo*npo).reshape(npo,npo) # I(x1;y)
mi2ar = np.zeros(npo*npo).reshape(npo,npo) # I(x2;y)

for q1i, q1 in enumerate(nli):
    for q2i, q2 in enumerate(nli):
        mis1 = []
        mis2 = []
        for trial in range(5):
            midat = np.loadtxt("data_Fig3/%s_%s_%s_%s/mi.txt"%(q1,q2,alpha,gamma,trial), delimiter=',')
            mis1.append(midat[2,0])  # Get MI from last 10 seconds
            mis2.append(midat[2,1])
    
        mi1ar[q2i,q1i] = np.mean(mis1)
        mi2ar[q2i,q1i] = np.mean(mis2)

mitot = mi1ar + mi2ar  # I(x1;y) + I(x2;y)

# All coordinates
X=[]
Y=[]
Z=[]

for q1i in range(npo):
    for q2i in range(npo):
        Z.append(mitot[q1i,q2i])
        X.append(nli[q1i])
        Y.append(nli[q2i])
    
# Type 1 coordinates
X1=[]
Y1=[]
Z1=[]

for q1i in range(npo):
    Z1.append(mitot[q1i,0])
    X1.append(nli[q1i])
    Y1.append(nli[0])

# Type 2 coordinates
X2=[]
Y2=[]
Z2=[]

for q1i in range(npo):
    Z2.append(mitot[q1i,q1i]+0.001)
    X2.append(nli[q1i])
    Y2.append(nli[q1i])

# Interpolate over data at increments of 3
Xint=[]
Yint=[]
Zint=[]

for q1i in range(0,npo,3):
    for q2i in range(0,npo,3):
        Zint.append(mitot[q1i,q2i])
        Xint.append(nli[q1i])
        Yint.append(nli[q2i])

f = interpolate.interp2d(Xint, Yint, Zint, kind = 'cubic')


#plot 3d interpolation and data points

xs = np.linspace(0,15,100)
ys = np.linspace(0,15,100)
zs = f(xs,ys) - 0.01

xs, ys = np.meshgrid(xs, ys)

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

ax.view_init(27, -110)

ax.grid(False)
ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))

surf = ax.plot_surface(xs, ys, zs, cmap=cm.coolwarm,
                       linewidth=0, zorder=-1, antialiased=False)

ax.set_xlabel("$n_1$",labelpad=10)
ax.set_ylabel("$n_2$",labelpad=10)
ax.set_zlabel("$I(x_1;y) + I(x_2;y)$", labelpad=20)
for tick in ax.zaxis.get_major_ticks():
    tick.set_pad(10)

ax.set_zlim(0,0.16)
ax.set_xlim(0,16)
ax.set_ylim(0,16)
    
ax.scatter(X, Y, Z, c=Z,cmap=cm.coolwarm)
ax.scatter(X1, Y1, Z1, c=(146/255,206/255,80/255),zorder=99)  # Type 1 
ax.scatter(X2, Y2, Z2, c=(234/255,79/255,240/255),zorder=99)  # Type 2

ax.set_axis_off()

ew = 0.6
ax.plot([0,15],[0,0],[0,0],c='black',linewidth=ew)
ax.plot([0,0],[0,15],[0,0],c='black',linewidth=ew)
ax.plot([0,15],[15,15],[0,0],c='black',linewidth=ew)
ax.plot([15,15],[0,15],[0,0],c='black',linewidth=ew)
ax.plot([0,15],[0,0],[0.16,0.16],c='black',linewidth=ew,zorder=100)
ax.plot([0,0],[0,15],[0.16,0.16],c='black',linewidth=ew,zorder=100)
ax.plot([15,15],[0,15],[0.16,0.16],c='black',linewidth=ew)
ax.plot([0,15],[15,15],[0.16,0.16],c='black',linewidth=ew)
ax.plot([0,0],[0,0],[0,0.16], c='black', linewidth=ew, zorder=100)
ax.plot([15,15],[0,0],[0,0.16], c='black', linewidth=ew, zorder=100)
ax.plot([0,0],[15,15],[0,0.16], c='black', linewidth=ew, zorder=100)
ax.plot([15,15],[15,15],[0,0.16], c='black', linewidth=ew)

# ticks
ax.plot([5,5],[0,0.5],[0,0],c='black',linewidth=ew)
ax.plot([10,10],[0,0.5],[0,0],c='black',linewidth=ew)
ax.plot([0,0.5],[5,5],[0,0],c='black',linewidth=ew)
ax.plot([0,0.5],[10,10],[0,0],c='black',linewidth=ew)
#zticks
ax.plot([15,15-0.5],[0,0],[0.08,0.08],c='black',linewidth=ew)

plt.savefig("Fig3B.png", dpi=250)
