import numpy as np
import individual_star_properties as isp
import matplotlib.pyplot as plt

from matplotlib.collections import LineCollection
from matplotlib.colors      import ListedColormap, BoundaryNorm

lw             = 3.0
normalize      = True
blackbody_only = False
use_ZAMS       = True
color_lines    = True
light_blue     = 'lightblue'
light_red      = 'lightcoral'
Z              = 0.00100001
Z = Z / isp.const.Zsolar_parsec

def generate_colored_lines(x,y,cval, color = ['red','black'], ls = '-', lw = lw):

    cmap     = ListedColormap([color[0], 'purple', color[1]])
    norm     = BoundaryNorm([-0.1, 0.1, 0.9, 1.1], cmap.N)

    points   =  np.array([x, y]).T.reshape(-1,1,2)
    segments =  np.concatenate([points[:-1], points[1:]], axis=1)



    lc = LineCollection(segments, cmap = cmap, norm = norm)
    lc.set_array(cval)
    lc.set_linewidth(lw)
    lc.set_linestyle(ls)
    

    return lc

# make a list of star masses
#M = np.logspace(1.0, 2.0, 100.0)
M = np.linspace(1.0, 100.0, 1000.0)

star_list = [None]*np.size(M)
T = np.zeros(np.size(M))
g = np.zeros(np.size(M))
L = np.zeros(np.size(M))
R = np.zeros(np.size(M))
flag = np.zeros(np.size(M))
N  = 1.0
for i in np.arange(np.size(M)):
    star_list[i] = isp.individual_star(M[i], Z=Z,sigma_factor =0.0, blackbody_only=blackbody_only, use_ZAMS=use_ZAMS)

    if normalize:
        N = 4.0 * np.pi * star_list[i]._R**2

    T[i] = star_list[i]._T 
    g[i] = star_list[i].g
    L[i] = star_list[i].L
    R[i] = star_list[i].R
    flag[i] = star_list[i].flag

fig, ax = plt.subplots(2,2)

ax[(0,0)].plot(M, T/1.0E3, lw = lw, color = 'black', ls = '-')

# plot OSTAR limits
ax[(0,0)].plot([np.min(M),np.max(M)], [57.5,57.5], color = 'black', ls = '--', lw=lw*0.8)
ax[(0,0)].plot([np.min(M),np.max(M)], [27.5,27.5], color = 'black', ls = '--', lw=lw*0.8)

ax[(0,1)].plot(M, np.log10(g), lw = lw, color = 'black', ls = '-')
ax[(0,1)].plot([np.min(M),np.max(M)], [3.0, 3.0], color = 'black', ls = '--', lw=lw*0.8)
ax[(0,1)].plot([np.min(M),np.max(M)], [3.5, 3.5], color = 'black', ls = '--', lw=lw*0.5)
ax[(0,1)].plot([np.min(M),np.max(M)], [4.0, 4.0], color = 'black', ls = '--', lw=lw*0.25)

ax[(1,0)].plot(T/1.0E3, np.log10(g), lw = lw, color = 'black', ls = '-')


# plot radius
ax[(1,1)].plot(M, R, lw = lw, color = 'black', ls = '-')


ax[(0,0)].set_xlabel(r'Mass (M$_{\odot}$)')
ax[(0,1)].set_xlabel(r'Mass (M$_{\odot}$)')
ax[(1,0)].set_xlabel(r'Temperature (kK)')
ax[(1,1)].set_xlabel(r'Mass (M$_{\odot}$)')

ax[(0,0)].set_ylabel(r'Temperature (kK)')
ax[(0,1)].set_ylabel(r'log(Surface Gravity [g cm$^{2}$ s$^{-2}$])')
ax[(1,0)].set_ylabel(r'log(Surface Gravity [g cm$^{2}$ s$^{-2}$])')
ax[(1,1)].set_ylabel(r'Radius (R$_{\odot}$)')


ax[(0,1)].set_ylim(3.0,5.0)
ax[(1,0)].set_ylim(2.9,5.1)
ax[(0,0)].set_ylim(0.0, 70)
ax[(1,0)].set_xlim(10, 65.0)
ax[(1,1)].set_ylim(0.0,15.0)
ax[(1,0)].set_xlim(ax[(1,0)].get_xlim()[::-1])
ax[(1,0)].set_ylim(ax[(1,0)].get_ylim()[::-1])

# plot ostar points:
data = np.genfromtxt('./data/OSTAR2002.txt',names=True)

ax[(1,0)].scatter( data['T'] / 1000.0, data['logg'], color = 'red')

#ax[2].set_xlabel(r'log(Surface Gravity [g cm$^{2}$ s$^{-2}$])')
#ax[2].set_ylabel(r'Temperature (K)')
for a in [(0,0),(0,1),(1,0),(1,1)]:
    ax[a].minorticks_on()
fig.set_size_inches(10,8)
plt.tight_layout()
outname = '%0.5fZ.png'%(Z*isp.const.Zsolar_parsec)
fig.savefig('star_properties/' + outname)
