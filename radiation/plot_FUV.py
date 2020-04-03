import numpy as np
import individual_star_properties as isp
import matplotlib.pyplot as plt

from astropy import units as u
from astropy import constants as const

from matplotlib.collections import LineCollection
from matplotlib.colors      import ListedColormap, BoundaryNorm

lw             = 3.0
normalize      = True
use_ZAMS       = True
blackbody_only = False
light_blue     = 'blue'
light_red      = 'orange'
Zsolar         = isp.const.Zsolar_parsec
Z              = 0.01
Z = Z / isp.const.Zsolar_parsec


#clight = const.c.to(u.cm/u.s).value
#distance = (100.0 * u.pc).to(u.cm).value
#n_H      = 1.0

def generate_colored_lines(x,y,cval, color = ['orange','black'], ls = '-', lw = lw):

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
#M = np.logspace(1.0, 2.0, 100)
M = np.linspace(1.0, 100.0, 1000)

star_list = [None]*np.size(M)
T   = np.zeros(np.size(M))
R   = np.zeros(np.size(M))
FUV = np.zeros(np.size(M))
IR = np.zeros(np.size(M))
LWf = np.zeros(np.size(M))
flag = np.zeros(np.size(M))
N  = 1.0
for i in np.arange(np.size(M)):
    star_list[i] = isp.individual_star( M[i], Z = Z, sigma_factor = 0.0, blackbody_only=blackbody_only, use_ZAMS=use_ZAMS)
    T[i]    = star_list[i]._T
    R[i]    = star_list[i]._R
    FUV[i]  = star_list[i].FUV
    IR[i]   = star_list[i].IR
    LWf[i]  = star_list[i].LW
    flag[i] = star_list[i].flag

N = 4.0 * np.pi * R[i]**2
plt.plot(M, FUV * N , lw = lw, color = 'C0', ls = '-', label = 'FUV')
plt.plot(M, IR  * N , lw = lw, color = 'C1', ls = '-', label = 'IR')
plt.plot(M, LWf  * N , lw = lw, color = 'C2', ls = '-', label = 'LW')

plt.legend(loc='lower right')
plt.xlabel(r'Mass (M$_{\odot})$')
plt.minorticks_on()

plt.ylabel(r'FUV Luminosity (erg s$^{-1}$)')

outname = '%.5fZ'%(Z*isp.const.Zsolar_parsec)


plt.semilogy()
plt.savefig('FUV_rad/' + outname + '.png')
