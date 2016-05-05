import numpy as np
import individual_star_properties as isp
import matplotlib.pyplot as plt

from matplotlib.collections import LineCollection
from matplotlib.colors      import ListedColormap, BoundaryNorm

lw = 3.0 
normalize = True 
Zsolar = isp.const.Zsolar_parsec 
Z = 0.005
Z = Z / isp.const.Zsolar_parsec

# make a list of star masses M = np.logspace(1.0, 2.0, 100.0)
M = np.linspace(1.0, 100.0, 1000.0)

star_list_bb = [None]*np.size(M)
star_list_ostar = [None]*np.size(M)
q0_ostar = np.zeros(np.size(M))
q1_ostar = np.zeros(np.size(M))
q0_bb = np.zeros(np.size(M))
q1_bb = np.zeros(np.size(M))

flag = np.zeros(np.size(M))
N  = 1.0
for i in np.arange(np.size(M)):
    star_list_ostar[i] = isp.individual_star(M[i], Z=Z,sigma_factor =0.0, blackbody_only=False, use_ZAMS=True)
    star_list_bb[i] = isp.individual_star(M[i], Z=Z,sigma_factor=0.0, blackbody_only=True, use_ZAMS=True)

    if normalize:

        N = 4.0 * np.pi * star_list_ostar[i]._R**2

    q0_ostar[i] = star_list_ostar[i].q0 / N
    q1_ostar[i] = star_list_ostar[i].q1 / N
    q0_bb[i] = star_list_bb[i].q0 /N
    q1_bb[i] = star_list_bb[i].q1 /N

    flag[i] = star_list_ostar[i].flag


plt.plot(M, q0_ostar / q0_bb, label = 'HI - OSTAR/BB', lw = lw, color = 'black', ls = '-')
plt.plot(M, q1_ostar / q1_bb, label = 'HeI - OSTAR/BB', lw = lw, color = 'black', ls = '--')



plt.legend(loc='lower right')
plt.xlabel(r'Mass (M$_{\odot})$')
plt.minorticks_on()
plt.ylabel(r'Ionizing Photon Flux OSTAR / BB (s$^{-1}$ cm$^{-2}$)')


outname = '%.5fZ'%(Z*isp.const.Zsolar_parsec)

plt.semilogy()
plt.savefig('q_rad/ratio/' + outname + '.png')
