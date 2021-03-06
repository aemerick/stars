import numpy as np
import individual_star_properties as isp
import matplotlib.pyplot as plt

from matplotlib.collections import LineCollection
from matplotlib.colors      import ListedColormap, BoundaryNorm

lw             = 3.0
normalize      = True
use_ZAMS       = True
blackbody_only = False
color_lines    = True
light_blue     = 'blue'
light_red      = 'orange'
Zsolar         = isp.const.Zsolar_parsec
Z              = 0.001
Z = Z / isp.const.Zsolar_parsec

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
#M = np.logspace(1.0, 2.0, 100.0)
M = np.linspace(1.0, 100.0, 1000)

star_list = [None]*np.size(M)
q0 = np.zeros(np.size(M))
q1 = np.zeros(np.size(M))
q2 = np.zeros(np.size(M))
ir = np.zeros(np.size(M))
fuv = np.zeros(np.size(M))
LWf  = np.zeros(np.size(M))
E0 = np.zeros(np.size(M))
E1 = np.zeros(np.size(M))
E2 = np.zeros(np.size(M))
flag = np.zeros(np.size(M))
N  = 1.0
for i in np.arange(np.size(M)):
    star_list[i] = isp.individual_star(M[i], Z=Z,sigma_factor =0.0, blackbody_only=blackbody_only, use_ZAMS=use_ZAMS)

    if normalize:
        N = 4.0 * np.pi * star_list[i]._R**2

    q0[i] = star_list[i].q0 / N
    q1[i] = star_list[i].q1 / N
    q2[i] = star_list[i].q2 / N
    ir[i] = star_list[i].IR / N
    LWf[i] = star_list[i].LW / N
    fuv[i] =star_list[i].FUV / N
    E0[i] = star_list[i].E0
    E1[i] = star_list[i].E1
    E2[i] = star_list[i].E2
    flag[i] = star_list[i].flag


#
# generate fit points:
#
index_low  = [i for i,x in enumerate(flag) if x == 0][0]
index_high = [i for i,x in enumerate(flag) if x == 0][-1]

#
#
#
names = ['HI','HeI','HeII','IR','FUV','LW']
for i, x in enumerate([q0,q1,q2,ir,fuv,LWf]):

    if index_low > 0:
        low_factor = x[index_low]/x[index_low-1]
        M_low = M[index_low]
    else:
        low_factor = 0.0
        M_low = np.min(M)

    if index_high < np.size(x)-1:
        high_factor = x[index_high]/x[index_high+1]
        M_high = M[index_high]
    else:
        high_factor = 0.0
        M_high = np.max(M)

    print(names[i], " M_low = %.2f M_high = %.2f -- ratios %5.5E %5.5E "%(M_low,M_high,low_factor, high_factor))

if color_lines:
    lc1 = generate_colored_lines(M, q0, flag, lw = lw, color = ['C0','black'], ls = '-')

    plt.gca().add_collection(lc1)

    lc2 = generate_colored_lines(M, q1, flag, lw = lw, color = ['C1','black'], ls = '--')

    plt.gca().add_collection(lc2)

    lc3 = generate_colored_lines(M, q2, flag, lw = lw, color = ['C2','black'], ls = ':')

    plt.gca().add_collection(lc3)

else:
    plt.plot(M, q0, label = 'HI', lw = lw, color = 'black', ls = '-')
    plt.plot(M, q1, label = 'HeI', lw = lw, color = 'black', ls = '--')
    plt.plot(M, q2, lable = 'HeII', lw = lw, color ='black', ls = ':')

if not use_ZAMS:
    # repeat the above with the spread:
    star_list = [None]*np.size(M)
    q0 = np.zeros(np.size(M))
    q1 = np.zeros(np.size(M))
    q2 = np.zeros(np.size(M))
    E0 = np.zeros(np.size(M))
    E1 = np.zeros(np.size(M))
    E2 = np.zeros(np.size(M))
    flag = np.zeros(np.size(M))
    for i in np.arange(np.size(M)):
        star_list[i] = isp.individual_star(M[i], Z=Z,sigma_factor =1.0, blackbody_only=blackbody_only, use_ZAMS=use_ZAMS)
        if normalize:
            N = 4.0 * np.pi * star_list[i]._R**2

        q0[i] = star_list[i].q0 / N
        q1[i] = star_list[i].q1 / N
        q2[i] = star_list[i].q2 / N
        E0[i] = star_list[i].E0
        E1[i] = star_list[i].E1
        E2[i] = star_list[i].E2
        flag[i] = star_list[i].flag

    if color_lines:
        lc3 = generate_colored_lines(M, q0, flag, lw = 0.8*lw, color = [light_red,'grey'], ls = '-')

        plt.gca().add_collection(lc3)
        lc4 = generate_colored_lines(M, q1, flag, lw = 0.8*lw, color = [light_blue,'grey'], ls = '--')
        plt.gca().add_collection(lc4)
        lc5 = generate_colored_lines(M, q2, flag, lw = 0.8*lw, color = [light_blue,'grey'],ls=':')
        plt.gca().add_collection(lc5)

    else:
        plt.plot(M, q0, lw = 0.8*lw, color = 'grey', ls = '-')
        plt.plot(M, q1, lw = 0.8*lw, color = 'grey', ls = '--')
        plt.plot(M, q2, lw = 0.8*lw, color = 'grey', ls = ':')


    # and again
    # repeat the above with the spread:
    star_list = [None]*np.size(M)
    q0 = np.zeros(np.size(M))
    q1 = np.zeros(np.size(M))
    q2 = np.zeros(np.size(M))
    E0 = np.zeros(np.size(M))
    E1 = np.zeros(np.size(M))
    E2 = np.zeros(np.size(M))
    flag = np.zeros(np.size(M))
    for i in np.arange(np.size(M)):
        star_list[i] = isp.individual_star(M[i], Z=Z, sigma_factor = -1.0, blackbody_only=blackbody_only, use_ZAMS=use_ZAMS)
        if normalize:
            N = 4.0 * np.pi * star_list[i]._R**2


        q0[i] = star_list[i].q0 / N
        q1[i] = star_list[i].q1 / N
        q2[i] = star_list[i].q2 / N
        E0[i] = star_list[i].E0
        E1[i] = star_list[i].E1
        E2[i] = star_list[i].E2
        flag[i] = star_list[i].flag

    if color_lines:
        lc5 = generate_colored_lines(M, q0, flag, lw = 0.8*lw, color = [light_red,'grey'], ls = '-')
        plt.gca().add_collection(lc5)
        lc6 = generate_colored_lines(M, q1, flag, lw = 0.8*lw, color = [light_blue,'grey'], ls = '--')
        plt.gca().add_collection(lc6)
        lc7 = generate_colored_lines(M, q2, flag, lw = 0.8*lw, color = [light_blue,'grey'], ls = '--')
        plt.gca().add_colections(lc7)

    else:
        plt.plot(M, q0, lw = 0.8*lw, color = 'grey', ls = '-')
        plt.plot(M, q1, lw = 0.8*lw, color = 'grey', ls = '--')
        plt.plot(M, q2, lw = 0.8*lw, color = 'grey', ls = ':')
# ---------- end if spread lines -------- #

plt.semilogy()
if normalize:
    plt.ylim([1.0E12,1.0E26])
else:
    plt.ylim([1.0E42, 4.0E50])


if color_lines:
    xlim = plt.xlim()
    ylim = plt.ylim()
    plt.plot([1000.0,1001.0],[1.0,2.0], ls ='-', color = 'C0', lw = lw, label ='HI')
    plt.plot([1000.0,1001.0],[1.0,2.0], ls ='-', color = 'C1', lw = lw, label ='HeI')
    plt.plot([1000.0,1001.0],[1.0,2.0], ls ='-', color = 'C2', lw = lw, label ='HeII')
    plt.xlim(xlim); plt.ylim(ylim)


plt.legend(loc='lower right')
plt.xlabel(r'Mass (M$_{\odot})$')
plt.minorticks_on()
if normalize:
    plt.ylabel(r'Ionizing Photon Flux (s$^{-1}$ cm$^{-2}$)')
else:
    plt.ylabel(r'Ionizing Photon Rate (s$^{-1}$)')


outname = '%.5fZ'%(Z*isp.const.Zsolar_parsec)
if normalize:
    outname += '_norm'

if use_ZAMS:
    outname += '_ZAMS'

plt.savefig('q_rad/' + outname + '.png')
