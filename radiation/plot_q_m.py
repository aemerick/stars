import numpy as np
import individual_star_properties as isp
import matplotlib.pyplot as plt

lw = 3.0

normalize = True




# make a list of star masses
M = np.logspace(0.0, 2.0, 100.0)

star_list = [None]*np.size(M)
q0 = np.zeros(np.size(M))
q1 = np.zeros(np.size(M))
E0 = np.zeros(np.size(M))
E1 = np.zeros(np.size(M))
flag = np.zeros(np.size(M))
N  = 1.0
for i in np.arange(np.size(M)):
    star_list[i] = isp.individual_star(M[i])

    if normalize:
        N = 4.0 * np.pi * star_list[i]._R**2

    q0[i] = star_list[i].q0 / N
    q1[i] = star_list[i].q1 / N
    E0[i] = star_list[i].E0
    E1[i] = star_list[i].E1
    flag[i] = star_list[i].flag

plt.plot(M, q0, label='HI - Blackbody', lw = lw, color = 'black', ls = '-')
plt.plot(M, q1, label='HeI - Blackbody', lw = lw, color = 'black', ls = '--')

# repeat the above with the spread:
star_list = [None]*np.size(M)
q0 = np.zeros(np.size(M))
q1 = np.zeros(np.size(M))
E0 = np.zeros(np.size(M))
E1 = np.zeros(np.size(M))

for i in np.arange(np.size(M)):
    star_list[i] = isp.individual_star(M[i], sigma_factor = 1.0)
    if normalize:
        N = 4.0 * np.pi * star_list[i]._R**2

    q0[i] = star_list[i].q0 / N
    q1[i] = star_list[i].q1 / N
    E0[i] = star_list[i].E0
    E1[i] = star_list[i].E1
    flag[i] = star_list[i].flag

plt.plot(M, q0, lw = 0.8*lw, color = 'grey', ls = '-')
plt.plot(M, q1, lw = 0.8*lw, color = 'grey', ls = '--')


# and again
# repeat the above with the spread:
star_list = [None]*np.size(M)
q0 = np.zeros(np.size(M))
q1 = np.zeros(np.size(M))
E0 = np.zeros(np.size(M))
E1 = np.zeros(np.size(M))

for i in np.arange(np.size(M)):
    star_list[i] = isp.individual_star(M[i], sigma_factor = -1.0)
    if normalize:
        N = 4.0 * np.pi * star_list[i]._R**2


    q0[i] = star_list[i].q0 / N
    q1[i] = star_list[i].q1 / N
    E0[i] = star_list[i].E0
    E1[i] = star_list[i].E1
    flag[i] = star_list[i].flag

plt.plot(M, q0, lw = 0.8*lw, color = 'grey', ls = '-')
plt.plot(M, q1, lw = 0.8*lw, color = 'grey', ls = '--')



plt.loglog()
plt.legend(loc='best')
plt.xlabel(r'Mass (M$_{\odot}$')
plt.minorticks_on()
if normalize:
    plt.ylabel(r'Ionizing Photon Flux (s$^{-1}$ cm$^{-2}$)')
else:
    plt.ylabel(r'Ionizing Photon Rate (s$^{-1}$)')

if normalize:
    plt.savefig('q_m_norm.png')
else:
    plt.savefig('q_m.png')

print flag
