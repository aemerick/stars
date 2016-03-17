import numpy as np
import individual_star_properties as isp
import matplotlib.pyplot as plt

lw = 3.0

normalize = False

# make a list of star masses
M = np.logspace(0.0, 2.0, 1000.0)

star_list = [None]*np.size(M)
q0 = np.zeros(np.size(M))
q1 = np.zeros(np.size(M))
E0 = np.zeros(np.size(M))
E1 = np.zeros(np.size(M))
N  = 1.0
for i in np.arange(np.size(M)):
    star_list[i] = isp.individual_star(M[i])

    if normalize:
        N = 4.0 * np.pi * star_list[i]._R**2

    q0[i] = star_list[i].q0 / N
    q1[i] = star_list[i].q1 / N
    E0[i] = star_list[i].E0
    E1[i] = star_list[i].E1

plt.plot(M, E0, label='HI', lw = lw, color = 'black', ls = '-')
plt.plot(M, E1, label='HeI', lw = lw, color = 'black', ls = '--')

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

plt.plot(M, E0, lw = 0.8*lw, color = 'grey', ls = '-')
plt.plot(M, E1, lw = 0.8*lw, color = 'grey', ls = '--')


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

plt.plot(M, E0, lw = 0.8*lw, color = 'grey', ls = '-')
plt.plot(M, E1, lw = 0.8*lw, color = 'grey', ls = '--')


plt.semilogx()
xlim = plt.xlim()
ylim = plt.ylim()

plt.plot(xlim, [13.6,13.6], lw = 2, ls = '-', color ='black')
plt.plot(xlim, [24.587, 24.587], lw = 2, ls = '--', color ='black')

plt.legend(loc='best')
plt.xlabel(r'Mass (M$_{\odot}$')

plt.ylim(12.0,35.0)
plt.minorticks_on()
plt.ylabel(r'Average Photon Energy (eV)')

plt.savefig('E_m.png')
