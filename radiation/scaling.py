import numpy as np
import matplotlib.pyplot as plt

lw = 3

# M
M = np.linspace(10.0, 100.0, 1000)

T = M**(0.2868)
g = M**(-0.5789)

# normalize to arbitrary units
T = T / np.max(T)
g = g / np.max(g)

plt.plot(M, T, label = 'Temperature', color = "black", ls = '-',lw=lw)
plt.plot(M, g, label = 'Surface Gravity', color = 'black', ls = '--',lw=lw)

plt.legend(loc='best')
plt.xlabel(r'Mass (M$_{\odot}$)')
plt.ylabel(r'Arbitrary units')
plt.savefig('T_g_scaling.png')

