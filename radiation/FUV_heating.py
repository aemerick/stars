"""
  Using OSTAR2002 SED's and CLOUDY dust models, compute the 
  normalized heating rate for a given star.

  Heating rate computed is:

      E_heat / (n_H**2 * c)

  where E_heat has units of (erg s^-1 cm^-3), n_H is hydrogen
  number density and c is the speed of light.

  Intent is to take this value as the heating rate density
  "at the star" and reduce this by 1/r**2 outward from the 
  stellar source.
"""

import numpy as np
from astropy import units as u
from astropy import constants as const
from scipy import integrate
from scipy import interpolate

micron    = (1.0 * u.micrometer).to(u.cm).value
Ryd       = (const.h * const.c * const.Ryd).to(u.erg).value

#szd_file   = './data/graphite_szd_cloudy.dat'
#ostar_file = './data/ostar2002_p03_extracted.dat'
#abs_file   = './data/graphite_ism_01.opc'

def integrate_szd(a, dNda, units='CLOUDY'):

    # if units == CLOUDY assume cloudy defaults
    # else assume the user knows what they are doing

    # convert to dNda from a^4 dNda and micron to cm
    if units == 'CLOUDY':
        a = a * micron
        dNda = dNda / a**4

    # now integrate from a min to a max
    return integrate.simps(dNda, x = a)

def integrate_SED(frequency, sed):
    nu_min = (( 6.0 * u.eV) / (const.h.cgs)).to(u.Hz).value
    nu_max = ((13.6 * u.eV) / (const.h.cgs)).to(u.Hz).value
    selection = (frequency <= nu_max) * (frequency >= nu_min)

    frequency = frequency[::-1]
    sed       = sed[::-1]

    return integrate.simps(sed[selection]/(const.c.to(u.cm/u.s).value), x = frequency[selection])

def integrate_absorption_general(E_sigma, sigma, frequency, sed, method='simps'):
    # takes in energy in Ryd
    # sigma in cm^2 / n_H
    # frequency in Hz
    # SED in _______
    nu_min = (( 6.0 * u.eV) / (const.h.cgs)).to(u.Hz).value
    nu_max = ((13.6 * u.eV) / (const.h.cgs)).to(u.Hz).value

    #
    # order SED in increasing nu
    #
    frequency = frequency[::-1]
    sed       = sed[::-1]


    selection = (frequency <= nu_max) * (frequency >= nu_min)
    #
    # set up an interpolation function for the SED data
    #
    SED_function = interpolate.interp1d(frequency[selection], 
                                            sed[selection] / (const.c.to(u.cm/u.s).value), 'linear')
   
    #
    # compute frequency from E_sigma
    #
    
    frequency_sigma = ( E_sigma * Ryd ) / (const.h.cgs.value)
    


    #
    # Sigma function
    #
    selection = (frequency_sigma <= nu_max) * (frequency_sigma >= nu_min)

    if (method == 'quad'):
        sigma_function = interpolate.interp1d(frequency_sigma[selection], sigma[selection], 'cubic')
    

        integrand = lambda x : sigma_function(x)*SED_function(x)
        return   integrate.quad(integrand, nu_min, nu_max)

    else:
        #
        # make new table of data
        #

        SED_data = SED_function(frequency_sigma[selection])
        integrand = SED_data * sigma[selection]
        return np.trapz(integrand, x = frequency_sigma[selection])



#        return integrate.romb(SED_data * sigma[selection], x = frequency_sigma[selection])
  

 #   selection = (frequency_sigma >= nu_min) * (frequency_sigma <= nu_max)
    

    #return 1.0

# do this:

#
# Load size distribution and integrate (returns N / n_H)
#
#szd_data     = np.genfromtxt(szd_file, skip_header=1)

#
# Load SED
#
#star_data     = np.genfromtxt(ostar_file, usecols=(0,10))

#
# Load cross section data
#
#crs_data      = np.genfromtxt(abs_file, skip_header=1)

#
# number of grains per number density of H
#
#N_grains    = integrate_szd(szd_data[:,0], szd_data[:,1])

# 
# integrate the SED and abs cross section
#
#abs_integral = integrate_absorption_general(crs_data[:,0], crs_data[:,1], star_data[:,0], star_data[:,1])

#print abs_integral
#
# heating raate
#
#norm_heating_rate = abs_integral * N_grains

#print "%8.8E"%(norm_heating_rate)

#n_H = 10.0 # cc
#c   = const.c.to(u.cm/u.s).value # speed of light

#E_heat = n_H*n_H*c * norm_heating_rate

#print "%8.8E"%(E_heat)

#print "%8.8E"%( integrate_SED(star_data[:,0], star_data[:,1]) )
