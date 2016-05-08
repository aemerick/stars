"""
  Using OSTAR2002 SED's and CLOUDY dust models, compute the 
  normalized heating rate for a given star.

  Heating rate computed is:

      E_heat / (n_H**2)

  where E_heat has units of (erg s^-1), n_H is hydrogen
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

micron    = (1.0 * u.micrometer).to(u.cm)
Ryd       = (const.h * const.c * const.Ryd).to(u.erg)

#szd_file   = './data/graphite_szd_cloudy.dat'
#ostar_file = './data/ostar2002_p03_extracted.dat'
#abs_file   = './data/graphite_ism_01.opc'

def integrate_szd(a, dNda, units='CLOUDY'):

    # if units == CLOUDY assume cloudy defaults
    # else assume the user knows what they are doing

    # convert to dNda from a^4 dNda and micron to cm
    if units == 'CLOUDY':
        a = a * u.micrometer
        dNda = dNda *( u.cm)**3 / (a*a*a*a)

    print np.min(dNda), np.max(dNda)
    # now integrate from a min to a max
    return np.trapz( dNda.to(u.cm**(-1)).value, x = a.to(u.cm).value)

def integrate_SED(frequency, sed, nu_min = None, nu_max = None):
    if nu_min == None:
        nu_min = (( 6.0 * u.eV) / (const.h.cgs)).to(u.Hz).value
    if nu_max == None:
        nu_max = ((13.6 * u.eV) / (const.h.cgs)).to(u.Hz).value
    
    selection = (frequency <= nu_max) * (frequency >= nu_min)

    if frequency[-1] < frequency[0]: # sort in ascending order
        frequency = frequency[::-1]
        sed       = sed[::-1]

    

    return 4.0 * np.pi * np.trapz(sed[selection], x = frequency[selection])

def integrate_absorption_general(E_sigma, sigma, frequency, sed, W, method='simps'):
    # W assumed in Ryd


    # takes in energy in Ryd
    # sigma in cm^2 / n_H
    # frequency in Hz
    # SED in _______
#    nu_min = (( 6.0 * u.eV) / (const.h.cgs)).to(u.Hz).value
    #
    # Integrate from work function up to H ionization energy
    #
    nu_min = (( W * Ryd ) / (const.h.cgs))
    nu_max = ((13.6 * u.eV) / (const.h.cgs))


    #
    # order SED in increasing nu
    #
    frequency = frequency[::-1] * u.Hz
    sed       = sed[::-1]       * u.erg / u.s / u.cm**2 / u.Hz

    if frequency[-1] < frequency[0]:
        print "need to re order SED frequency"
    if E_sigma[-1] < E_sigma[0]:
        print "need to re order energy"

    #
    # a little slop on either end for interpoaltion purposes
    #
    selection = (frequency <= nu_max*1.02) * (frequency >= nu_min*0.98)
    #
    # set up an interpolation function for the SED data
    #
    SED_function = interpolate.interp1d((frequency.cgs.value)[selection], (sed.cgs.value)[selection], 'linear')
   
    #
    # compute frequency from E_sigma
    #
    
    W               = ( W * Ryd ) # now in erg
#    frequency_sigma = ( E_sigma * Ryd ) / (const.h.cgs.value)
    frequency_sigma = (E_sigma * (Ryd / const.h.cgs) )
   
    #
    # Sigma function
    #
    selection = (frequency_sigma <= nu_max) * (frequency_sigma >= nu_min)

    sigma = sigma * u.cm * u.cm

    if (method == 'quad'):
        selection = (frequency_sigma <= nu_max*1.02) * (frequency_sigma >= nu_min*0.98)

        sigma_function = interpolate.interp1d((frequency_sigma.cgs.value)[selection], (sigma.cgs.value)[selection], 'linear')
    

        integrand = lambda x : sigma_function(x)*SED_function(x) * 0.5
        return   integrate.quad(integrand, nu_min.cgs.value, nu_max.cgs.value)[0]

    else:
        #
        # make new table of data
        #

        SED_data  = SED_function(frequency_sigma[selection].cgs.value)  # Flux

        # ratio of ejected electron energy to photon energy
        energy_ratio = (const.h.cgs * frequency_sigma[selection].cgs - W.cgs) / (const.h.cgs*frequency_sigma[selection].cgs)
        

                    # abs / n_H       *  F_nu     * (h*nu - W) / (h*nu)
        integrand = sigma[selection]  * SED_data  * energy_ratio


      
        norm_heating = np.trapz(integrand.cgs.value, x = frequency_sigma[selection].cgs.value)

        return norm_heating


        #/    nst.h.cgs.value * frequency_sigma[selection])
        #integrand = SED_data *  (const.h.cgs.value * frequency_sigma[selection] - W)

         # integrand  = integrand * sigma[selection]
        #return np.trapz(integrand, x = frequency_sigma[selection])



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
