"""
    Author: A. Emerick
    Date  : 04/05/16
    
    Derived fields for yt meant for postprocessing analysis on individual star
    formation in Enzo with the ability to reconstruct stellar properties
    computed that are not written to file. This requires access to additional
    analysis modules and data.
"""

from imf import imf
from radiation import individual_star_properties as isp
import yt
import numpy as np

def _stellar_luminosity(field, data):
    M = data[('io','particle_mass')].convert_to_units('Msun')
    M = M.value

    Z = data[('io','metallicity_fraction')]
    Z = Z.value

    data_shape = np.shape(M)

    M = M.flatten()
    Z = Z.flatten()

    L = np.zeros(np.size(M))
    for i in range(np.size(M)):
        L[i] = isp.ZAMSData.interpolate(M[i], Z[i], return_only = 'luminosity')

    L = L.reshape(data_shape)

    L = L * isp.const.Lsun
    L = L * yt.units.erg / yt.units.s

    return L
yt.add_field('star_luminosity', function=_stellar_luminosity, units='erg*s**(-1)')

def _stellar_radius(field, data):
    # need to know star mass
    M = data['particle_mass'].convert_to_units('Msun')
    M = M.value

    Z = data[('io','metallicity_fraction')]
    Z = Z.value

    data_shape = np.shape(M)
 
    # reshape here
    M = M.flatten()
    Z = Z.flatten()
    #
    R = np.zeros(np.size(M))
    for i in range(np.size(M)):
        R[i] = isp.ZAMSData.interpolate(M[i], Z[i], return_only = 'radius')

    R = R.reshape(data_shape)

    R = R * yt.units.cm

    return R
yt.add_field('star_radius', function=_stellar_radius, units='cm')

def _stellar_temperature(field, data):
    # need to know star mass
    M = data[('io','particle_mass')].convert_to_units('Msun')
    M = M.value

    Z = data[('io','metallicity_fraction')]
    Z = Z.value

    data_shape = np.shape(M)

    M = M.flatten()
    Z = Z.flatten()

    T = np.zeros(np.size(M))
    for i in range(np.size(M)):
        T[i] = isp.ZAMSData.interpolate(M[i], Z[i], return_only = 'Teff')

    T = T.reshape(data_shape)

    T = T * yt.units.K

    return T
yt.add_field('star_Teff', function=_stellar_temperature, units='K')
