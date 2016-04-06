"""

 Author: Andrew Emerick
 affil : Columbia University - Dept. of Astronomy
         American Museum of Natural History - Dept. of Astrophysics
 date  : March 2016
 

 Description:

     This is a rough translation of C++ code I've put into Enzo 
     to compute stellar properties and their ionizing radiation
     rates. This code is purposefully non-pythonic to be as 
     close of a copy to the C++ code as possible. Purpose of this 
     is for testing and easily plotting of stellar properties without
     having to save all properties at each output in Enzo (this would be
     memory intensive and a pain to code given how the star objects work).

"""

import numpy as np

class constants:
    """ 
    Helpful contants. In cgs or cgs conversions 
    except ionization energies (eV)
    """

    def __init__(self):
        self.eV_erg  = 6.24150934326E11
        self.k_boltz = 1.380658E-16
        self.c       = 2.99792458E10
        self.h       = 6.6260755E-27
        self.G       = 6.6743E-8
        self.Msun = 1.998E33 # solar masses in cgs
        self.Rsun = 69.63E9  # solar radius in cgs
        self.Lsun = 3.9E33   # solar luminosity in cgs
        self.Tsun = 5777.0   # solar temperature in cgs
        self.tau_sun = 10.0E9 * 3.1536E7 # solar lifetime in cgs
        self.E_HI    = 13.6 # eV
        self.E_HeI   = 24.587
        self.Zsolar_ostar  = 0.01700 # Grevesse & Sauval 1998 - used in OSTAR2002
        self.Zsolar_parsec = 0.01524     # Caffau et. al. 2009 / 2011 - used in PARSEC SE code

        return None

# making a global object to be used below
const = constants()
###########################################################
RADDATADIR  = '/home/emerick/Research/stars/radiation/data/'
ZAMSDATADIR = '/home/emerick/Research/stars/radiation/sissa_data/'

class RadData:
    """
    RadData object to mimick both the IndividualStarRadData struct in Enzo,
    but also includes the routine for interpolating along the tabulated
    OSTAR2002 data.

    Tabulated data dimensions and T, g, Z bin values are hardcoded
    """

    def __init__(self):
        self.NumberOfTemperatureBins = 12
        self.NumberOfSGBins          = 8
        self.NumberOfMetallicityBins = 10
        self._array_size = self.NumberOfTemperatureBins * self.NumberOfSGBins * self.NumberOfMetallicityBins

        self.T = np.arange(27500.0, 57500.0, 2500.0)
        self._logg = np.arange(3.0, 5.0, 0.25)
        self.g = 10.0**(self._logg)
        self.Z = np.array([0.0,0.001,0.01,1/50.0,1/30.0,0.1,0.2,0.5,1.0,2.0])

        self.read_star_data()

        return None

    def read_star_data(self):
        """
        Read the data from file in 3D array (as is done in C++ code)
        """

        self.q0 = np.zeros(self._array_size)
        self.q1 = np.zeros(self._array_size)
        self.q0 = self.q0.reshape((self.NumberOfTemperatureBins, self.NumberOfSGBins,self.NumberOfMetallicityBins))
        self.q1 = self.q1.reshape((self.NumberOfTemperatureBins, self.NumberOfSGBins,self.NumberOfMetallicityBins))


        i = 0; j = 0; k = 0
        data =  np.genfromtxt(RADDATADIR + 'q0_rates.in', usecols=(2,3,4,5,6,7,8,9,10,11))
        for line in data:

            for k in np.arange(np.size(line)):
                self.q0[i][j][np.size(line) - k - 1] = line[k]
            
            j = j + 1
            
            if (j >= self.NumberOfSGBins):
                j = 0
                i = i + 1
            
        i = 0; j = 0; k = 0
        data =  np.genfromtxt(RADDATADIR + 'q1_rates.in', usecols=(2,3,4,5,6,7,8,9,10,11))
        for line in data:

            for k in np.arange(np.size(line)):
                self.q1[i][j][np.size(line) - k - 1] = line[k]

            j = j + 1

            if (j >= self.NumberOfSGBins):
                j = 0
                i = i + 1
      
        # un-log q values
        self.q0 = 10.0**(self.q0)
        self.q1 = 10.0**(self.q1)


        return None

    def interpolate(self, T, g, Z):        
        """
        Tri-linear interpolation of data given desired point in T, g, Z space

        Returning 0 signifies a failure and need to do blackbody curve integration.
        Returning 1 signifies a successful interpolation

        """

        # make sure T, g, and Z are in bounds of tabulated data
        if ( T < np.min(self.T) or T > np.max(self.T)):
            return 0, 0, 0
        if ( g < np.min(self.g) or g > np.max(self.g)):
            return 0, 0 ,0
        if ( Z < np.min(self.Z) or Z > np.max(self.Z)):
            return 0, 0, 0

        # find the nearest bin in each dimension (cheating from C++ version)
        i     = (np.abs(self.T - T)).argmin()
        if (T < self.T[i]):
            i = i - 1
        
        j = (np.abs(self.g - g)).argmin()
        if (g < self.g[j]):
            j = j - 1
     
        k = (np.abs(self.Z - Z)).argmin()
        if (Z < self.Z[k]):
            k = k - 1

        # compute the ratios from nearest neighber
        t = (T - self.T[i])/(self.T[i+1] - self.T[i])
        u = (g - self.g[j])/(self.g[j+1] - self.g[j])
        v = (Z - self.Z[k])/(self.Z[k+1] - self.Z[k])
 
        # a q value of 1.0 means a logq value of 0.00 in the table,
        # which means that this combination of T, g is untabulated
        if( self.q0[i  ][j  ][k  ] == 1.0 or
            self.q0[i  ][j+1][k  ] == 1.0 or
            self.q0[i+1][j+1][k  ] == 1.0 or
            self.q0[i+1][j  ][k  ] == 1.0 or
            self.q0[i  ][j  ][k+1] == 1.0 or
            self.q0[i  ][j+1][k+1] == 1.0 or
            self.q0[i+1][j+1][k+1] == 1.0 or
            self.q0[i+1][j  ][k+1] == 1.0   ):
       
            return 0,0,0

        # compute and return values
#        q0 = (1.0 - t)*(1.0 - u)*(1.0 - v) * self.q0[i  ][j  ][k  ] +\
#             (      t)*(1.0 - u)*(1.0 - v) * self.q0[i  ][j+1][k  ] +\
#             (      t)*(      u)*(1.0 - v) * self.q0[i+1][j+1][k  ] +\
#             (1.0 - t)*(      u)*(1.0 - v) * self.q0[i+1][j  ][k  ] +\
#             (1.0 - t)*(1.0 - u)*(      v) * self.q0[i  ][j  ][k+1] +\
#             (      t)*(1.0 - u)*(      v) * self.q0[i  ][j+1][k+1] +\
#             (      t)*(      u)*(      v) * self.q0[i+1][j+1][k+1] +\
#             (1.0 - t)*(      u)*(      v) * self.q0[i+1][j  ][k+1] 
 
        q0 = (1.0 - t)*(1.0 - u)*(1.0 - v) * self.q0[i  ][j  ][k  ] +\
             (1.0 - t)*(      u)*(1.0 - v) * self.q0[i  ][j+1][k  ] +\
             (      t)*(      u)*(1.0 - v) * self.q0[i+1][j+1][k  ] +\
             (      t)*(1.0 - u)*(1.0 - v) * self.q0[i+1][j  ][k  ] +\
             (1.0 - t)*(1.0 - u)*(      v) * self.q0[i  ][j  ][k+1] +\
             (1.0 - t)*(      u)*(      v) * self.q0[i  ][j+1][k+1] +\
             (      t)*(      u)*(      v) * self.q0[i+1][j+1][k+1] +\
             (      t)*(1.0 - u)*(      v) * self.q0[i+1][j  ][k+1] 

        q1 = (1.0 - t)*(1.0 - u)*(1.0 - v) * self.q1[i  ][j  ][k  ] +\
             (1.0 - t)*(      u)*(1.0 - v) * self.q1[i  ][j+1][k  ] +\
             (      t)*(      u)*(1.0 - v) * self.q1[i+1][j+1][k  ] +\
             (      t)*(1.0 - u)*(1.0 - v) * self.q1[i+1][j  ][k  ] +\
             (1.0 - t)*(1.0 - u)*(      v) * self.q1[i  ][j  ][k+1] +\
             (1.0 - t)*(      u)*(      v) * self.q1[i  ][j+1][k+1] +\
             (      t)*(      u)*(      v) * self.q1[i+1][j+1][k+1] +\
             (      t)*(1.0 - u)*(      v) * self.q1[i+1][j  ][k+1] 


#        q1 = (1.0 - t)*(1.0 - u)*(1.0 - v) * self.q1[i  ][j  ][k  ] +\
#             (      t)*(1.0 - u)*(1.0 - v) * self.q1[i  ][j+1][k  ] +\
#             (      t)*(      u)*(1.0 - v) * self.q1[i+1][j+1][k  ] +\
#             (1.0 - t)*(      u)*(1.0 - v) * self.q1[i+1][j  ][k  ] +\
#             (1.0 - t)*(1.0 - u)*(      v) * self.q1[i  ][j  ][k+1] +\
#             (      t)*(1.0 - u)*(      v) * self.q1[i  ][j+1][k+1] +\
#             (      t)*(      u)*(      v) * self.q1[i+1][j+1][k+1] +\
#             (1.0 - t)*(      u)*(      v) * self.q1[i+1][j  ][k+1]


        return 1, q0, q1

# make the global variable
RadiationData = RadData()
###########################################################


class ZAMS_data:


    def __init__(self):
        self.M = np.array( [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0,
                            10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 24.0, 28.0, 30.0, 
                            35.0, 40.0, 50.0, 55.0, 60.0, 70.0, 90.0, 100.0, 120.0])
        self.Z = np.array( [0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.004, 0.006, 0.008,
                            0.01, 0.014, 0.017] )
 
        self.NumberOfMassBins = np.size(self.M)
        self.NumberOfMetallicityBins = np.size(self.Z)    
        self._array_size = self.NumberOfMassBins * self.NumberOfMetallicityBins


        self.read_data()
        return None

    def read_data(self):

        
        self.L    = np.zeros(self._array_size)
        self.Teff = np.zeros(self._array_size)
        self.R    = np.zeros(self._array_size)
        self.L    = self.L.reshape((self.NumberOfMassBins, self.NumberOfMetallicityBins))        
        self.Teff = self.Teff.reshape((self.NumberOfMassBins,self.NumberOfMetallicityBins))
        self.R    = self.R.reshape((self.NumberOfMassBins, self.NumberOfMetallicityBins))

        i = 0; j = 0
        data = np.genfromtxt(ZAMSDATADIR + 'ZAMS_data.in')
        for line in data:
            self.L[i][j]    = 10**(line[2])
            self.Teff[i][j] = 10**(line[3])
            self.R[i][j]    = 10**(line[4])
            
            j = j + 1
            if ( j >= self.NumberOfMetallicityBins):
                j = 0
                i = i +1
         


        return None

    def interpolate(self, M, Z, return_only = None):
        """
        Bi-linear interpolation of data given desired point in M, Z space

        Returning 0 signifies a failure 
        Returning 1 signifies a successful interpolation
        """

        # make sure T, g, and Z are in bounds of tabulated data
        if ( M < np.min(self.M) or M > np.max(self.M)):
            print "Mass out of bounds"
            if return_only == None:
                return 0, 0, 0, 0
            else:
                return 0
        if ( Z < np.min(self.Z)):
            print "metallicity too low"
            if return_only == None:
                return 0, 0, 0, 0
            else:
                return 0
        if ( Z > np.max(self.Z)):
            print "warning: making ZAMS interpolation metallicity just under max"
            Z = 0.99 * np.max(self.Z)

        # find the nearest bin in each dimension (cheating from C++ version)
        i     = (np.abs(self.M - M)).argmin()
        if (M <= self.M[i]):
            i = i - 1

        j = (np.abs(self.Z - Z)).argmin()
        if (Z <= self.Z[j]):
            j = j - 1

        # compute the ratios from nearest neighber
        t = (M - self.M[i])/(self.M[i+1] - self.M[i])
        u = (Z - self.Z[j])/(self.Z[j+1] - self.Z[j])

        L    = (1.0 - t)*(1.0 - u) * self.L[i  ][j  ] +\
               (1.0 - t)*(      u) * self.L[i  ][j+1] +\
               (      t)*(      u) * self.L[i+1][j+1] +\
               (      t)*(1.0 - u) * self.L[i+1][j  ]

        Teff = (1.0 - t)*(1.0 - u) * self.Teff[i  ][j  ] +\
               (1.0 - t)*(      u) * self.Teff[i  ][j+1] +\
               (      t)*(      u) * self.Teff[i+1][j+1] +\
               (      t)*(1.0 - u) * self.Teff[i+1][j  ]

        R    = (1.0 - t)*(1.0 - u) * self.R[i  ][j  ] +\
               (1.0 - t)*(      u) * self.R[i  ][j+1] +\
               (      t)*(      u) * self.R[i+1][j+1] +\
               (      t)*(1.0 - u) * self.R[i+1][j  ]

        if return_only == None:

            return 1, L, Teff, R
        elif return_only == 'radius':
            return R
        elif return_only == 'luminosity':
            return L
        elif return_only == 'Teff':
            return Teff


ZAMSData = ZAMS_data()

###############################################################
class individual_star:
    """
    Mimicking code to compute stellar properties.
    """

    def __init__(self, M, Z = 1.0, use_ZAMS = False, sigma_factor = 0.0, blackbody_only = False):
        """
        NOTE: M, Z, L, R, T are in solar units. CGS values 
        are stored as _M, _Z, _R, _L, _T
        """

        self.M = M
        self._M = self.M * const.Msun
        self.Z  = Z
        self._Z = self.Z * const.Zsolar_parsec

        if use_ZAMS:
           self.compute_ZAMS_parameters()

        else:

            self.sigma_factor = sigma_factor

            self.set_luminosity()
            self._L = self.L*const.Lsun
            self.set_radius()
            self._R = self.R * const.Rsun
            self.set_temperature()
            self._T = self.T * const.Tsun


        self.compute_self_gravity()        
        self.compute_lifetime()
        self._lifetime = self.lifetime * const.tau_sun
        self.blackbody_only = blackbody_only
        self.compute_ionizing_rate()

        return None

    def compute_ZAMS_parameters(self):
        """
        Use ZAMS data to find stellar parameters
        """
        
        err, L, T, R = ZAMSData.interpolate(self.M, self._Z)

        self._L  = L
        self.L = self._L / const.Lsun

        self._T = T
        self.T  = self._T / const.Tsun
 
        self._R = R
        self.R  = self._R / const.Rsun

        return None

    def set_luminosity(self):
        """
        NOTE: Major difference between here and Enzo is we do not do random
        number call to account for spread in L-M relation. Rather, make
        deterministic and specify spread at initialization in factors of
        one standard deviation.
        Ekel et. al. 2015
        """

        if ( self.M <= 1.05):
            alpha = 4.841
            A     = -0.026
            sigma = 0.121
        elif ( 1.05 <= self.M and self.M <=2.40):
            alpha = 4.328
            A = -0.002
            sigma = 0.108
        elif (2.40 <= self.M and self.M <= 7.00):
            alpha = 3.962
            A = 0.120
            sigma = 0.165
        else:
            alpha = 2.726
            A     = 1.237
            sigma = 0.158

        self.L = 10.0**(alpha * np.log10(self.M) + A + sigma*self.sigma_factor)

        return None

    def set_radius(self):
        # make functions continious

        slopes = [0.9, 3.0/7.0, 15.0/19.0, 9.0/19.0]
        if (self.M <= 0.5):
            self.R = self.M**(slopes[0])
        elif (self.M < 2.0): # normalize everything relative to this (solar)
            self.R = self.M **(slopes[1]) # solar
        elif (self.M < 20.0):
            A = 2.0**(slopes[1] - slopes[2])
            self.R = (A) * self.M**(slopes[2])
        else:
            slope = 9.0/19.0
            A = 2.0**(slopes[1] - slopes[2])
            B = A * 20.0**(slopes[2] - slopes[3])
            self.R = B * self.M**(slopes[3])

        return None

    def set_temperature(self):
        self.T = (self.L/self.R / self.R)**0.25
        return None

    def compute_self_gravity(self):
        self.g = const.G * (self._M)/ ((self._R)**2)

        return None

    def compute_lifetime(self):
        self.lifetime =  self.M / (self.L)
        return None


    def compute_ionizing_rate(self):
        """
        Tries to interpolate over radiation data. Does black body if this does not work
        """
        blackbody_flag = 1
        table_flag     = 0

        if self.blackbody_only:
            self.ionizing_rate_blackbody()
            self.flag = blackbody_flag
        else:
            error, self.q0, self.q1 = RadiationData.interpolate(self._T,self.g,
                                                                self._Z/const.Zsolar_ostar)
            if (error == 0): # use black body instead
                self.ionizing_rate_blackbody()
                self.flag = blackbody_flag
                if False:
                    if (self.M > 20.0):
                        self.q0 = self.q0 * 2.89
                        self.q1 = self.q1 * 5.2
                    if (self.M < 20.0):
                       self.q0 = self.q0 / 10.0
                       self.q1 = self.q1 / 100.0
            else:
                self.flag = table_flag

        # convert to photons / s
        surface_area = 4.0 * np.pi * self._R**2
        self.q0 = self.q0 * surface_area
        self.q1 = self.q1 * surface_area

        # now compute the average energy using a black body
        self.E0 = average_energy(const.E_HI/const.eV_erg, self._T)
        self.E1 = average_energy(const.E_HeI/const.eV_erg, self._T)
        self.E0 = self.E0 * const.eV_erg
        self.E1 = self.E1 * const.eV_erg



        return None

    def ionizing_rate_blackbody(self):

        x = (const.E_HI / const.eV_erg) / (const.k_boltz * self._T)
        self.q0 = photon_radiance(x)
        x = (const.E_HeI / const.eV_erg) / (const.k_boltz * self._T)
        self.q1 = photon_radiance(x)

        A=  2.0 * const.k_boltz**3 * self._T**3 / (const.h**3 *const.c**2)
        self.q0 = self.q0 * A
        self.q1 = self.q1 * A

        return None

def photon_radiance(x):
    max_iter = 513
    min_iter = 4
    tolerance = 1.0E-10

    difference = 1.0E10

    sum = 0.0; old_sum = 0.0
    i = 1
    while((difference > tolerance and i < max_iter) or i < min_iter):
        old_sum = sum
        sum += ( x*x/(1.0*i) + 2.0*x/(1.0*i*i) + 2.0/(1.0*i*i*i))*np.exp(-i*x)
        difference = sum - old_sum
        i = i + 1

    return sum


def average_energy(E_i, T):
    max_iter = 513 
    min_iter = 4
    tolerance = 1.0E-10
    x = E_i / (const.k_boltz * T)


    difference = 1.0E10;
    sum = old_sum = 0.0;
    i = 1;
    while((difference > tolerance and i < max_iter) or i < min_iter):
        old_sum = sum
        sum += ( x*x*x/(1.0* i) + 3.0*x*x/(1.0* i*i)
                     + 6.0*x/(1.0* i*i*i) + 6.0 / (1.0*i*i*i*i)) * np.exp(-i*x)
        difference = sum - old_sum
        i = i + 1

    u_dens_sum = sum

    sum = 0.0; old_sum = 0.0; i = 1; difference = 1.0E10
 
    while( (difference > tolerance and i < max_iter) or i < min_iter):
        old_sum = sum
        sum += ( x*x/(1.0*i) + 2.0*x/(1.0*i*i) + 2.0/(1.0*i**3))*np.exp(-i*x)
        difference = sum - old_sum
        i = i + 1

    return (const.k_boltz * T)*(u_dens_sum / sum)

def compute_blackbody_rate(T):

    x = (const.E_HI / const.eV_erg) / (const.k_boltz * T)
    q0 = photon_radiance(x)
    x = (const.E_HeI / const.eV_erg) / (const.k_boltz * T)
    q1 = photon_radiance(x)

    A=  2.0 * const.k_boltz**3 * T**3 / (const.h**3 *const.c**2)
    q0 = q0 * A
    q1 = q1 * A

    return q0, q1