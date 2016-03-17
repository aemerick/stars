import numpy as np

class constants:

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
        self.tau_sun = 1.0E9 * 3.1536E13 # solar lifetime in cgs
        self.E_HI    = 13.6 # eV
        self.E_HeI   = 24.587

const = constants()
###########################################################


class RadData:

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

    def read_star_data(self):

        self.q0 = np.zeros(self._array_size)
        self.q1 = np.zeros(self._array_size)
        self.q0 = self.q0.reshape((self.NumberOfTemperatureBins, self.NumberOfSGBins,self.NumberOfMetallicityBins))
        self.q1 = self.q1.reshape((self.NumberOfTemperatureBins, self.NumberOfSGBins,self.NumberOfMetallicityBins))


        i = 0; j = 0; k = 0
        data =  np.genfromtxt('./data/q0_rates.in', usecols=(2,3,4,5,6,7,8,9,10,11))
        for line in data:

            for k in np.arange(np.size(line)):
                self.q0[i][j][np.size(line) - k - 1] = line[k]
            
            j = j + 1
            
            if (j >= self.NumberOfSGBins):
                j = 0
                i = i + 1
            
        i = 0; j = 0; k = 0
        data =  np.genfromtxt('./data/q1_rates.in', usecols=(2,3,4,5,6,7,8,9,10,11))
        for line in data:

            for k in np.arange(np.size(line)):
                self.q1[i][j][np.size(line) - k - 1] = line[k]

            j = j + 1

            if (j >= self.NumberOfSGBins):
                j = 0
                i = i + 1


    def interpolate(self, T, g, Z):        
        

        if ( T < np.min(self.T) or T > np.max(self.T)):
        
            return 0, 0, 0
        if ( g < np.min(self.g) or T > np.max(self.g)):
            return 0, 0 ,0
        if ( Z < np.min(self.Z) or T > np.max(self.Z)):
            return 0, 0, 0

        # find the nearest bin
        i     = (np.abs(self.T - T)).argmin()
        if (T < self.T[i]):
            i = i - 1
        
        j = (np.abs(self.g - g)).argmin()
        if (g < self.g[j]):
            j = j - 1
     
        k = (np.abs(self.Z - Z)).argmin()
        if (Z < self.Z[k]):
            k = k - 1

        # compute the differences
        t = (T - self.T[i])/(self.T[i+1] - self.T[i])
        u = (g - self.g[i])/(self.g[j+1] - self.g[j])
        v = (Z - self.Z[k])/(self.Z[k+1] - self.Z[k])
 
        if( self.q0[i  ][j  ][k  ] == 1.0 or
            self.q0[i  ][j+1][k  ] == 1.0 or
            self.q0[i+1][j+1][k  ] == 1.0 or
            self.q0[i+1][j  ][k  ] == 1.0 or
            self.q0[i  ][j  ][k+1] == 1.0 or
            self.q0[i  ][j+1][k+1] == 1.0 or
            self.q0[i+1][j+1][k+1] == 1.0 or
            self.q0[i+1][j  ][k+1] == 1.0   ):
       
            return 0,0,0

        q0 = (1.0 - t)*(1.0 - u)*(1.0 - v) * self.q0[i  ][j  ][k  ] +\
             (      t)*(1.0 - u)*(1.0 - v) * self.q0[i  ][j+1][k  ] +\
             (      t)*(      u)*(1.0 - v) * self.q0[i+1][j+1][k  ] +\
             (1.0 - t)*(      u)*(1.0 - v) * self.q0[i+1][j  ][k  ] +\
             (1.0 - t)*(1.0 - u)*(      v) * self.q0[i  ][j  ][k+1] +\
             (      t)*(1.0 - u)*(      v) * self.q0[i  ][j+1][k+1] +\
             (      t)*(      u)*(      v) * self.q0[i+1][j+1][k+1] +\
             (1.0 - t)*(      u)*(      v) * self.q0[i+1][j  ][k+1] 

        q1 = (1.0 - t)*(1.0 - u)*(1.0 - v) * self.q1[i  ][j  ][k  ] +\
             (      t)*(1.0 - u)*(1.0 - v) * self.q1[i  ][j+1][k  ] +\
             (      t)*(      u)*(1.0 - v) * self.q1[i+1][j+1][k  ] +\
             (1.0 - t)*(      u)*(1.0 - v) * self.q1[i+1][j  ][k  ] +\
             (1.0 - t)*(1.0 - u)*(      v) * self.q1[i  ][j  ][k+1] +\
             (      t)*(1.0 - u)*(      v) * self.q1[i  ][j+1][k+1] +\
             (      t)*(      u)*(      v) * self.q1[i+1][j+1][k+1] +\
             (1.0 - t)*(      u)*(      v) * self.q1[i+1][j  ][k+1]


        return 1, q0, q1

RadiationData = RadData()
###########################################################
class individual_star:

    def __init__(self, M, Z = 1.0, sigma_factor = 0.0):
        self.M = M
        self._M = self.M * const.Msun
        self.Z = Z
        self.sigma_factor = sigma_factor

        self.set_luminosity()
        self._L = self.L*const.Lsun
        self.set_radius()
        self._R = self.R * const.Rsun
        self.set_temperature()
        self._T = self.T * const.Tsun
        self.compute_lifetime()
        self._lifetime = self.lifetime * const.tau_sun
        self.compute_self_gravity()
        

        self.compute_ionizing_rate()

    def set_luminosity(self):
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

    def set_radius(self):
        if (self.M <= 0.5):
            self.R = self.M**(0.9)
        elif (self.M < 2.0):
            self.R = self.M **(3.0/7.0)
        else:
            self.R = self.M**(15.0/19.0)

        self.R = self.R

    def set_temperature(self):
        self.T = (self.L/self.R / self.R)**0.25

    def compute_self_gravity(self):
        self.g = const.G * (self._M)/ ((self._R)**2)

    def compute_lifetime(self):
        self.lifetime =  self.M / (self.L)



    def compute_ionizing_rate(self):

        error, self.q0, self.q1 = RadiationData.interpolate(self._T,self.g,self.Z)
        if (error == 0):
            self.ionizing_rate_blackbody()
            self.flag = 1
        else:
            self.flag = 0


    def ionizing_rate_blackbody(self):

        x = (const.E_HI / const.eV_erg) / (const.k_boltz * self._T)
        self.q0 = 4.0 * np.pi * self._R**2 * photon_radiance(x)
        x = (const.E_HeI / const.eV_erg) / (const.k_boltz * self._T)
        self.q1 = 4.0 * np.pi * self._R**2 * photon_radiance(x)

        A=  2.0 * const.k_boltz**3 * self._T**3 / (const.h**3 *const.c**2)
        self.q0 = self.q0 * A
        self.q1 = self.q1 * A


        self.E0 = average_energy(const.E_HI/const.eV_erg, self._T)
        self.E1 = average_energy(const.E_HeI/const.eV_erg, self._T)
        self.E0 = self.E0 * const.eV_erg
        self.E1 = self.E1 * const.eV_erg


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

