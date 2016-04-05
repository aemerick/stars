"""
    Author: A. Emerick
    Date  : 4/4/16

    Class meant for evolution of a single stellar population 
    formed at the same time with no ongoing star formation. Desire
    to fix this in the future. This is meant to do simple tests 
    of ionizing radiation and feedback properties and associated 
    energies over time for a given stellar mass galaxy without 
    having to run a full Enzo simulation
"""

from __future__ import division

import numpy as np
from imf import imf
from radiation import individual_star_properties as isp


class population:
    """
    Population of stars
    """

    def __init__(self, IMF = None, IMF_kwargs = {},
                       IMF_function = None,
                       M_min = 1.0, M_max = 100.0   ):

        self.t = 0.0 # current time
        self.M_min = M_min
        self.M_max = M_max

        if IMF_function != None:
            # user defined function
            if IMF != None:
                self.IMF = IMF
            else:
                self.IMF = 'user'
            self.IMF_function = IMF_function

            return

        
        self.IMF = IMF

        if self.IMF == None:
            self.IMF_function = None
        elif self.IMF == 'salpeter':
            self.IMF_function =\
                    lambda x: imf.salpeter(x, **IMF_kwargs)
        elif self.IMF == 'diet_salpeter':
            self.IMF_function =\
                    lambda x: imf.diet_salpeter(x, **IMF_kwargs)



        return

    def generate_population(self, Nstars, metallicity = None):
        """
        Sample from IMF
        """

        self.Nstars = Nstars
        self.M0     = imf.sample_IMF(Nstars, self.IMF_function,
                                        M_min =self.M_min, M_max =self.M_max)

        if np.size(metallicity) < self.Nstars:
            self.Z0     = metallicity * np.ones(self.Nstars)
        else:
            self.Z0     = metallicity

        self.populate_stellar_properties()
        self._define_status_lookup()
        return 

    def populate_stellar_properties(self):
        """
        compute stellar properties
        """

        self.L = np.zeros(self.Nstars)
        self.T = np.zeros(self.Nstars)
        self.R = np.zeros(self.Nstars)

        for i in np.arange(self.Nstars):
            err, L, T, R = isp.ZAMSData.interpolate(self.M0[i], self.Z0[i])

            self.L[i] = L 
            self.T[i] = T / isp.const.Tsun
            self.R[i] = R / isp.const.Rsun

        self.g        = isp.const.G * (self.M0 * isp.const.Msun)/(self.R*isp.const.Rsun)**2
        self.lifetime = (self.M0 / self.L) * isp.const.tau_sun

        self.status0   = np.zeros(self.Nstars)

        return None

    def _define_status_lookup(self):

        self.status_codes = {0 : 'MS',
                             1 : 'MS_rad',
                             2 : 'SNII',
                             3 : 'WD',
                             4 : 'WD_SN1a',
                             5 : 'SN1a'}
        self._inv_status_codes = {v: k for k, v in self.status_codes.items()}

        return

    def _promote_dimensions(self, nsteps):
        """
        Change 1D arrays to 2D
        """
        self.M = np.zeros((nsteps, self.Nstars))
        self.Z = np.zeros((nsteps, self.Nstars))
        self.status = np.zeros((nsteps, self.Nstars))

        return None

    def evolve(self, t_start, t_end, dt):
        """

        """
        Myr = 3.1536E13
        yr  = Myr / 1.0E6
        t_start = t_start * Myr
        t_end   = t_end   * Myr
        dt      = dt * Myr


        self.t = np.arange(t_start, t_end + dt, dt)

        nsteps = np.size(self.t)
        # make all properties 2D arrays
        self._promote_dimensions(nsteps)

        # now assign initial values
        self.M[0] = self.M0
        self.Z[0] = self.Z0
        self.status[0] = self.status0

        for i in range(1, nsteps):

            self.status[i] = self.status[i-1]
            self.M[i]      = self.M[i-1]
            self.Z[i]      = self.Z[i-1]

            # does the star die?
            new_dead = (self.t[i] > self.lifetime) * (self.t[i-1] < self.lifetime)
                
            SNII   = new_dead * (self.M[i-1] > 8.0)                        * (self.status[i-1] <= 1)
            WD     = new_dead * (self.M[i-1] <= 3.0) * (self.M[i-1] > 1.0) * (self.status[i-1] <= 1)
            WDSN1a = new_dead * (self.M[i-1] <=8.0) * (self.M[i-1] > 3.0) * (self.status[i-1] <= 1)

            # Salaris et. al. 2009	
            self.M[i][WD * (self.M[i]<4.0)] = 0.134 * self.M[i][WD * (self.M[i]<4.0)] + 0.331
            self.M[i][WD * (self.M[i]>4.0)] = 0.047 * self.M[i][WD * (self.M[i]>4.0)] 

            self.M[i][SNII] = 0.0 # gone

            self.status[i][WD    ] = self._inv_status_codes['WD']
            self.status[i][WDSN1a] = self._inv_status_codes['WD_SN1a']
            self.status[i][SNII  ] = self._inv_status_codes['SNII']


            # do stellar wind junk here

            # check for SN1a
            if any( s == self._inv_status_codes['WD_SN1a'] for s in self.status[i-1]):

                WD = (self.status[i-1] == self._inv_status_codes['WD_SN1a'])

                # progenitor masses
                prog     = self.M[0][WD]

                IMF_norm = prog**(-2.35)
                
                msf = np.zeros(np.size(prog))

                for j in np.arange(np.size(prog)):
           #         if prog[i] < 0.0: # set to zero, makes probability zero
            #            msf[i] = 0.0
             #       else:
                    msf[j] = IMF_norm[j] * ( imf.integrate_imf(self.IMF_function, self.M_min, self.M_max) )

                time_elapsed = self.t[i] - self.t[0]

                probability = imf.DTD(time_elapsed) * msf * (dt / yr)

	

                # choose random numbers
                rnum = np.random.rand(np.size(probability))
                
                new_status = np.zeros(np.size(rnum))
                new_status[rnum <= probability] = self._inv_status_codes['SN1a']
                new_status[rnum > probability]  = self._inv_status_codes['WD_SN1a']

                self.status[i][WD] = new_status
                self.M[i][self.status[i] == self._inv_status_codes['SN1a']] = 0.0

            # endif white dwarfs SN1a check

        self.t = self.t /Myr
        # evolution complete
        return None

def accumulate(starpop, particle_id, scalar=False):

    cum_dist = np.zeros(np.size( starpop.t))
    dist = np.zeros(np.size(cum_dist))

    for i in np.arange(0, np.size(starpop.t)):
        dist[i] = np.size( starpop.status[i][ starpop.status[i] == starpop._inv_status_codes[particle_id] ])
        

    if particle_id == 'SN' or particle_id == 'SN1a':
        cum_dist    = dist * 1.0
        for i in np.arange(1, np.size(starpop.t)):
            dist[i]     = cum_dist[i] - np.sum(cum_dist[:i])


    if scalar:
        if np.sum(cum_dist) == 0.0:
            return np.sum(dist)
        else:
            return cum_dist[-1]
    else:
        return dist, cum_dist
