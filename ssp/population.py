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
                       M_min = 1.0, M_max = 100.0,
                       SFR   = 0.0, M_gas = 0.0, Zgas = 0.01, max_stars = 1.0E6):

        self.t = 0.0 # current time
        self.M_min = M_min
        self.M_max = M_max
        self.Nstars = 0
        self.max_stars = 1.0E6
        # allow for continual star formation from a gas resivoir
        self.SFR     = SFR
        self.M0_gas  = M_gas
        self.Zgas    = Zgas
      
        
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

    def generate_population(self, X, sample_method, metallicity = None):
        """
        Sample from IMF
        X is number of stars or star mass depending on sample method
        sample method is either "star_count" or "stellar_mass"
        """

        if sample_method == "star_count":
            
            self.Nstars = X
            self.M0     = imf.sample_IMF(self.IMF_function, N = self.Nstars,
                                         M_min =self.M_min, M_max =self.M_max)

        elif sample_method == "stellar_mass":
            self.M0     = imf.sample_IMF(self.IMF_function, M = X,
                                         M_min =self.M_min, M_max =self.M_max)
            self.Nstars = np.size(self.M0)
            
        
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
        self.formation_time = np.zeros(self.Nstars)
        
        self.q0 = np.zeros(self.Nstars)
        self.q1 = np.zeros(self.Nstars)
        
        for i in np.arange(self.Nstars):
            err, L, T, R = isp.ZAMSData.interpolate(self.M0[i], self.Z0[i])

            self.L[i] = L 
            self.T[i] = T #/ isp.const.Tsun
            self.R[i] = R #/ isp.const.Rsun

        self.g        = isp.const.G * (self.M0 * isp.const.Msun)/(self.R)**2
        self.lifetime = (self.M0 / self.L) * isp.const.tau_sun

        
        # compute radiation properties:
        for i in np.arange(self.Nstars):
            error, q0, q1 = isp.RadiationData.interpolate(self.T[i], self.g[i],
                                                 self.Z0[i]/isp.const.Zsolar_ostar)
            
            if error == 0:
                # ionizing rate blackbody
                self.q0[i], self.q1[i] = isp.compute_blackbody_rate(self.T[i])
                
        self.q0 = self.q0 * 4.0 * np.pi * self.R**2
        self.q1 = self.q1 * 4.0 * np.pi * self.R**2
        
        self.status0   = np.zeros(self.Nstars)

        return None
    
    def _new_star_properties(self, t_i, imin, imax):
        
        for i in np.arange(imin, imax):
            err, L, T, R = isp.ZAMSData.interpolate(self.M[t_i][i], self.Z[t_i][i])

            self.L[i] = L 
            self.T[i] = T / isp.const.Tsun
            self.R[i] = R / isp.const.Rsun

        self.g[imin:imax]        = isp.const.G * (self.M[t_i][imin:imax] *
                                                  isp.const.Msun)/(self.R[imin:imax]*isp.const.Rsun)**2
        self.lifetime[imin:imax] = (self.M[t_i][imin:imax] / self.L[imin:imax]) * isp.const.tau_sun            

        return None
    
    def _define_status_lookup(self):

        self.status_codes = {999 : 'NOSTAR',
                             0   : 'MS',
                             1   : 'MS_rad',
                             2   : 'SNII',
                             3   : 'WD',
                             4   : 'WD_SN1a',
                             5   : 'SN1a'}
        self._inv_status_codes = {v: k for k, v in self.status_codes.items()}

        return

    def _promote_dimensions(self, nsteps, stars):
        """
        Change 1D arrays to 2D
        """
        self.M = np.zeros((nsteps, stars))
        self.Z = np.zeros((nsteps, stars))
        self.status = np.zeros((nsteps, stars))

        return None

    def _resize_for_active_star_formation(self, max_stars):
        """
        Resizes arrays to allow for star formation but cap it at a maximum
        number, otherwise we will have memory issues
       
        This handles all stellar properties
        """
        
        self.L = np.resize( self.L, max_stars)
        self.T = np.resize( self.T,  max_stars)
        self.R = np.resize( self.R, max_stars)
        self.g = np.resize( self.g, max_stars)
        self.lifetime = np.resize( self.lifetime, max_stars)
        self.formation_time = np.resize( self.formation_time, max_stars)

        # check if these exist
        if self.Nstars > 0:        
            self.M0 = np.resize(self.M0, max_stars)
            self.Z0 = np.resize(self.Z0, max_stars)
            self.status0 = np.resize(self.status0, max_stars)
            
            self.M0[self.Nstars:] = 0.0
            self.Z0[self.Nstars:] = 0.0
            
            self.status0[self.Nstars:] = 999
         
        
            self.L[self.Nstars:] = 0.0
            self.T[self.Nstars:] = 0.0
            self.R[self.Nstars:] = 0.0
            self.g[self.Nstars:] = 0.0
            self.lifetime[self.Nstars] = 0.0
            self.formation_time[self.Nstars:] = 0.0
    
    def evolve(self, t_start, t_end, dt):
        """

        if SFR is given, then do something a little bit more interesting

        """
        Myr = 3.1536E13
        yr  = Myr / 1.0E6
        
        t_start = t_start * Myr
        t_end   = t_end   * Myr
        dt      = dt * Myr
        self.t = np.arange(t_start, t_end + dt, dt)
        nsteps = np.size(self.t)
        
        # now, if SFR is on, need to do something a little bit weirder
        if self.SFR > 0.0:
            self._resize_for_active_star_formation( self.max_stars)
            self._promote_dimensions(nsteps, self.max_stars)
            
            self.M_gas    = np.zeros(nsteps)
            self.M_gas[0] = self.M0_gas
        else:
            self.formation_time = np.ones(self.Nstars) * self.t[0]
            self._promote_dimensions(nsteps, self.Nstars)
        
        # make all properties 2D arrays

        # now assign initial values if they exist
        if self.Nstars > 0:
            self.M[0]      = self.M0
            self.Z[0]      = self.Z0
            self.status[0] = self.status0
        else:
            self.status[0] = -1.0 * np.ones(np.size(self.status[0]))

        for i in range(1, nsteps):

            # copy over previous first, then modify if something happens to it below
            self.status[i] = self.status[i-1]
            self.M[i]      = self.M[i-1]
            self.Z[i]      = self.Z[i-1]
            
            # form stars if gas is available and SFR is ON
            if self.SFR > 0.0 and self.M_gas[i-1] > 0.0:
                # SFR in units of solar masses per year
                # M_gas in units of solar masses
                
                SF = self.SFR * (dt / yr) # desired mass to form this timestep
                
                mass_sample = imf.sample_IMF(self.IMF_function, M = SF,
                                             M_min =self.M_min, M_max =self.M_max)
                
                SF = np.sum(mass_sample) # actual mass of stars
                
                if SF > self.M_gas[i-1]:
                    SF = self.M_gas[i-1]
                    self.M_gas[i] = 0.0
                else:
                    self.M_gas[i] = self.M_gas[i-1] - SF
                    
                new_stars = np.size(mass_sample)
                
                if self.Nstars + new_stars > self.max_stars:
                    print "Making too many stars", self.Nstars, new_stars, self.max_stars
                    return None
                                            
                self.M[i,self.Nstars: self.Nstars + new_stars] = mass_sample
                self.Z[i,self.Nstars: self.Nstars + new_stars] = self.Zgas
                self.status[i,self.Nstars : self.Nstars + new_stars] = self._inv_status_codes['MS']
                self.formation_time[self.Nstars:self.Nstars + new_stars] = self.t[i]
                # need to calculate properties:
                self._new_star_properties(i, self.Nstars, self.Nstars + new_stars)
                self.Nstars = self.Nstars + new_stars
                
            # ----------------------------------------------------------------------------------
            # does the star die?
            new_dead = (self.t[i] > self.lifetime) * (self.t[i-1] < self.lifetime)
            
            MS_stars = (self.status[i-1] == 1) + (self.status[i-1]==0)
            
            SNII   = new_dead * (self.M[i-1] > 8.0)                        * MS_stars
            WD     = new_dead * (self.M[i-1] <= 3.0) * (self.M[i-1] > 1.0) * MS_stars
            WDSN1a = new_dead * (self.M[i-1] <=8.0) * (self.M[i-1] > 3.0)  * MS_stars

            # Salaris et. al. 2009	
            self.M[i][(WD + WDSN1a) * (self.M[i]<4.0)] = 0.134 * self.M[i][(WD + WDSN1a) * (self.M[i]<4.0)] + 0.331
            self.M[i][(WD + WDSN1a) * (self.M[i]>4.0)] = 0.047 * self.M[i][(WD + WDSN1a) * (self.M[i]>4.0)] 

            self.M[i][SNII] = 0.0 # gone

            self.status[i][WD    ] = self._inv_status_codes['WD']
            self.status[i][WDSN1a] = self._inv_status_codes['WD_SN1a']
            self.status[i][SNII  ] = self._inv_status_codes['SNII']

            # shut of radiation for dead stars
            #self.q0[new_dead] = 0.0
            #self.q1[new_dead] = 0.0

            # -----------------------------
            # do stellar wind junk here

            # -----------------------------
            
            # check for SN1a, but not for the newly formed WD_SN1a
            if any( s == self._inv_status_codes['WD_SN1a'] for s in self.status[i-1]):

                WD = (self.status[i-1] == self._inv_status_codes['WD_SN1a'])

                eta   = 0.043 # maoz
                beta  = 1.2 # maoz
                power = -beta + 1
                t_h   = 4.55114687E17 # in seconds... using Planck H = 67.8
                # t_min should be lifetime of most massive star that can form a supernova
                # t_min = self.lifetime[ np.argmin(np.abs(self.M0 - 8.0 )) ]
                
                # if np.size(t_min) > 1:
                #    t_min = t_min[0]
                
                t_o          = self.lifetime[WD] #+ self.lifetime[WD] # time the WD formed
                time_elapsed = self.t[i] - self.formation_time[WD]
                #self.lifetime[WD]# - (t_o)             # time since WD has formed
                
                
                probability = eta * (-beta+1) / (((t_h + t_o)**(-beta+1) - (t_o)**(-beta+1))) * time_elapsed**(-beta)
                probability = probability * dt
                print probability
                #probability = np.ones(np.size(self.M[i][WD])) * probability * dt
                # progenitor masses
                #prog     = self.M[0][WD]

                #IMF_norm = prog**(-2.35)
                
                #msf = np.zeros(np.size(prog))

                #for j in np.arange(np.size(prog)):
           #         if prog[i] < 0.0: # set to zero, makes probability zero
            #            msf[i] = 0.0
             #       else:
                #    msf[j] = IMF_norm[j] * ( imf.integrate_imf(self.IMF_function, self.M_min, self.M_max) )

                #time_elapsed = self.t[i] - self.t[0]

                #probability = imf.DTD(time_elapsed) * msf * (dt / yr)

	

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

def unique(starpop, particle_id, scalar=False):
    
    # loop over every particle and sum
    number = 0
    for i in np.arange(0,int(starpop.Nstars)):
        if any(starpop.status[:,i] == starpop._inv_status_codes[particle_id]):
            number = number + 1
            
    return number
    
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
