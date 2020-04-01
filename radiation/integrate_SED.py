import numpy as np
from astropy import units as u
from astropy import constants as const
from scipy import integrate
from scipy import interpolate

import glob
import sys, os

IR_emin   = 0.76  # eV
IR_emax   = 5.6   # eV

FUV_emin  = 5.6   # eV
FUV_emax  = 11.2  # eV

LW_emin   = 11.2  # eV
LW_emax   = 13.6  # eV

class OSTAR_SED:

    def __init__(self, filepath, name=None):
        """
        Class to integrate OSTAR 2002 SED's and construct tables
        in general format that can be used readily in other applications.
        Handles reading in OSTAR tables (assuming processed with extract_sed.py)
        and compiling models at all available T, g, and Z values into a single
        table in a given band.
        """

        # OSTAR metallicity file names
        self.metallicity_names = { # 'm99' : 0.000, # ignore z = 0 file
                                   'm30' : 0.001, 'm20' : 0.01,
                                   'm17' : 1.0/50.0, 'm15' : 1.0/30.0,
                                   'm10' : 0.10, 'm07' : 0.20, 'm03' : 0.50,
                                   'p00' : 1.00, 'p03' : 2.00 }
        self._z_ids = list(self.metallicity_names.keys())
        self._z_vals = list(self.metallicity_names.values())

        for k in self._z_ids:
            self.metallicity_names[ self.metallicity_names[k] ] = k


        #
        self.filepath = filepath

        if not os.path.isdir(self.filepath):
            print(self.filepath, " not valid")
            raise ValueError

        # check to make sure files exist
        for zname in self._z_ids:
            filename = 'ostar2002_'+zname+'_extracted.dat'
            if not os.path.isfile(self.filepath + filename):
                print(filename, " not found in ", self.filepath)
                raise ValueError

        self.name = name

        # full Teff and g available in grid
        # however, some exclusions apply
        self.Teff = np.arange(27500.0, 55001.0, 2500.0)
        self.g    = np.arange(3.0, 4.75001, 0.25)

        # this is a bit of a hack, but given exclusions
        # maps column numbers in extracted SED files to for
        # each Teff,g combination.
        i = 0
        self._grid_mask = {}
        for Tval in self.Teff:
            for gval in self.g:
                if self.check_excluded(Tval,gval):
                    self._grid_mask[(Tval,gval)] = -1
                else:
                    self._grid_mask[(Tval,gval)] = i
                    i = i + 1

        return

    def check_excluded(self, Tval, gval):
        """
        The OSTAR2002 grid has some T and g exlusion regions
        Check to see if this point is in the grid or not.
        If it is, returns true, otherwise false
        """

        is_excluded =  ((gval==3.00) and (Tval>=32500.0)) or\
                       ((gval==3.25) and (Tval>=37500.0)) or\
                       ((gval==3.50) and (Tval>=42500.0)) or\
                       ((gval==3.75) and (Tval>=50000.0))

        return is_excluded


    def load_data(self, metallicity, Tval=None, gval=None,
                        ignore_error = False):
        """
        Load one of the OSTAR2002 extracted files. First column
        is the frequency, while the rest give the SED for a
        certain pair of T and g values. If Tval and gval
        are specified, this will only load that column.

        Parameters:
        -----------
        metallicity : float or string
            Metallicity (either the OSTAR2002 ID name or actual value)
            of desired table to load. Must be exactly equal to an
            available table

        Tval : float, optional
            The whole table is loaded unless Tval and gval are specified.
            Tval is the desired effective temperature of model. Default : None

        gval : float, optional
            The whole table is loaded unless Tval and gval are specified.
            gval is the desired surface gravity (in log(g)). Default : None

        ignore_error : bool, optional
            By default, routine throws an error if Tval and gval are not
            a valid grid point in the OSTAR2002 grids. Set to True to
            ignore this and return None instead of the data table. Default : False

        Returns:
        --------
        ostar_data : np.ndarray
            2D array containing frequencies and SED (flux) in cgs

        """

        if metallicity in self._z_vals:
            mname = self.metallicity_names[metallicity]
        elif metallicity in self._z_ids:
            mname = metallicity
        else:
            raise ValueError

        ostar_file = self.filepath + 'ostar2002_' + mname + '_extracted.dat'

        if (Tval is None) and (gval is None):
            ostar_data = np.genfromtxt(ostar_file)
        else:
            colnum     = self._grid_mask[(Tval,gval)]

            if (colnum < 0) and (not ignore_error):
                print("T = %5.5E  g = %5.5E are not available"%(Tval,gval))
                raise ValueError
            elif (colnum < 0) and ignore_error:
                return None

            ostar_data = np.genfromtxt(ostar_file, usecols=(0,colnum+1))

        return ostar_data

    def construct_table(self, emin, emax, name = None):
        """
        Build a single table containing the total flux
        in a given energy band.

        Parameters
        ----------
        emin : float
            Minimum energy in desired band (in eV)

        emax : float
            Maximum energy in desired band (in eV)

        name : string, optional
            Output naming string (FUV, IR, LW, etc.) to
            describe output file. output file will be of 
            format: 'ostar2002_'+name+'_flux_all_models.dat' Default : None

        Returns
        -------
        void
        """

        if name is None:
            name =''
        else:
            name = '_' + name

        flux_data = {}
        for metal_id in self._z_ids:

            flux_data[metal_id] = np.zeros(np.size(self.Teff)*np.size(self.g))

            # make (unfortunately) a bunch of intermediate files
            # f = open(self.filepath + 'ostar2002_'+metal_id+name+'_flux.dat'
            # f.write("# T logg FUV\n")
            i = 0
            for Tval in self.Teff:
                for gval in self.g:

                    ostar_data = self.load_data(metal_id, Tval, gval,
                                                ignore_error=True)

                    if ostar_data is None:
                        flux_data[metal_id][i] = 0.0
                    else:
                        flux_data[metal_id][i] = self.integrate_SED(ostar_data[:,0],
                                             ostar_data[:,1],
                                             emin=emin, emax=emax)
                    i = i + 1
            #
            # f.write("%5.5f %3.3f %6.6E\n"%(Tval,gval, flux))


        # now grab and merge these tables
        nrows = np.size(self.Teff) * np.size(self.g)
        ncol  = 9 + 2

        data_array = np.zeros((nrows,ncol))
        data_array[:,0] = np.array(np.sort(list(self.Teff)*np.size(self.g)))
        data_array[:,1] = np.array(list(self.g)*np.size(self.Teff))

        for i,metal_id in enumerate(self._z_ids):
            #flux_file = self.filepath + 'ostar2002_'+metal_id+name+'_flux.dat'
            data_array[:,i+2] = flux_data[metal_id]

        np.savetxt(self.filepath + 'ostar2002' + name +'_flux_all_models.dat',
                   data_array,
                   fmt = "%5.f %3.3f %6.6E %6.6E %6.6E %6.6E %6.6E %6.6E %6.6E %6.6E %6.6E")

        return

    def integrate_SED(self, freq, SED, emin, emax):
        """
        Integrates an SED.

        Parameters
        ----------
        freq : np.ndarray
            Frequency (in Hz). Check is made to ensure this is 
            in ascending order

        SED : np.ndarray
            SED flux (in cgs)

        emin : float
            Minimum energy of desired band to integrate (in eV)

        emax : fload
            Maximum energy of desired band to integrate (in eV)
        """

        nu_min = ((emin * u.eV) / (const.h.cgs)).to(u.Hz).value
        nu_max = ((emax * u.eV) / (const.h.cgs)).to(u.Hz).value

        if freq[-1] < freq[0]: # sort in ascending order
            freq = freq[::-1]
            SED  = SED[::-1]

        selection = (freq<=nu_max)*(freq>nu_min)

        return 4.0 * np.trapz(SED[selection], x = freq[selection])


def compute_all_tables(filepath):


    ostar = OSTAR_SED(filepath)

    print("Building tables needed for individual star model in Enzo")
    ostar.construct_table(FUV_emin, FUV_emax, name = 'FUV')
    print("FUV complete")
    ostar.construct_table(LW_emin, LW_emax, name = 'LW')
    print("LW complete")
    ostar.construct_table(IR_emin, IR_emax, name = 'IR')
    print("IR complete")


    return

if __name__ == '__main__':


    if len(sys.argv) > 1:
        compute_all_tables(sys.argv[1])
    else:
        compute_all_tables('./')
