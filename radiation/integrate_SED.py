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

HI_emin = 13.6  # eV
HI_emax = np.inf

HeI_emin = 24.6   # eV
HeI_emax = np.inf

HeII_emin = 54.42  # eV
HeII_emax = np.inf

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
        #   this needs to be in high-to-low order
        self.metallicity_names = { 'p03' : 2.00,'p00':1.0,'m03':0.50,
                                   'm07' : 0.20, 'm10':0.10, 'm15' : 1.0/30.0,
                                   'm17' : 1.0/50.0, 'm20' : 0.01, 'm30' : 0.001,
                                   'm99' : 0.000}

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

    def construct_table(self, emin, emax, name = None, precision = 4,
                        flux_type = 'flux', output_log = False,
                        output_name = None):
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

        flux_type : string, optional
            Output table as either energy flux (erg/s/cm^2) or photon
            flux (1/s/cm^2). Default is energy flux. See 
            'integrate_SED' for more details. Default : 'flux'

        output_log : bool, optional
            Output logged fluxes. Default: False

        output_precision : int, optional
            Precision of output fluxes. Default : 4

        output_name : string, optional
            If provided, overrides default naming scheme with this string
            for the output file. Default : None

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
                                             emin=emin, emax=emax, flux_type = flux_type)
                    i = i + 1
            #
            # f.write("%5.5f %3.3f %6.6E\n"%(Tval,gval, flux))


        # now grab and merge these tables
        nrows = np.size(self.Teff) * np.size(self.g)
        ncol  = np.size(self._z_ids) + 2 # num of metallicities, Teff, and g

        data_array = np.zeros((nrows,ncol))
        data_array[:,0] = np.array(np.sort(list(self.Teff)*np.size(self.g)))
        data_array[:,1] = np.array(list(self.g)*np.size(self.Teff))

        for i,metal_id in enumerate(self._z_ids):
            #flux_file = self.filepath + 'ostar2002_'+metal_id+name+'_flux.dat'
            data_array[:,i+2] = flux_data[metal_id]

        fmt = "%5.f %3.3f "
        prec_string = "%" + "%i.%i"%(precision+2,precision)

        if output_log:
            fmt = fmt + (" " + prec_string + "f")*np.size(self._z_ids)
        else:
            fmt = fmt + (" " + prec_string + "E")*np.size(self._z_ids)


        if output_name is None:
            output_name = 'ostar2002' + name + '_flux_all_models.dat'

        np.savetxt(self.filepath + output_name, data_array, fmt = fmt)

        return

    def integrate_SED(self, freq, SED, emin, emax, flux_type = 'flux'):
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

        emax : float
            Maximum energy of desired band to integrate (in eV)

        flux_type : string, optional
            Output type. Defaults to returning the energy flux (erg/s/cm^2)
            if 'flux' or 'energy flux' or 'energy' is provided. Otherwise
            'photon' or 'photon flux' will return the photon counts
            (1/s/cm^2). Default : 'flux'
        """

        nu_min = ((emin * u.eV) / (const.h.cgs)).to(u.Hz).value
        nu_max = ((emax * u.eV) / (const.h.cgs)).to(u.Hz).value

        if freq[-1] < freq[0]: # sort in ascending order
            freq = freq[::-1]
            SED  = SED[::-1]

        selection = (freq<=nu_max)*(freq>nu_min)


        if flux_type in ['flux','energy','energy flux']:

            flux = np.trapz(SED[selection], x = freq[selection])

        elif flux_type in ['photon','photon flux']:

            photon_energy = ((freq[selection]*u.Hz)*const.h).to(u.erg).value

            flux = np.trapz(SED[selection]/photon_energy, x = freq[selection])

        return 4.0 * np.pi * flux


def compute_all_tables(filepath):
    """
    Build all tables
    """

    ostar = OSTAR_SED(filepath)

    precision   = 4
    output_log  = False

    print("Building tables needed for individual star model in Enzo")
    print("by default, these tables compute the IR, FUV, and LW band fluxes (erg/s/cm^2)")
    print("and the HI, HeI, HeII ionizing photon fluxes (1/s/cm^2). This may take a few minutes.")

    ostar.construct_table(FUV_emin, FUV_emax, output_name = "FUV_energy_rates.in",
                          precision=precision,output_log=output_log)
    print("FUV complete")

    ostar.construct_table(LW_emin, LW_emax, output_name = "LW_energy_rates.in",
                          precision=precision,output_log=output_log)
    print("LW complete")

    ostar.construct_table(IR_emin, IR_emax, output_name = "IR_energy_rates.in",
                          precision=precision,output_log=output_log)
    print("IR complete")

    ostar.construct_table(HI_emin, HI_emax, output_name = "q0_photon_rates.in",
                          flux_type = 'photon',
                          precision=precision,output_log=output_log)
    print("HI complete")

    ostar.construct_table(HeI_emin, HeI_emax, output_name = "q1_photon_rates.in",
                          flux_type = 'photon',
                          precision=precision,output_log=output_log)
    print("HeI complete")

    ostar.construct_table(HeII_emin, HeII_emax, output_name = "q2_photon_rates.in",
                          flux_type = 'photon',
                          precision=precision,output_log=output_log)
    print("HeII complete")

    return

if __name__ == '__main__':


    if len(sys.argv) > 1:
        compute_all_tables(sys.argv[1])
    else:
        compute_all_tables('./')
