"""

  Author : A. Emerick
  Assoc  : Columbia University - American Museum of Natural History
  Date   : April 2016


  Code to extract (most) of the SED from binned OSTAR2002 data
  and repack with first column as the frequency, and each subsequent
  column a separate OSTAR2002 model. Due to the silliness of the 
  SED files, this cuts out the last three frequencey bins in every
  case, corresponding to very low energy photons (well outside the UV
  region we care about)

  This should be run on the OSTAR2002 Cloudy-binned models

"""

import numpy as np
import sys                                     # Z/Z_sun

def extract_sed(filepath = './ostar2002_sed'):
    """
    Extract the OSTAR2002 SED's for the listed identifier
    names in a format thats a bit easier to handle
    """
    identifier_names = ['p03', # 2
                        'p00', # 1
                        'm03', # 1/2
                        'm07', # 1/5
                        'm10', # 1/10
                        'm15', # 1/30
                        'm17', # 1/50
                        'm20', # 1/100
                        'm30'] # 1/1000
#                    'm99'] -- don't use z = 0

    header_count   = 29 
    full_rows      = 4000 # (4000 for all, but 1000 for Z=0) each model takes up this many rows in the file
    rows_per_model = full_rows - 1 # cut out last row since it only has 3 col

    nbins  = rows_per_model * 5   # number of bins
    nmodel = 69                   # number of OSTAR2002 models


    for identifier in identifier_names:
        OSTAR_file  = filepath + 'ostar2002_' + identifier + '.ascii'
        outname     = filepath + 'ostar2002_' + identifier + '_extracted.dat'

        SED_data = np.zeros( (nmodel + 1, nbins) ) # SED for each star, first bin is freq

        # load frequency
        data = np.genfromtxt( OSTAR_file, max_rows = rows_per_model, skip_header = header_count)

        SED_data[0] = data.flatten()


        print("extracting data for ostar file", OSTAR_file)
        # loop over everything else:
        for i in np.arange(nmodel):
            data = np.genfromtxt(OSTAR_file, skip_header = (full_rows)*(i+1) + header_count,
                                             max_rows    = (rows_per_model)  )

            SED_data[i + 1] = data.flatten()

        # now write out to file
        np.savetxt(outname, np.transpose(SED_data), fmt = "%.5E")

    return

if __name__ == "__main__":


    if len(sys.argv) > 1:
        extract_sed(filepath=sys.argv[1])
    else:
        extract_sed()


