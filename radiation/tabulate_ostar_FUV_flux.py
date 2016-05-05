import numpy as np
import FUV_heating as FUV

#
# define OSTAR2002 radiation files and file properties
#
ostar_dir  = "/home/emerick/Research/stars/radiation/data/ostar2002_sed/"

# stars are binned by teff and g like follows, but not all 
# g are computed for each temperature !!!
Teff = np.arange(27500.0, 55001.0, 2500.0)
g    = np.arange(3.0000, 4.75001, 0.25)

# dictionary of metallicities and SED file identifiers
Z    = { 'm99' : 0.000, 'm30' : 0.001, 'm20' : 0.01, 'm17' : 1.0/50.0,
         'm15' : 1.0/30.0, 'm10' : 0.10, 'm07' : 0.20, 'm03' : 0.50, 'p00' : 1.00, 'p03' : 2.00 }

for metal_id in Z:
    ostar_file  = ostar_dir + 'ostar2002_' + metal_id + '_extracted.dat'

    output_file = ostar_dir + 'ostar2002_' + metal_id + '_FUV_flux.dat'


    f = open(output_file, 'w')
    f.write("# T logg FUV\n")

    # loop over each model
    i = 0
    for T in Teff:
        for gval in g:
            # skip models that don't exist:
            if ( gval == 3.0  and T >= 32500) or\
               ( gval == 3.25 and T >= 37500) or\
               ( gval == 3.50 and T >= 42500) or\
               ( gval == 3.75 and T >= 50000):
               FUV_flux = 0.0 # flag that model is not on grid                
            else:
                ostar_data = np.genfromtxt(ostar_file, usecols = (0, i + 1) )

                # compute normalized FUV heating rate from both silicate and graphite grains
                FUV_flux = FUV.integrate_SED(ostar_data[:,0], ostar_data[:,1])

                i = i + 1


            print "T = %5.5f -- g = %3.3f -- Flux = %8.8E\n"%(T, gval, FUV_flux)
            f.write("%5.5f %3.3f %6.6E\n"%(T, gval, FUV_flux))



    f.close()
    # star alllll over again


# now grab and merge all the files
nrows = np.size(Teff) * np.size(g)
ncol  = 10 + 2 # 10 Z's, 1 for each T and g


data_array = np.zeros((nrows,ncol))

data_array[:,0] = np.array(np.sort(list(Teff) * np.size(g)))
data_array[:,1] = np.array(list(g) * np.size(Teff))

i = 2
for metal_id in Z:

    FUV_file        = ostar_dir + 'ostar2002_' + metal_id + '_FUV_flux.dat'
    data_array[:,i] = np.genfromtxt(FUV_file, usecols=2)

    i = i + 1



np.savetxt(ostar_dir + 'ostar2002_FUV_flux_all_models.dat', data_array,
                 fmt = "%5.f %3.3f %6.6E %6.6E %6.6E %6.6E %6.6E %6.6E %6.6E %6.6E %6.6E %6.6E")
    
