import numpy as np
import FUV_heating as FUV

#
# work function values from cloudy (in Ryd)
#
W_graphite = 3.235E-1
W_silicate = 5.883E-1

#
# cloudy dust files (size grain distribution and absorption crs)
#
cloudy_dir = "/home/emerick/Research/stars/radiation/data/cloudy_dust/"
ostar_dir  = "/home/emerick/Research/stars/radiation/data/ostar2002_sed/"

graphite_file_szd = cloudy_dir + "graphite_szd_cloudy.dat"
graphite_file_abs = cloudy_dir + "graphite_ism_01.opc"

silicate_file_szd = cloudy_dir + "silicate_szd_cloudy.dat"
silicate_file_abs = cloudy_dir + "silicate_ism_01.opc"

#
# load szd data once
#
szd_graphite = np.genfromtxt(graphite_file_szd, skip_header = 1)
crs_graphite = np.genfromtxt(graphite_file_abs, skip_header = 1)

szd_silicate  = np.genfromtxt(silicate_file_szd, skip_header = 1)
crs_silicate  = np.genfromtxt(silicate_file_abs, skip_header = 1)

#
# define OSTAR2002 radiation files and file properties
#
# stars are binned by teff and g like follows, but not all 
# g are computed for each temperature !!!
Teff = np.arange(27500.0, 55001.0, 2500.0)
g    = np.arange(3.0000, 4.75001, 0.25)

# dictionary of metallicities and SED file identifiers
Z    = { 'm99' : 0.000, 'm30' : 0.001, 'm20' : 0.01, 'm17' : 1.0/50.0,
         'm15' : 1.0/30.0, 'm10' : 0.10, 'm07' : 0.20, 'm03' : 0.50, 'p00' : 1.00, 'p03' : 2.00 }


#
# compute N grains once for graphite and silicates
#
N_grains_silicate = FUV.integrate_szd(szd_silicate[:,0], szd_silicate[:,1])
N_grains_graphite = FUV.integrate_szd(szd_graphite[:,0], szd_graphite[:,1])

#N_grains_silicate = 2.412114E-21
#N_grains_graphite = 2.084718E-21

print "N_grains_si = %5.5E ---- N_grains_graph = %5.5E"%(N_grains_silicate, N_grains_graphite)


#
# each metallicity gets a column, each T , g combination a row
# values < 0 for heating rate represent models NOT on the OSTAR grid
#

for metal_id in Z:
    ostar_file  = ostar_dir + 'ostar2002_' + metal_id + '_extracted.dat'

    output_file = ostar_dir + 'ostar2002_' + metal_id + '_FUV_heating.dat' 

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
               norm_heating = -9999 # flag that model is not on grid                
            else:
                ostar_data = np.genfromtxt(ostar_file, usecols = (0, i + 1) )

                # compute normalized FUV heating rate from both silicate and graphite grains
                norm_heating =\
                    N_grains_silicate * FUV.integrate_absorption_general(crs_silicate[:,0],
                                                                         crs_silicate[:,1],
                                                                         ostar_data[:,0],
                                                                         ostar_data[:,1],
                                                                         W_graphite, method='simps') +\
                    N_grains_graphite * FUV.integrate_absorption_general(crs_graphite[:,0],
                                                                         crs_graphite[:,1],
                                                                         ostar_data[:,0],
                                                                         ostar_data[:,1],
                                                                         W_silicate, method='simps')

                i = i + 1


            #
            # write this to the output file
            #
            print "T = %5.5f -- g =  %3.3f -- Heat =  %6.6E"%(T,gval,norm_heating)
            f.write("%5.5f %3.3f %6.6E\n"%(T, gval, norm_heating))



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

    FUV_file        = ostar_dir + 'ostar2002_' + metal_id + '_FUV_heating.dat'
    data_array[:,i] = np.genfromtxt(FUV_file, usecols=2)

    i = i + 1



np.savetxt(ostar_dir + 'ostar2002_FUV_all_models.dat', data_array,
                 fmt = "%5.f %3.3f %6.6E %6.6E %6.6E %6.6E %6.6E %6.6E %6.6E %6.6E %6.6E %6.6E")
    
