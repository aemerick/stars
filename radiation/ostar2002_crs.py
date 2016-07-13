import numpy as np
import FUV_heating as FUV
import phfit2

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

#
# Tabulate cross sections for each frequency first
#
c = 2.998E10        # cm / s
h = 4.135667662E-15 # ev * s

Mb_to_cgs = 1.0E-18 # cm**2

v_hi         = 13.59840 / h
v_hi         = v_hi
v_hei        = 24.58740 / h 
v_heii       = 54.51776 / h 
def phfit_wrapper(z, ne, ishell, photon_energy):
    return phfit2.phfit2(z, ne, ishell, photon_energy)

def total_crs(z, ne, photon_energy):
    crs = 0.0

    for ishell in np.arange(1, 7 + 1): # total of 7 shells
        crs += phfit_wrapper(z, ne, ishell, photon_energy)

    return crs * Mb_to_cgs



def avg_e(freq, SED, crs, nu_min):

    selection = (freq >= nu_min)

    num_integrand   = h * freq * SED * crs
    denom_integrand = SED * crs

    num = np.trapz(num_integrand[selection], x = freq[selection])
    denom = np.trapz(denom_integrand[selection], x = freq[selection])

    if denom == 0.0:
        E_avg = 0.0
    else:
        E_avg = num / denom

    return E_avg
    

for metal_id in Z:
    ostar_file  = ostar_dir +  'ostar2002_' + metal_id + '_extracted.dat'

    output_file = ostar_dir + '/E_avg/' + 'ostar2002_' + metal_id + '_average_energy.dat'


    freq = np.genfromtxt(ostar_file, usecols=0)
    freq = freq.flatten()

    nfreq = np.size(freq)

    nrow = np.size(Teff) * np.size(g)
    hi_crs   = np.zeros(nfreq) ; E_avg_hi   = np.zeros(nrow)
    hei_crs  = np.zeros(nfreq) ; E_avg_hei  = np.zeros(nrow)
    heii_crs = np.zeros(nfreq) ; E_avg_heii = np.zeros(nrow)


    ii = 0
    for nu in freq:
        hi_crs[ii]   = total_crs(1, 1, nu * h)
        hei_crs[ii]  = total_crs(2, 2, nu * h)
        heii_crs[ii] = total_crs(2, 1, nu * h)
        ii = ii + 1
    print np.min(hi_crs), np.max(hi_crs), np.min(hei_crs), np.max(hei_crs), np.min(heii_crs), np.max(heii_crs)

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
                print "T = %5.5f -- g = %3.3f -- %.4f %.4f %.4f"%(T, gval, 0.0, 0.0, 0.0)
                f.write("%5.5f %3.3f %.4f %.4f %.4f\n"%(T, gval, 0.0, 0.0, 0.0))
            else:
                ostar_data = np.genfromtxt(ostar_file, usecols = (0, i + 1) )
                
                #
                # compute integrands to pass to numpy
                #
                SED = ostar_data[:, 1]
                E_avg_hi[i] = avg_e(freq, SED, hi_crs, v_hi)
                E_avg_hei[i] = avg_e(freq, SED, hei_crs, v_hei)
                E_avg_heii[i] = avg_e(freq, SED, heii_crs, v_heii)

                print "T = %5.5f -- g = %3.3f -- %.4f %.4f %.4f"%(T, gval, E_avg_hi[i], E_avg_hei[i], E_avg_heii[i])
                f.write("%5.5f %3.3f %.4f %.4f %.4f\n"%(T, gval, E_avg_hi[i], E_avg_hei[i], E_avg_heii[i]))
                i = i + 1



    f.close()
    # star alllll over again


# now grab and merge all the files
nrows = np.size(Teff) * np.size(g)
ncol  = 10 + 2 # 10 Z's, 1 for each T and g


data_array = np.zeros((nrows,ncol))

data_array[:,0] = np.array(np.sort(list(Teff) * np.size(g)))
data_array[:,1] = np.array(list(g) * np.size(Teff))

label = ['E_avg_hi', 'E_avg_hei', 'E_avg_heii']

for j in [2,3,4]:
    i = 2
    for metal_id in Z:

        data_file        = ostar_dir + '/E_avg/ostar2002_' + metal_id + '_average_energy.dat'
        data_array[:,i] = np.genfromtxt(data_file, usecols = j)

        i = i + 1

    np.savetxt(ostar_dir + 'ostar2002_' + label[j - 2] + '_all_models.dat', data_array,
                     fmt = "%5.f %3.3f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f")
    

