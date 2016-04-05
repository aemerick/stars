from scipy import integrate

import numpy as np

# just a bunch of IMF functions

def kroupa(M, alpha_1 = 0.3, alpha_2 = 1.3, alpha_3 = 2.3, xi_o = 1.0):
    M = np.asarray(M)
    scalar_input = False
    if M.ndim == 0:
        M = M[None]
        scalar_input = True

    low_mass = M[ (M <= 0.08) ]
    mid_mass = M[ (M  > 0.08)*(M <= 0.5) ]
    salpeter = M[ (M  > 0.5 ) ]

    dNdM = np.zeros(np.shape(M))

    dNdM[low_mass] = M[low_mass]**(-alpha_1)
    dNdM[mid_mass] = M[mid_mass]**(-alpha_2)
    dNdM[salpeter] = M[salpeter]**(-alpha_3)

    if scalar_input:
        return np.squeeze(dNdM)
    else:
        return dNdM

def salpeter(M, alpha = 2.35, xi_o=1.0):
    return xi_o * M**(-alpha)


def sample_IMF(N, IMF,  M_min = 1.0, M_max = 100.0, npoints = 1000, **kwargs):

    # bin IMF in logspace
    dm = np.log10(M_max / M_min) / (1.0*(npoints - 1))

    # m_o
    m_o = np.log10(M_min)

    i = np.arange(0, npoints)
    # cumulative probability density
    m        = 10.0**(m_o + i*dm)
    
    IMF_vals = IMF(m, **kwargs)
    IMF_vals = np.cumsum(IMF_vals)
    IMF_vals = IMF_vals / (IMF_vals[-1] * 1.0)

    random_numbers = np.random.rand(N)
    # now sample

    # do a bisect search for each number
    mass_sample = np.zeros(N)
    for i in np.arange(N):

        bin_number = npoints / 2
        width      = npoints / 2

        while ( width > 1):
            width = width / 2
            if (random_numbers[i] > IMF_vals[bin_number]):
                bin_number = bin_number + width
            elif (random_numbers[i] < IMF_vals[bin_number]):
                bin_number = bin_number - width
            else:
                break
        mass_sample[i] = 10.0**(bin_number *dm)
            

    return mass_sample
    
    
     
def integrate_imf(function, a, b):

    # integrate a function from a to be

    return integrate.quad(function, a, b)[0]


def prob_SN1a(M_proj, t, t_form, dt,
              t_min = 40.0, t_max = 14.6E3 , N_1a_per_msun = 3.0E-3,
              alpha = 2.35, beta  = 1.06, M_min = 1.0, M_max = 100.0):

    t = t - t_form

    norm_DTD = (-beta+1)*(t_max**(-beta + 1.0) - t_min**(-beta+1.0))**(-1.0)*N_1a_per_msun

    M_SF = M_proj**(alpha) * (M_max**(-alpha+2.0) - M_min**(-alpha+2.0)) / (-alpha + 2.0)

    
    return norm_DTD * M_SF * t**(-beta) * dt

def DTD(t, beta = 1.2, NSN = 0.00130):

    Gyr = 3.1556E13 * 1.0E3
    yr  = Gyr / 1.0E9

    return NSN * (t / Gyr)**(-beta) 
