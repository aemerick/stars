import numpy as np


def salpeter(M, alpha=-2.35, xi_o=1.0):
    return xi_o * M**(alpha)


def sample_IMF(N, IMF,  M_min = 0.1, M_max = 100.0, npoints = 1000, **kwargs):

    # bin IMF in logspace
    dm = np.log10(M_max / M_min) / (1.0*(npoints - 1))

    # m_o
    m_o = np.log10(M_min)

    i = np.arange(0, npoints)
    # cumulative probability density
    m        = 10.0**(m_o + i*dm)
    total_fn = IMF(m) 
    

    # determine mass sample points
    M_sample = np.logspace(M_min, M_max, npoints)
 
    # sample the IMF 
    IMF_vals = IMF(M_sample, **kwargs)
    IMF_vals = IMF_vals / (1.0*np.sum(IMF_vals))
    IMF_vals = np.cumsum(IMF_vals)

    # normalize to one


    random_numbers = np.random.rand(N)
    # now sample



    # do a bisect search for each number
    mass_sample = np.zeros(N)
    bin_number = npoints / 2
    width = npoints / 2 
    for i in np.arange(N):

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
    
    
     



def prob_SN1a(M_proj, t, t_form, dt,
              t_min = 40.0, t_max = 14.6E3 , N_1a_per_msun = 3.0E-3,
              alpha = 2.35, beta  = 1.06, M_min = 1.0, M_max = 100.0):

    t = t - t_form

    norm_DTD = (-beta+1)*(t_max**(-beta + 1.0) - t_min**(-beta+1.0))**(-1.0)*N_1a_per_msun

    M_SF = M_proj**(alpha) * (M_max**(-alpha+2.0) - M_min**(-alpha+2.0)) / (-alpha + 2.0)

    
    return norm_DTD * M_SF * t**(-beta) * dt

