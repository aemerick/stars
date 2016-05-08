"""
    Script to extract stellar yields from a given
    NuGrid data set to translate them into an Enzo-readable
    table. This relies on Christian Ritter's read_yields.py

    Author: Andrew Emerick
    Date  : 04/2016
"""
from chemistry.NuGrid import read_yields as ry
import numpy as np


#
master_element_list = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na',
                       'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti',
                       'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As',
                       'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru',
                       'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs',
                       'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy',
                       'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir',
                       'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi'] # 

# make it so they are unique
for i in range(np.size(master_element_list)):
    if len(master_element_list[i]) == 1:
        master_element_list[i] = master_element_list[i] + '-'

#
# need to make a table containing following columns
# M    Z    combined yields
#

NuGrid_dir  = '/home/emerick/code/chemistry/NuGrid/'
NuGrid_dset = 'set2'
yield_table = 'isotope_yield_table_MESA_only_fryer12_delay'

table = ry.read_nugrid_yields(NuGrid_dir + NuGrid_dset + '/' + yield_table + '.txt')
wind_table = ry.read_nugrid_yields(NuGrid_dir + NuGrid_dset + '/' + yield_table + '_winds.txt')

# load available metallicities and masses
# have to hard code this because they aren't stored... dumb
allM = np.array([1.0, 1.65, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 12.0, 15.0, 20.0, 25.0])
allZ = table.metallicities[::-1]


# now loop through and extract metallicities for each

# what isotopes are availabel?
isotopes = table.get(allM[0], allZ[0], 'Isotopes')
# get a list of combined isotopes:
#abridged = [s[:2] for s in isotopes]
#indexes  = np.unique( abridged, return_index=True )[1]
#elements = [abridged[i] for i in sorted(indexes)]
#elements = np.array(['H-','Mg','Fe'])

elements  = master_element_list

fmt_string = "%5.5E  %5.5E " + "%5.5E "*np.size(elements) + "\n"
fmt = "%5.5E "

f     = open('stellar_yields.dat','w')
fwind = open('stellar_yields_winds.dat','w')

for file in [f, fwind]:
    file.write("# M Z mtot")
    for e in elements:
        if '-' in e:
            file.write(e[0] + ' ')
        else:
            file.write(e + ' ')
    file.write('\n')

for m in allM:

    for z in allZ:

        yields      = table.get(  m, z, 'Yields')
        winds_yield = wind_table.get( m, z, 'Yields') 

        # make a dictionary
        yield_dict = dict(zip(isotopes, yields))
        winds_dict  = dict(zip(isotopes, winds_yield))

        # now combine them all 
        yield_combined  = dict(zip(elements, np.zeros(np.size(elements))))
        winds_combined  = dict(zip(elements, np.zeros(np.size(elements))))
        sn_combined     = dict(zip(elements, np.zeros(np.size(elements))))

        for e in elements:
            selected_isotopes = [yield_dict[s] for s in isotopes if e in s]
            yield_combined[e]  = np.sum(np.array(selected_isotopes))
            selected_isotopes = [winds_dict[s]  for s in isotopes if e in s]
            winds_combined[e]  = np.sum(np.array(selected_isotopes))

            sn_combined[e] = yield_combined[e] - winds_combined[e]
        # find total stellar mass
        
       
        # we now have a dictionary with keys elements
        # and values summed yields for each elemental species
        total_mass = np.sum(sn_combined.values())
        metal_mass = total_mass - sn_combined['H-'] - sn_combined['He']

        f.write("%.3f %.5f "%(m,z)); f.write(fmt%(total_mass)); f.write(fmt%(metal_mass));

        total_mass = np.sum(winds_combined.values())
        metal_mass = total_mass - winds_combined['H-'] - winds_combined['He']

        fwind.write("%.3f %.5f "%(m,z)); fwind.write(fmt%(total_mass)); fwind.write(fmt%(metal_mass));
        for e in elements:
            f.write(    fmt%(sn_combined[e]))
            fwind.write(fmt%(winds_combined[e]))

            if (( yield_combined[e] - winds_combined[e] ) < 0.0):
                print "WARNING: Negative Yields", m, z, e

        f.write('\n')
        fwind.write('\n')

f.close()
fwind.close()        
            
f = open('column_key.dat','w')

f.write("# element col\n")
i = 3
for e in elements:
    if '-' in e:
        name = e[0] + ' '
    else:
        name = e

    f.write(name + " %i\n"%(i))
    i = i + 1
f.close()
