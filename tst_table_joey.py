import swiftsimio as sw
import numpy as np
import interpolate_X_Ray_seperate_lines_cont as ips
from unyt import mp, cm

def load_snapshot():
    filename = "/cosma6/data/dp004/dc-kuge1/200LH_bestfit_tests/FLAMFIDL200N360Nearest_x_ray/bahamas_0008.hdf5"
    mask = sw.mask(filename)
    # The full metadata object is available from within the mask
    boxsize = mask.metadata.boxsize
    # load_region is a 3x2 list [[left, right], [bottom, top], [front, back]]
    # We load a small region
    load_region = [[0.0 * b, 0.001 * b] for b in boxsize]
    # Constrain the mask
    mask.constrain_spatial(load_region)
    # Now load the snapshot with this mask
    data = sw.load(filename, mask=mask)
    return data

def update_data_structure(data, part_lc):
    data.gas.densities = np.power(10, np.array([-1.6, -1.6])) * cm**-3 * mp / data.gas.smoothed_element_mass_fractions.hydrogen[0:2]
    data.gas.temperatures = np.power(10, np.array([6.0, 6.0]))
    data.gas.masses = data.gas.masses[0:2]
    data.gas.smoothed_element_mass_fractions.hydrogen = data.gas.smoothed_element_mass_fractions.hydrogen[0:2]
    data.gas.smoothed_element_mass_fractions.helium = data.gas.smoothed_element_mass_fractions.helium[0:2]
    data.gas.smoothed_element_mass_fractions.carbon = data.gas.smoothed_element_mass_fractions.carbon[0:2]
    data.gas.smoothed_element_mass_fractions.nitrogen = data.gas.smoothed_element_mass_fractions.nitrogen[0:2]
    data.gas.smoothed_element_mass_fractions.oxygen = data.gas.smoothed_element_mass_fractions.oxygen[0:2]
    data.gas.smoothed_element_mass_fractions.neon = data.gas.smoothed_element_mass_fractions.neon[0:2]
    data.gas.smoothed_element_mass_fractions.magnesium = data.gas.smoothed_element_mass_fractions.magnesium[0:2]
    data.gas.smoothed_element_mass_fractions.silicon = data.gas.smoothed_element_mass_fractions.silicon[0:2]
    data.gas.smoothed_element_mass_fractions.iron = data.gas.smoothed_element_mass_fractions.iron[0:2]



    # set solar metallicity
    data.gas.densities = np.power(10, np.array([-1.6, -1.6])) * cm**-3 * mp / np.array([7.37873810e-01, 7.37873810e-01])
    data.gas.smoothed_element_mass_fractions.hydrogen = np.array([7.37873810e-01, 7.37873810e-01])
    data.gas.smoothed_element_mass_fractions.helium = np.array([2.49406147e-01, 2.49406147e-01])
    data.gas.smoothed_element_mass_fractions.carbon = np.array([2.36627149e-03, 2.36627149e-03])
    data.gas.smoothed_element_mass_fractions.nitrogen = np.array([6.93361461e-04, 6.93361461e-04])
    data.gas.smoothed_element_mass_fractions.oxygen = np.array([5.73642799e-03, 5.73642799e-03])
    data.gas.smoothed_element_mass_fractions.neon = np.array([1.25731103e-03, 1.25731103e-03])
    data.gas.smoothed_element_mass_fractions.magnesium = np.array([7.08445625e-04, 7.08445625e-04])
    data.gas.smoothed_element_mass_fractions.silicon = np.array([6.65385454e-04, 6.65385454e-04])
    data.gas.smoothed_element_mass_fractions.iron = np.array([1.29283988e-03, 1.29283988e-03]) 

    data.gas.redshift = np.zeros(2)
    return data

data = load_snapshot()

data = update_data_structure(data, 0)

for table_type in ['all', 'lines', 'cont']:
    lum, ens = ips.interpolate_X_Ray(data.gas.densities, data.gas.temperatures, data.gas.smoothed_element_mass_fractions, data.gas.redshift, data.gas.masses,\
                        fill_value = None, bin_energy_lims = [0.6, 0.7], table_type = table_type)

    print(table_type, lum)