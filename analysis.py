from base import *
from ParticleTracking.generate_haloid_file import snapnums_cptmarvel
def load_ejected_particles(key):
    """
    Load the ejected particles data from an HDF5 file.

    Parameters:
    key (str): The key of desired [simname]_[halo number]
        Ex. cptmarvel_1

    Returns:
    pd.DataFrame: The ejected particles data.
    """
    return pd.read_hdf('Data/SNe/ejected_particles.hdf5', key=key)

def load_expelled_particles(key):
    """
    Load the ejected particles data from an HDF5 file.

    Parameters:
    key (str): The key of desired [simname]_[halo number]
        Ex. cptmarvel_1

    Returns:
    pd.DataFrame: The ejected particles data.
    """
    return pd.read_hdf('Data/SNe/expelled_particles.hdf5', key=key)

def load_tracked_particles(key):
    """
    Load the tracked particles data from an HDF5 file.

    Parameters:
    key (str): The key of desired [simname]_[halo number]
        Ex. cptmarvel_1

    Returns:
    pd.DataFrame: The tracked particles data.
    """
    return pd.read_hdf('Data/tracked_particles.hdf5', key=key)

def load_sim(sim,snapnum):
    """
    Load the simulation the data/Sims folder

    Parameters:
    snapnum (int): The snapshot number of the simulation.

    Returns:
    pynbody.snapshot.Snapshot: The loaded simulation snapshot.
    pynbody.halo.Halo: The halos in the simulation.
    """
    path = f'{SIM_FOLDER_PATH}/{sim}.cosmo25cmb/{sim}.cosmo25cmb.4096g5HbwK1BH/snapshots_200crit_{sim}/{sim}.cosmo25cmb.4096g5HbwK1BH.00{snapnum}'
    sim = pynbody.load(path)
    sim.physical_units()

    halos = sim.halos(halo_numbers='v1')
    return sim, halos

def get_haloids(sim, z0halo):
    """
    Get the halo IDs from the simulation.

    Parameters:
    sim (pynbody.snapshot.Snapshot): The loaded simulation snapshot.
    z0halo (int): The halo number at z=0.

    Returns:
    list: A list of halo IDs.
    """
    
    d = pickle.load(open('Data/simpaths_haloids.pickle','rb'))
    ids = d['haloids'][sim][z0halo]
    paths = d['filepaths'][sim][z0halo]
    return ids, paths

def get_halo_mlf(key, snap_num):
    """
    Get the mass loading factor of a halo at z=0.

    Parameters:
    key (str): The key of desired [simname]_[halo number]
        Ex. cptmarvel_1
    snap_num (int): The snapshot number of the simulation.
    delta_T (float): The time interval over which to calculate the star formation rate, in Gyr.
    
    Returns:
    float: The mass loading factor of the halo.
    """
    # loading in coressponding halos
    sim, halos = load_sim(key.split('_')[0], snap_num)
    halo = halos[int(key.split('_')[1])]
    # extracting stellar information
    stars = halo.star
    star_tform = stars['tform'].in_units('Gyr')

    # loading in tracked expelled particles
    all_expelled = load_expelled_particles(key)

    # loading in snapshot numbers and current snapshot time
    snap_numbers = np.array([int(num) for num in globals()[f'snapnums_{key.split("_")[0]}']])
    snap_time = sim.properties['time'].in_units('Gyr')
    
    # deriving expelled mass
    selected_expelled = all_expelled[np.isclose(all_expelled['time'], snap_time, atol=1e-5)]
    expelled_mass = pynbody.units.Unit(f"{selected_expelled['mass'].sum()} Msol")

    # calcuating time from the last snapshot: delta_T
    current_snapnum_index = np.where(snap_numbers == snap_num)[0][0]
    if current_snapnum_index == len(snap_numbers) - 1:
        logger.debug(f"Current snapshot {snap_num} is the last snapshot in the simulation. Cannot calculate mass loading factor.")
        return float('nan')
    else:
        prev_snapnum = snap_numbers[current_snapnum_index + 1]
        delta_T = snap_time - load_sim(key.split('_')[0], prev_snapnum)[0].properties['time'].in_units('Gyr')
        delta_T = pynbody.units.Unit(f"{delta_T} Gyr")
        # checking if delta_T is zero
        if np.isclose(delta_T.in_units('Gyr'), float(0), atol=1e-9):
            logger.debug(f"Delta T is zero for snapshot {snap_num}. Cannot calculate mass loading factor.")
            return float('nan')

    # calculating the mass outflow rate
    mor = np.divide(expelled_mass, delta_T)

    # selecting formed stars within the time interval
    formed_stars = stars[(star_tform <= snap_time) & (star_tform > (snap_time - delta_T.in_units('Gyr')))]
    formed_star_mass = pynbody.units.Unit(f"{formed_stars['mass'].sum().in_units('Msol')} Msol")

    # calculating the star formation rate
    sfr = np.divide(formed_star_mass, delta_T) 

    # calculating the mass loading factor
    if np.isclose(sfr.in_units('Msol Gyr**-1'), float(0), atol=1e-9):
        logger.debug(f"Star formation rate is zero for snapshot {snap_num}. Cannot calculate mass loading factor.")
        return float('nan')
    else:
        # returns dimensionless mass loading factor as raw float
        return np.divide(mor, sfr).ratio(pynbody.units.Unit(1))

def get_halo_mlf_anulus(key, snap_num, r, r_width):
    """
    Get the mass loading factor of a halo at z=0 within a specific radius.

    Parameters:
    key (str): The key of desired [simname]_[halo number]
        Ex. cptmarvel_1

    snap_num (int): The snapshot number of the simulation.
    r (float): The multiple of virial radius within which to calculate the mass loading factor.
    r_width (float): The width of the anulus within which to calculate the mass loading factor, in kpc.

    returns:
    float: The dimensionless loading factor of the halo within the specified anulus.
    """

       # loading in coressponding halos
    sim, halos = load_sim(key.split('_')[0], snap_num)
    halo = halos[int(key.split('_')[1])]
    # extracting stellar information
    stars = halo.star
    star_tform = stars['tform'].in_units('Gyr')

    # loading in tracked expelled particles
    all_expelled = load_expelled_particles(key)

    # loading in snapshot numbers and current snapshot time
    snap_numbers = np.array([int(num) for num in globals()[f'snapnums_{key.split("_")[0]}']])
    snap_time = sim.properties['time'].in_units('Gyr')
    
    # deriving expelled mass
    selected_expelled = all_expelled[np.isclose(all_expelled['time'], snap_time, atol=1e-5)]
    expelled_mass = pynbody.units.Unit(f"{selected_expelled['mass'].sum()} Msol")
    # selecting expelled mass within the specified anulus

    # calcuating time from the last snapshot: delta_T
    current_snapnum_index = np.where(snap_numbers == snap_num)[0][0]
    if current_snapnum_index == len(snap_numbers) - 1:
        logger.debug(f"Current snapshot {snap_num} is the last snapshot in the simulation. Cannot calculate mass loading factor.")
        return float('nan')
    else:
        prev_snapnum = snap_numbers[current_snapnum_index + 1]
        delta_T = snap_time - load_sim(key.split('_')[0], prev_snapnum)[0].properties['time'].in_units('Gyr')
        delta_T = pynbody.units.Unit(f"{delta_T} Gyr")
        # checking if delta_T is zero
        if np.isclose(delta_T.in_units('Gyr'), float(0), atol=1e-9):
            logger.debug(f"Delta T is zero for snapshot {snap_num}. Cannot calculate mass loading factor.")
            return float('nan')

    # calculating the mass outflow rate
    mor = np.divide(expelled_mass, delta_T)

    # selecting formed stars within the time interval
    formed_stars = stars[(star_tform <= snap_time) & (star_tform > (snap_time - delta_T.in_units('Gyr')))]
    formed_star_mass = pynbody.units.Unit(f"{formed_stars['mass'].sum().in_units('Msol')} Msol")

    # calculating the star formation rate
    sfr = np.divide(formed_star_mass, delta_T) 

    # calculating the mass loading factor
    if np.isclose(sfr.in_units('Msol Gyr**-1'), float(0), atol=1e-9):
        logger.debug(f"Star formation rate is zero for snapshot {snap_num}. Cannot calculate mass loading factor.")
        return float('nan')
    else:
        # returns dimensionless mass loading factor as raw float
        return np.divide(mor, sfr).ratio(pynbody.units.Unit(1))

def get_halo_masses(key, snap_num):
    """
    Get the mass of a halo at z=0.

    Parameters:
    key (str): The key of desired [simname]_[halo number]
        Ex. cptmarvel_1
    snap_num (int): The snapshot number of the simulation.
    
    Returns:
    float: The mass of the halo in solar masses.
    float: The stellar mass of the halo in solar masses.
    """
    _, halos = load_sim(key.split('_')[0], snap_num)
    halo = halos[int(key.split('_')[1])]
    mass = pynbody.units.Unit(halo['mass'].sum().in_units('Msol')).ratio(pynbody.units.Unit(1))
    solar_mass = pynbody.units.Unit(halo.star['mass'].sum().in_units('Msol')).ratio(pynbody.units.Unit(1))
    return mass, solar_mass
