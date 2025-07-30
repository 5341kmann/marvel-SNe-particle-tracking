from base import *
import h5py
import matplotlib.pyplot as plt
import snapnums
from tqdm.notebook import tqdm
from enum import Enum
from sim_enums import MarvelSim, JusticeLeagueSim


class ParticleType(Enum):
    PERMANENT = "permanent"
    DISK = "disk"
    HALO = "halo"
    EJECTED = "ejected"


def list_items(items):
    """
    List all items in the given list with a comma separator.

    Parameters:
    items (list): The list of items to be listed.

    Returns:
    str: A string of items separated by commas.
    """
    if not items:
        return "No items"
    return ", ".join(str(item) for item in items)


def reload():
    """
    Reload the analysis module to reflect any changes made to the code.
    This is useful during development to see changes without restarting the kernel.
    """
    import importlib
    import analysis

    # from analysis import ParticleType, MarvelSim, JusticeLeagueSim
    # importlib.reload(ParticleType)
    # importlib.reload(MarvelSim)
    # importlib.reload(JusticeLeagueSim)
    importlib.reload(analysis)
    logger.debug("Analysis module reloaded.")


def save_recompute(recompute):
    with open("Data/analysis/recompute.json", "w") as f:
        json.dump(recompute, f, indent=4)


def hdf5_contains_key(hdf5_file, key):
    """
    Check if the given key exists in the HDF5 file.

    Parameters:
    hdf5_file (str): The path to the HDF5 file.
    key (str): The key to check for existence.

    Returns:
    bool: True if the key exists, False otherwise.
    """
    try:
        with h5py.File(hdf5_file, "r") as f:
            return key in f.keys()
    except Exception as e:
        logger.error(f"Error checking key '{key}' in HDF5 file {hdf5_file}: {e}")
        return False


def save_hdf5_to_excel(hdf5_file, row_limit=10000, save_original_path=True):
    """
    Save a pandas DataFrame to an Excel file with multiple sheets.

    parameters:
    hdf5_file (str): The path to the HDF5 file containing the DataFrame.
    """
    if row_limit <= 0:
        logger.debug("Row limit must be a positive integer.")
        return

    if row_limit > 1048576:
        logger.debug(
            "Row limit exceeds Excel's maximum row limit of 1,048,576. Setting to maximum."
        )
        row_limit = 1048576

    if not os.path.exists(hdf5_file):
        logger.debug(f"Error: HDF5 file not found at {hdf5_file}")
        return

    file_keys = []
    with h5py.File(hdf5_file, "r") as f:
        for key in f.keys():
            if isinstance(f[key], h5py.Dataset) or isinstance(f[key], h5py.Group):
                file_keys.append(key)

    if not file_keys:
        logger.debug(f"No keys (datasets) found in the file: {hdf5_file}")
        return

    if save_original_path:
        # Save the Excel file with the same name as the HDF5 file
        excel_file_path = hdf5_file.replace(".hdf5", ".xlsx")
    else:
        # Save the Excel file to the Data directory with original name
        excel_file_path = os.path.join(
            "Data", os.path.basename(hdf5_file).replace(".hdf5", ".xlsx")
        )
    with pd.ExcelWriter(excel_file_path, engine="xlsxwriter") as writer:

        for key in tqdm(file_keys):
            try:
                data = pd.read_hdf(hdf5_file, key=key, stop=row_limit)

                sheet_name = key[:31]

                data.to_excel(writer, sheet_name=sheet_name, index=True)
                logger.debug(f"Saved key '{key}' to sheet '{sheet_name}'")
            except Exception as e:
                logger.debug(f"Error saving key '{key}' to Excel: {e}")


def calculate_massive_haloids(justice_league=False):
    """
    Calculate the IDs of massive halos in the simulations based on stellar mass.
    Massive halos are defined as those with stellar masses between 1e5 and 1e9 solar masses.
    The results are saved to a JSON file.
    If the file already exists, it will load the existing IDs instead of recalculating.

    Note:
    This function runs quite quickly, but at the cost of memory usage.
    Parameters:
    justice_league (bool): If True, use Justice League simulations; otherwise, use Marvel simulations.

    """

    if os.path.exists("Data/analysis/massive_haloids.json"):
        logger.debug("Massive halo IDs file exists. loading in existing IDs.")
        with open("Data/analysis/massive_haloids.json", "r") as f:
            massive_halo_ids = json.load(f)
    else:
        logger.debug(
            "Massive halo IDs file does not exist. Initializing empty dictionary."
        )
        massive_halo_ids = {}

    for sim in tqdm(
        JusticeLeagueSim if justice_league else MarvelSim,
        desc=f"Calculating massive halo IDs",
    ):

        # checks if the simulation already has massive halo IDs calculated
        if sim.value in massive_halo_ids:
            logger.debug(
                f"Massive halo IDs for {sim.value} already exist. Skipping calculation."
            )
            continue

        s, h = load_sim(sim, 4096)

        halos = h[1:]  # Exclude the first halo which is the whole simulation

        # loading in all the star particle IDs and their masses
        all_sim_star_iords = s.star["iord"]
        all_sim_star_masses = s.star["mass"].in_units("Msol")

        # Create a Series to map star particle IDs to their masses
        star_mass_by_iord = pd.Series(all_sim_star_masses, index=all_sim_star_iords)

        # prefilling stellar masses array to avoid appending repeatedly
        calculated_stellar_masses_Msol = np.empty(len(halos))
        calculated_stellar_masses_Msol.fill(0.0)

        # Iterate over each halo and sum the stellar masses
        for i, halo in enumerate(
            tqdm(halos, desc=f"Summing stellar masses for {sim.value}")
        ):

            # loading in star particle IDs and iords for the current halo
            halo_star_iords = np.array(halo.star["iord"])

            # check if halo has star particles
            if len(halo_star_iords) > 0:

                # filtering for only star particles that are in the global star particle IDs (i.e. have masses)
                mask_in_global = np.isin(halo_star_iords, star_mass_by_iord.index)

                if np.any(mask_in_global):
                    # Sum the masses of the star particles found in this halo
                    calculated_stellar_masses_Msol[i] = star_mass_by_iord.loc[
                        halo_star_iords[mask_in_global]
                    ].sum()

        # stellar mass cuts
        min_mass_cut = 1e5
        max_mass_cut = 1e9

        mass_filter_bool_array = (calculated_stellar_masses_Msol >= min_mass_cut) & (
            calculated_stellar_masses_Msol <= max_mass_cut
        )

        massive_halos_filtered = [
            halo
            for i, halo in tqdm(
                enumerate(halos), desc=f"selecting massive halos for {sim.value}"
            )
            if mass_filter_bool_array[i]
        ]
        # storing the halo IDs of the massive halos in a dictionary
        massive_halo_ids[sim.value] = [
            int(h.properties["halo_number"])
            for h in tqdm(
                massive_halos_filtered, desc=f"collecting halo numbers for {sim.value}"
            )
        ]

    with open("Data/analysis/massive_haloids.json", "w") as f:
        json.dump(massive_halo_ids, f, indent=4)


def get_massive_haloids(sim_enum_member):
    """
    Returns a list of halo IDs for the Marvel simulations with the massive halo mass cut
    """
    logger = set_logger(sim_enum_member, "get_massive_haloids")

    if not os.path.exists("Data/analysis/massive_haloids.json"):
        logger.debug(
            "Massive halo IDs file not found. calculating massive halos first."
        )
        calculate_massive_haloids()

    with open("Data/analysis/massive_haloids.json", "r") as f:
        massive_halos = json.load(f)

    if sim_enum_member.value not in massive_halos:
        logger.debug(
            f"No massive halo ID for simulation {sim_enum_member.value} Running calculate_massive_haloids()"
        )
        calculate_massive_haloids()
        # reloading the massive halos after calculation
        with open("Data/analysis/massive_haloids.json", "r") as f:
            massive_halos = json.load(f)

    if sim_enum_member == MarvelSim.CAPTAINMARVEL:
        return massive_halos[sim_enum_member.value]
    elif sim_enum_member == MarvelSim.ELEKTRA:
        return massive_halos[sim_enum_member.value]
    elif sim_enum_member == MarvelSim.ROGUE:
        return massive_halos[sim_enum_member.value]
    elif sim_enum_member == MarvelSim.STORM:
        return massive_halos[sim_enum_member.value]
    elif sim_enum_member == JusticeLeagueSim.SANDRA:
        return massive_halos[sim_enum_member.value]
    elif sim_enum_member == JusticeLeagueSim.RUTH:
        return massive_halos[sim_enum_member.value]
    elif sim_enum_member == JusticeLeagueSim.SONIA:
        return massive_halos[sim_enum_member.value]
    elif sim_enum_member == JusticeLeagueSim.ELENA:
        return massive_halos[sim_enum_member.value]
    else:
        logger.debug(
            f"Unknown simulation: {sim_enum_member.name}. Cannot get massive halo IDs."
        )
        return []


def get_columns_of_hdf5(hdf5_path):
    """
    Gets the column names of the DataFrame stored under the first
    top-level key in an HDF5 file, without loading all the data.

    Parameters:
    file_path (str): The path to the HDF5 file.

    Returns:
    list: A list of column names, or None if no keys are found or an error occurs.
    """
    logger.debug(
        f"Attempting to get columns from the first key in HDF5 file: {hdf5_path}"
    )
    try:
        with h5py.File(hdf5_path, "r") as f:
            all_keys = list(f.keys())

        # empty file
        if not all_keys:
            logger.debug(f"No keys (datasets) found in the HDF5 file: {hdf5_path}")
            return None

        # get the first key
        first_key = all_keys[0]
        with pd.HDFStore(hdf5_path, mode="r") as store:
            if first_key in store:
                node = store.get_node(first_key)
                if hasattr(node, "description") and hasattr(
                    node.description, "_v_names"
                ):
                    return list(node.description._v_names)
                else:
                    # Fallback for other formats or if structure isn't direct
                    logger.warning(
                        f"Could not get columns from PyTables description. Falling back to store.get_storer()."
                    )
                    # for HDFStore tables not generic h5py datasets
                    return list(store.get_storer(first_key).non_index_axes[0][1])
            else:
                logger.error(
                    f"First key '{first_key}' was listed but not found by HDFStore."
                )
                return None
    # No HDF5 file found
    except Exception as e:
        logger.error(f"Error getting columns from {hdf5_path}: {e}")
        return None


def load_ejected_particles(sim, halo_num):
    """
    Load the ejected particles data from an HDF5 file.

    Parameters:
    key (str): The key of desired [simname]_[halo number]
        Ex. cptmarvel_1

    Returns:
    pd.DataFrame: The ejected particles data.
    """
    key = f"{sim.value}_{halo_num}"
    # if sim in JusticeLeagueSim:
    #     path = JL_SNE_DATA_PATH + 'ejected_particles.hdf5'
    # else:
    path = "Data/SNe/ejected_particles.hdf5"

    try:
        return pd.read_hdf(path, key=key)
    except KeyError:
        logger.debug(
            f"No ejected particles data found for key {key}. Verify calc_ejected_expelled() logs for details."
        )
        try:
            logger.warning(
                f"Presumably zero ejected particles found. Returning empty DataFrame with columns from ejected_particles.hdf5."
            )
            return pd.DataFrame(columns=get_columns_of_hdf5(path))
        except Exception as e:
            logger.error(f"No found columns in ejected particles HDF5 file: {e}")
            return None
    except FileNotFoundError:
        logger.error(f"Ejected particles HDF5 file not found at {path}.")
        return None


def load_expelled_particles(sim, halo_num):
    """
    Load the ejected particles data from an HDF5 file.

    Parameters:
    sim (MarvelSim or JusticeLeagueSim): The simulation enum.
    halo_num (int): The halo number to load data for.

    Returns:
    pd.DataFrame: The ejected particles data.
    """
    key = f"{sim.value}_{halo_num}"
    # if sim in JusticeLeagueSim:
    #     path = JL_SNE_DATA_PATH + 'expelled_particles.hdf5'
    # else:
    path = "Data/SNe/expelled_particles.hdf5"

    try:
        return pd.read_hdf(path, key=key)
    except KeyError:
        logger.debug(
            f"No expelled particles data found for key {key}. Verify calc_ejected_expelled() logs for details."
        )
        try:
            logger.warning(
                f"Presumably zero expelled particles found. Returning empty DataFrame with columns from expelled_particles.hdf5."
            )
            return pd.DataFrame(columns=get_columns_of_hdf5(path))
        except Exception as e:
            logger.error(f"No found columns in expelled particles HDF5 file: {e}")
            return None
    except FileNotFoundError:
        logger.error(f"Expelled particles HDF5 file not found at {path}.")
        return None


def load_tracked_particles(sim, halo_num):
    """
    Load the tracked particles data from an HDF5 file.

    Parameters:
    key (str): The key of desired [simname]_[halo number]
        Ex. cptmarvel_1

    Returns:
    pd.DataFrame: The tracked particles data.
    """
    if sim in JusticeLeagueSim:
        path = JL_TRACKED_PARTICLES_PATH
    else:
        path = "Data/tracked_particles.hdf5"
    return pd.read_hdf(path, key=f"{sim.value}_{halo_num}")


def load_permanently_expelled_particles(sim, halo_num):
    """
    Load the permantly expelled particles data from an HDF5 file.

    Parameters:
    key (str): The key of desired [simname]_[halo number]
        Ex. cptmarvel_1

    Returns:
    pd.DataFrame: The permantly expelled particles data.
    """
    return pd.read_hdf(
        "Data/SNe/permanent_expelled_particles.hdf5", key=f"{sim.value}_{halo_num}"
    )


def load_sim(sim, snapnum):
    """
    Load the simulation the data/Sims folder

    Parameters:
    snapnum (str): The snapshot number of the simulation.

    Returns:
    pynbody.snapshot.Snapshot: The loaded simulation snapshot.
    pynbody.halo.Halo: The halos in the simulation.
    """
    if sim in JusticeLeagueSim:
        path = sim.get_path(snapnum)
        if not os.path.exists(path):
            logger.error(
                f"Simulation path {path} does not exist. Check simulation name and snapshot number."
            )
            raise FileNotFoundError(f"Simulation path {path} does not exist.")
    else:
        path = f"{SIM_FOLDER_PATH}/{sim.value}.cosmo25cmb/{sim.value}.cosmo25cmb.4096g5HbwK1BH/snapshots_200crit_{sim.value}/{sim.value}.cosmo25cmb.4096g5HbwK1BH.00{snapnum}"

    sim = pynbody.load(path)
    sim.physical_units()
    halos = sim.halos(halo_numbers="v1")
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

    d = pickle.load(open("Data/simpaths_haloids.pickle", "rb"))
    ids = d["haloids"][sim][z0halo]
    paths = d["filepaths"][sim][z0halo]
    return ids, paths


def halo_mlf_by_tracked(sim, z0halo, snap_num, identify=None):
    """
    Get the mass loading factor of a halo at z=0.

    Parameters:
    key (str): The key of desired [simname]_[halo number]
        Ex. cptmarvel_1
    snap_num (int): The snapshot number of the simulation.
    delta_T (float): The time interval over which to calculate the star formation rate, in Gyr.

    Returns:
    float: The mass loading factor of the halo.
    str: The snapshot number of the the calculated mass loading factor.
    """
    halo_num = z0halo
    halo_logger = set_logger(sim, halo_num, "halo_mlf_by_tracked")

    halo_logger.debug("Starting mass loading factor calculation for halo.")

    # check identify
    if identify is not None:
        # single particle type
        if isinstance(identify, ParticleType):
            identify_list = [identify]
        # list of particle types
        elif isinstance(identify, list):
            if not all(isinstance(item, ParticleType) for item in identify):
                halo_logger.debug(
                    f"Invalid identify list: {identify}. All items must be a list of ParticleType enum members."
                )
            identify_list = identify
        # invalid type
        else:
            halo_logger.debug(
                f"Invalid identify type: {type(identify)}. Must be ParticleType enum or a list of ParticleType enums."
            )
            return float("nan")
    else:
        # default to all particle types
        identify_list = [ParticleType.DISK, ParticleType.HALO]

    # loading in snashot numbers
    try:
        snap_numbers = np.array(
            [num[-4:] for num in getattr(snapnums, f"snapnums_{sim.value}")]
        )
    except AttributeError:
        halo_logger.debug(
            f"Could not find snap numbers for simulation {sim.value}. Ensure generate_haloid_file.py has been correctly \
                    populated with the snapshot numbers."
        )
        return float("nan")

    # Checking for most recent snapshot containing star forming particles
    has_formed_stars = False

    if sim in JusticeLeagueSim:
        # loading in quenching tables for begining snap shot number
        with open("Data/analysis/justice_league_quench_table.pickle", "rb") as f:
            quenching_tables = pickle.load(f)

        # looking up the snap number for sim and z0 halo
        # if sim.value is not in quenching_tables
        if sim.value not in quenching_tables["sim"].unique():
            halo_logger.debug(
                f"Simulation {sim.value} not found in quenching tables. please populate Data/analysis/justiceleague_quenching_tables.pickle with {sim.value} data."
            )
            halo_logger.debug(
                f"Starting search for star forming particles at snap {snap_num} for {sim.value}_{z0halo}."
            )
        # if sim.value is in quenching tables
        else:
            # getting the quenching data for the current sim and halo
            quench_data = quenching_tables.loc[
                (quenching_tables["sim"] == sim.value)
                & (quenching_tables["haloid"] == z0halo),
            ][["haloid_snap", "quenched", "quench_snap"]]

            # if quench_data is empty or not quenched
            if quench_data.empty or quench_data["quenched"].iloc[0] == False:
                halo_logger.debug(
                    f"No quenching data found for {sim.value}_{z0halo}. Starting search for star forming particles at snap {snap_num}."
                )
            else:
                # getting the snap number of the last snapshot before quenching
                snap_num = quench_data["quench_snap"].iloc[0]
                halo_num = int(quench_data["haloid_snap"].iloc[0])
                halo_logger.debug(
                    f"Using {sim.value}_{snap_num} at {snap_num} for {sim.value}_{z0halo} "
                )
    elif sim in MarvelSim:
        if sim == MarvelSim.CAPTAINMARVEL:
            if z0halo == 13:
                # Special case for Captain Marvel halo 13, 
                snap_num = '3200'
                halo_logger.debug(
                    f"Skipping to snap {snap_num} for {sim.value}_{z0halo} due to known lack of star formation."
                )
            if z0halo == 10:
                # Special case for Captain Marvel halo 10,
                snap_num = '2304'
                halo_logger.debug(
                    f"Skipping to snap {snap_num} for {sim.value}_{z0halo} due to known lack of star formation."
                )
            if z0halo == 11:
                # Special case for Captain Marvel halo 11,
                snap_num = '0896'
                halo_logger.debug(
                    f"Skipping to snap {snap_num} for {sim.value}_{z0halo} due to known lack of star formation."
                )
        if sim in [MarvelSim.ELEKTRA, MarvelSim.ROGUE, MarvelSim.STORM]:
            # special case for Elektra, Rogue, and Storm halos
            snap_num = '4096'

    while not has_formed_stars:
        print(
            f"Checking for formed stars in halo {sim.value}_{halo_num} at snap {snap_num}."
        )

        # loading in coressponding halos
        s, h = load_sim(sim, snap_num)

        # get index of current timestep
        current_snapnum_index = np.where(snap_numbers == snap_num)[0][0]

        # get halo_num from trace back at current snapshot

        halo = h[halo_num]

        # extracting stellar information
        stars = halo.star
        star_tform = stars["tform"].in_units("Gyr")

        # loading in current snapshot time, This takes an oddly long time
        snap_time = s.properties["time"].in_units("Gyr")

        # calcuating time from the last snapshot: delta_T
        if current_snapnum_index == len(snap_numbers) - 1:
            logger.debug(
                f"Current snapshot {snap_num} is the last snapshot in the available snapshots. Cannot calculate mass loading factor."
            )
            return float("nan")
        else:
            prev_snapnum = snap_numbers[current_snapnum_index + 1]
            if sim in [MarvelSim.ELEKTRA, MarvelSim.ROGUE, MarvelSim.STORM]:
                dt = 4096/snap_time
                delta_T = pynbody.units.Unit(f"{prev_snapnum*dt} Gyr")
            delta_T = snap_time - load_sim(sim, prev_snapnum)[0].properties[
                "time"
            ].in_units("Gyr")
            delta_T = pynbody.units.Unit(f"{delta_T} Gyr")
            # checking if delta_T is zero
            if np.isclose(delta_T.in_units("Gyr"), float(0), atol=1e-9):
                halo_logger.debug(
                    f"Delta T is zero between {snap_num} and {prev_snapnum}. Cannot calculate mass loading factor."
                )
                return float("nan")

        # selecting formed stars within the time interval
        formed_stars = stars[
            (star_tform <= snap_time)
            & (star_tform > (snap_time - delta_T.in_units("Gyr")))
        ]

        if len(formed_stars) == 0:
            if sim in [MarvelSim.ELEKTRA, MarvelSim.ROGUE, MarvelSim.STORM]:
                halo_logger.debug(
                    f'Cannot find formed stars for {sim.value}_{halo_num} at snap {snap_num}. Stopping due to no accessible previous snapshots.'
                )
                print('Cannot find formed stars. Stopping due to no accessible previous snapshots.')
                return float("nan"), snap_num
            halo_logger.debug(
                f"No formed stars found for {sim.value}_{halo_num} at snap {snap_num}. Trying previous snapshot."
            )
            snap_num = prev_snapnum
            print(
                f"Star formation zero. Trying previous snapshot {snap_num} for formed stars."
            )
            continue
        else:
            halo_logger.debug(
                f"Found {len(formed_stars)} formed stars for {sim.value}_{halo_num} at snap {snap_num}."
            )
            has_formed_stars = True

    # loading in tracked expelled particles
    if ParticleType.PERMANENT in identify_list:
        expelled_particles = load_permanently_expelled_particles(sim, z0halo)
        print(
            "Using permanently expelled particles for mass loading factor calculation."
        )
    else:
        expelled_particles = load_expelled_particles(sim, z0halo)
        print("Using expelled particles for mass loading factor calculation.")

    # time selected expelled particles
    expelled_in_time = expelled_particles[
        np.isclose(expelled_particles["time"], snap_time, atol=1e-2)
    ]

    # Determining outflow mass by selected idetification of outrow particles
    totoal_outflow_mass = 0.0

    if ParticleType.DISK in identify_list:
        # time selected ejected particles
        expelled_mass = pynbody.units.Unit(
            f"{expelled_in_time[expelled_in_time['prev_gal_disk'] == True]['mass'].sum()} Msol"
        )

        # only add if non-zero mass to avoid divide by zero errors of pynbody's CompositeUnit.expand() handling
        if not np.isclose(expelled_mass.in_units("Msol"), float(0), atol=1e-9):
            totoal_outflow_mass += expelled_mass.in_units("Msol")
        else:
            halo_logger.debug(
                f"No disk-expelled particles found for {sim.value}_{halo_num} at snap {snap_num}."
            )

    if ParticleType.HALO in identify_list:
        # time selected ejected particles
        expelled_mass = pynbody.units.Unit(
            f"{expelled_in_time[expelled_in_time['prev_gal_halo'] == True]['mass'].sum()} Msol"
        )

        if not np.isclose(expelled_mass.in_units("Msol"), float(0), atol=1e-9):
            totoal_outflow_mass += expelled_mass.in_units("Msol")
        else:
            halo_logger.debug(
                f"No halo-expelled particles found for {sim.value}_{halo_num} at snap {snap_num}."
            )

    if ParticleType.EJECTED in identify_list:
        ejected_particles = load_ejected_particles(sim, halo_num)
        selected_ejected = ejected_particles[
            np.isclose(ejected_particles["time"], snap_time, atol=1e-2)
        ]
        ejected_mass = pynbody.units.Unit(f"{selected_ejected['mass'].sum()} Msol")

        if not np.isclose(ejected_mass.in_units("Msol"), float(0), atol=1e-9):
            totoal_outflow_mass += ejected_mass.in_units("Msol")
        else:
            halo_logger.debug(
                f"No ejected particles found for {sim.value}_{halo_num} at snap {snap_num}."
            )

    totoal_outflow_mass = pynbody.units.Unit(f"{totoal_outflow_mass} Msol")

    # calculating the mass outflow rate
    mor = np.divide(totoal_outflow_mass, delta_T)

    # getting mass of formed stars
    formed_star_mass = pynbody.units.Unit(
        f"{formed_stars['mass'].sum().in_units('Msol')} Msol"
    )

    # calculating the star formation rate
    sfr = np.divide(formed_star_mass, delta_T)

    # calculating the mass loading factor
    if np.isclose(sfr.in_units("Msol Gyr**-1"), float(0), atol=1e-5):
        halo_logger.debug(
            f"Star formation rate ({sfr}) is zero for {sim.value}_{halo_num} at snapshot{snap_num}. Cannot calculate mass loading factor."
        )
        return float("nan")
    else:
        # returns dimensionless mass loading factor as raw float
        return np.divide(mor, sfr).ratio(pynbody.units.Unit(1)), snap_num


def halo_mlf_by_flux(key, snap_num, r_center, r_width):
    """
    Get the mass flux of a halo at z=0 within a specific radius.

    Note:
    depends on the existence of the snapnums_[simname] (Ex. snapnums_captmarvel) variable being populated in descending order in
    generate_haloid_file.py, which is used to calculate the time interval over which to calculate the star formation rate.

    Parameters:
    key (str): The key of desired [simname]_[halo number]
        Ex. cptmarvel_1

    snap_num (int): The snapshot number of the simulation.
    r_center (float): center of anulus through which to calculate the mass loading factor, a multiple of rvir.
    r_width (float): width of the anulus within which to calculate the mass loading factor, a multiple of rvir.

    returns:
    Ty: A dimensionless mass loading using the mass flux of the halo within the specified anulus.
    #"""

    sim, halos = load_sim(key.split("_")[0], snap_num)
    halo = halos[int(key.split("_")[1])]
    halo.physical_units()
    pynbody.analysis.halo.center(halo)
    snap_time = sim.properties["time"].in_units("Gyr")
    snap_numbers = np.array(
        [int(num) for num in globals()[f'snapnums_{key.split("_")[0]}']]
    )

    # extracting virial radius
    rvir = halo.properties["Rvir"]
    # check for non-positive virial radius
    if rvir <= 0:
        logger.debug(
            f"Virial radius is non-positive ({rvir} kpc) for halo {key} at snap {snap_num}. Cannot define annulus."
        )
        return float("nan")
    r_low = (r_center - (r_width / 2)) * rvir
    # check for negative lower radius
    if r_low < 0:
        logger.debug(
            f"Lower radius {r_low} kpc is negative when using a annulus width of {r_width} for halo {key} at snap {snap_num}. Cannot define annulus."
        )
        return float("nan")
    r_high = (r_center + (r_width / 2)) * rvir
    # extracting stellar information
    stars = halo.star
    star_tform = stars["tform"].in_units("Gyr")

    # calculating the time from the last snapshot: delta_T
    current_snapnum_index = np.where(snap_numbers == snap_num)[0][0]
    if current_snapnum_index == len(snap_numbers) - 1:
        logger.debug(
            f"Current snapshot {snap_num} is the last snapshot in the simulation. Cannot calculate mass loading factor."
        )
        return float("nan")
    else:
        prev_snapnum = snap_numbers[current_snapnum_index + 1]
        delta_T = snap_time - load_sim(key.split("_")[0], prev_snapnum)[0].properties[
            "time"
        ].in_units("Gyr")
        delta_T = pynbody.units.Unit(f"{delta_T} Gyr")
        # checking if delta_T is zero
        if np.isclose(delta_T.in_units("Gyr"), float(0), atol=1e-9):
            logger.debug(
                f"Delta T is zero for snapshot {snap_num}. Cannot calculate mass loading factor."
            )
            return float("nan")
    # calculating formed stars within the time interval
    formed_stars = stars[
        (star_tform <= snap_time) & (star_tform > (snap_time - delta_T.in_units("Gyr")))
    ]
    formed_star_mass = pynbody.units.Unit(
        f"{formed_stars['mass'].sum().in_units('Msol')} Msol"
    )
    # checking if formed star mass is zero
    if np.isclose(formed_star_mass.in_units("Msol"), float(0), atol=1e-9):
        logger.debug(
            f"Star formation mass is zero (or near zero) for {key}. Cannot calculate star formation rate."
        )
        return float("nan")
    # calculating the star formation rate
    sfr = np.divide(formed_star_mass, delta_T)

    # calculate the mass outflow
    gas = halo.gas
    # filtering gas particles in the annulus and with positive radial velocity
    gas = gas[(gas["r"] >= r_low) & (gas["r"] < r_high) & (gas["vr"] > 0)]

    outlfow_Msol_per_sec = pynbody.units.Unit(
        f"{np.divide(np.multiply(gas['mass'].in_units('Msol'),gas['vr'].in_units('kpc s**-1')),r_width).sum()} Msol s**-1"
    )  # Msol s**-1

    outflow_rate = pynbody.units.Unit(
        f"{outlfow_Msol_per_sec.in_units('Msol Gyr**-1')} Msol Gyr**-1"
    )

    mlf = np.divide(outflow_rate, sfr)
    return mlf.ratio(pynbody.units.Unit(1))


def compute_sim_mlfs(
    func,
    recompute_dict,
    sim,
    z0halos,
    snap_num,
):
    """
    Compute the mass loading factors for all halos in a simulation at z=0.

    Parameters:
    func (function): The function to compute the mass loading factor.
    sim_enum_name (str): The name of the simulation enum (e.g., 'cptmarvel').
    z0halos (list): List of halo IDs at z=0.
    snap_num (str): The snapshot number of the simulation.

    Returns:
    pandas.DataFram: A DataFrame containing the mass loading factors and their corresponding snapshot numbers.
    """
    mlfs = []
    mlf_snaps = []

    for z0halo in tqdm(z0halos, desc=f"Computing MLFs for {sim.value}"):
        key = f"{sim.value}_{z0halo}"

        # checking if the halo needs to be recomputed
        if recompute_dict[key]["expelled_mlf"]:
            print(f"Recomputing halo {z0halo} MLF")

            # recomputing the MLF for the halo
            mlf, mlf_snap = func(sim, z0halo, snap_num)
            print(f"Finished halo {z0halo}\nmlf: {mlf} at snapshot {mlf_snap}")

            # updating recompute settings
            recompute_dict[key]["expelled_mlf"] = False

            # appending the MLF and snapnumbers
            mlfs.append(mlf)
            mlf_snaps.append(mlf_snap)

        # if the halo does not need to be recomputed based on the recompute_dict
        else:
            print(f"Attempting to load {key} from disk...")

            # Loading in MLF data old format
            try:
                df = pd.read_hdf(f"Data/analysis/{key}.hdf5", key="expelled_mlf")

                # appending the MLF to the list
                if len(df) == 1:
                    mlfs.append(df["mlf"].iloc[0])
                    mlf_snaps.append(df["snap_num"].iloc[0])
                else:
                    # If multiple entries, find the one for the current halo
                    i = z0halos.index(z0halo)
                    mlfs.append(df["mlf"].iloc[i])
                    mlf_snaps.append(df["snap_num"].iloc[i])

                print(
                    f"Loaded {key} from disk. MLF: {mlfs[-1]} at snapshot {mlf_snaps[-1]}"
                )

            except FileNotFoundError:
                # loading in new file format
                try:
                    df = pd.read_hdf(
                        f"Data/analysis/{sim.value}.hdf5", key="expelled_mlf"
                    )
                    # appending the MLF to the list
                    mlfs.append(df["mlf"].loc[z0halo])
                    mlf_snaps.append(df["snap_num"].loc[z0halo])
                    recompute_dict[key]["expelled_mlf"] = True

                except KeyError:
                    # if the key is not found, it means the halo has no MLF data
                    logger.error(
                        f"No MLF data found for {key}. Please recompute the MLF."
                    )
                    mlfs.append(float("nan"))
                    mlf_snaps.append(float("nan"))
                    recompute_dict[key]["expelled_mlf"] = True

    # writie data to dictionary
    df = pd.DataFrame({
        "mlf": mlfs,
        "snap_num": mlf_snaps,
    },
        index=z0halos,
    )
    df.index.name = "halo_id"
    # saving the MLFs to an HDF5 file
    save_with_lock(df, f"Data/analysis/{sim.value}.hdf5", key="expelled_mlf")

    # saving recompute settings
    save_recompute(recompute_dict)

    return df


def compute_sim_masses(sim, z0halos, recompute_dict, snap_numbers):
    """
    Compute the masses of all halos in a simulation at z=0.

    Parameters:
    sim (MarvelSim or JusticeLeagueSim): The simulation enum.
    z0halos (list): List of halo IDs at z=0.
    snap_num (int): The snapshot number of the simulation.
    recompute_dict (dict): Dictionary containing recompute settings for halos.
    snap_numbers (list): List of snapshot numbers for the simulation corresponding to the halos.

    Returns:
    pandas.DataFrame: A DataFrame containing the halo IDs and their corresponding masses.
    """
    masses = []
    stellar_masses = []

    for z0halo, snap_num in tqdm(zip(z0halos, snap_numbers), desc=f"Computing masses for {sim.value}", total=len(z0halos)):
        key = f"{sim.value}_{z0halo}"
        if recompute_dict[key]["halo_mass"] or recompute_dict[key]["stellar_mass"]:
            print(f"Recomputing halo {z0halo} masses")
            mass, stellar_mass = get_halo_masses(sim, z0halo, snap_num)
            masses.append(mass)
            stellar_masses.append(stellar_mass)
            recompute_dict[f"{sim.value}_{z0halo}"]["halo_mass"] = False
            recompute_dict[f"{sim.value}_{z0halo}"]["stellar_mass"] = False
            save_recompute(recompute_dict)
        else:
            print(f"Attempting to load {key} from disk...")
            try:
                df = pd.read_hdf(f"Data/analysis/{sim.value}.hdf5", key="mass_properties")
                # appending the MLF to the list
                masses.append(df["halo_mass"].loc[z0halo])
                stellar_masses.append(df["stellar_mass"].loc[z0halo])
            except:
                print(f"No mass data found for {key}")
                logger.error(
                    f"No mass data found for {key}. Please recompute the masses."
                )
                recompute_dict[key]["halo_mass"] = True
                recompute_dict[key]["stellar_mass"] = True
                save_recompute(recompute_dict)
                masses.append(float("nan"))
                stellar_masses.append(float("nan"))
                continue
            #

    # creating a DataFrame with the results
    df = pd.DataFrame(
        {
            "halo_mass": masses,
            "stellar_mass": stellar_masses,
        },
        index=z0halos,
    )
    df.index.name = "halo_id"

    # saving the masses to an HDF5 file
    save_with_lock(df, f"Data/analysis/{sim.value}.hdf5", key="mass_properties")
    return df


def get_halo_masses(sim, z0halo, snap_num):
    """
    Get the mass of a halo at z=0.

    Parameters:
    sim (MarvelSim or JusticeLeagueSim): The simulation enum.
    z0halo (int): The halo number at z=0.
    snap_num (int): The snapshot number of the simulation.

    Returns:
    float: The mass of the halo in solar masses.
    float: The stellar mass of the halo in solar masses.
    """
    _, halos = load_sim(sim, snap_num)
    halo = halos[z0halo]
    mass = pynbody.units.Unit(halo["mass"].sum().in_units("Msol")).ratio(
        pynbody.units.Unit(1)
    )
    solar_mass = pynbody.units.Unit(halo.star["mass"].sum().in_units("Msol")).ratio(
        pynbody.units.Unit(1)
    )
    return mass, solar_mass


def plot_frequency(hdf5_file_path, key, column_name):

    try:
        # Efficiently load only the specified column
        # Pandas can optimize this to read just that column's data from disk.
        series = pd.read_hdf(hdf5_file_path, key=key, columns=[column_name])

        # Access the Series for the specific column
        # (read_hdf with columns returns a DataFrame, so select the column)
        target_column_data = series[column_name]

        if pd.api.types.is_float_dtype(target_column_data.dtype):
            target_column_data = target_column_data.round(7)

        # Perform frequency count on unique values
        # .value_counts() is highly optimized for this.
        # It returns a Series with unique values as index and their counts as values.
        frequency_counts = target_column_data.value_counts(dropna=False)

        # Sort the counts by index (the unique values) for cleaner plotting if they are numeric/categorical
        frequency_counts = frequency_counts.sort_index()
        return frequency_counts
        # Plotting
        # plt.figure(figsize=(10, 6))

        # # If the unique values are categorical or a few numbers
        # if frequency_counts.index.dtype == 'object' or len(frequency_counts) < 20: # Heuristic for categories
        #     # Use bar plot for discrete categories
        #     frequency_counts.plot(kind='bar')
        #     plt.xticks(rotation=45, ha='right') # Rotate labels if they overlap
        # else:
        #     # Use a histogram or scatter plot for continuous/many unique values
        #     # For frequency counts of unique values, a bar plot is still appropriate if index is discrete.
        #     # If the data is truly continuous, consider a histogram directly on `target_column_data`
        #     # or binning `frequency_counts` for a cleaner bar plot.
        #     frequency_counts.plot(kind='bar') # Defaulting to bar for value_counts output

        # plt.xlabel(column_name)
        # plt.ylabel('Frequency')
        # plt.title(f'Frequency of Unique Values in {column_name} (Key: {key})')
        # plt.grid(axis='y', linestyle='--', alpha=0.7)
        # plt.tight_layout()
        # plt.show()

    except KeyError:
        print(f"Error: Key '{key}' or column '{column_name}' not found in HDF5 file.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")


def plot_Dataframe(
    data,
    xlabel,
    ylabel,
    title,
    save=False,
):
    cols_to_plot = [col for col in data.columns]

    plt.figure(figsize=(24, 12))  # Adjust figure size as needed

    # Plot each selected column as a separate line
    # Pandas will automatically:
    # - Use the DataFrame's index for the x-axis
    # - Plot each column in `cols_to_plot` as a line
    # - Assign unique colors
    # - Create a legend using column names
    data[cols_to_plot].plot(
        ax=plt.gca(),  # Get current axes to ensure title/labels apply correctly
        kind="line",
        marker="o",  # Circle markers for data points
        # linestyle='-', # Solid lines
        alpha=0.8,  # Transparency
    )

    plt.xlabel(xlabel)  # Label for the categories on x-axis
    plt.ylabel(ylabel)  # Label for the frequencies on y-axis
    plt.title(title)
    plt.grid(True, linestyle="--", alpha=0.6)  # Add a grid
    plt.legend(title="Series/Key")  # Title the legend

    # Rotate x-axis labels if they are categorical or long to prevent overlap
    plt.xticks(rotation=45, ha="right")

    plt.tight_layout()  # Adjust layout to prevent labels/legend from overlapping
    plt.savefig(f"Data/plots/{title.replace(' ', '_')}.png") if save else None
    plt.show()
