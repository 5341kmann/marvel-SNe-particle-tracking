import os
import sys
import pickle
import pandas as pd

# Add the project root to sys.path so Python can find file imports
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from sim_enums import MarvelSim, JusticeLeagueSim
from snapnums import *
from base import *

# This file takes the snapshot numbers and main progenitor halo IDs and writes them to a file
# That file can then be read in via the get_stored_simpaths_haloids(sim,z0haloid) function defined in base.py
# Note: The snapnums lists are defined in the snapnums.py file, which is not shown here.


def traceback_to_dictionary(traceback_path):
    """
    Converts a traceback HDF5 file to a dictionary mapping halo IDs to their progenitor IDs.

    Parameters:
    traceback_path (str): Path to the traceback HDF5 file.

    Returns:
    dict: A dictionary where keys are halo IDs and values are lists of progenitor IDs.
    """

    df = pd.read_hdf(traceback_path).apply(list, axis=1).to_dict()
    return {k: [k] + v[: len(df) - 1] for k, v in df.items()}


# File paths for the Marvel simulations
filepaths_cptmarvel = [
    MarvelSim.CAPTAINMARVEL.get_path(snapnum[-4:]) for snapnum in snapnums_cptmarvel
]
filepaths_elektra = [
    MarvelSim.ELEKTRA.get_path(snapnum[-4:]) for snapnum in snapnums_elektra
]
filepaths_rogue = [MarvelSim.ROGUE.get_path(snapnum[-4:]) for snapnum in snapnums_rogue]
filepaths_storm = [MarvelSim.STORM.get_path(snapnum[-4:]) for snapnum in snapnums_storm]

# File paths for theJustice League simulations
filepaths_h148 = [
    JusticeLeagueSim.SANDRA.get_path(snapnum[-4:]) for snapnum in snapnums_h148
]
filepaths_h229 = [JusticeLeagueSim.RUTH.get_path(snapnum[-4:]) for snapnum in snapnums_h229]
filepaths_h242 = [JusticeLeagueSim.SONIA.get_path(snapnum[-4:]) for snapnum in snapnums_h242]
filepaths_h329 = [JusticeLeagueSim.ELENA.get_path(snapnum[-4:]) for snapnum in snapnums_h329]

# File path for traceback data
#
# Traceback paths for the Marvel simulations
traceback_cptmarvel = MarvelSim.CAPTAINMARVEL.get_traceback_path()
traceback_elektra = MarvelSim.ELEKTRA.get_traceback_path()
traceback_rogue = MarvelSim.ROGUE.get_traceback_path()
traceback_storm = MarvelSim.STORM.get_traceback_path()

# Traceback paths for the Justice League simulations
traceback_h148 = JusticeLeagueSim.SANDRA.get_traceback_path()
traceback_h229 = JusticeLeagueSim.RUTH.get_traceback_path()
traceback_h242 = JusticeLeagueSim.SONIA.get_traceback_path()
traceback_h329 = JusticeLeagueSim.ELENA.get_traceback_path()

# converting traceback to dictionary for the Marvel simulations
haloids_cptmarvel = traceback_to_dictionary(traceback_cptmarvel)
# haloids_cptmarvel:
#  1: [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 42, 49, 3],
#  2: [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 4, 4, 4, 4, 4, 4, 9, 2, 2, 2],
#  3: [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 5, 5, 8, 8, 8, 8, 8, 6, 6, 6, 5, 7, 8, 7, 4, 3, 3, 7],
#  5: [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 3, 5, 13, 27],
#  6: [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1],
#  7: [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 6, 6, 6, 6, 6, 5, 5, 5, 5, 5, 7, 8, 8, 8, 8, 7, 8, 13, 12, 18, 93],
#  8: [8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 9, 9, 9, 10, 10, 10, 9, 9, 12, 13, 11, 12, 13, 48, 41, 61, -1, -1, -1],
#  10: [10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 9, 9, 9, 9, 9, 8, 8, 8, 7, 7, 7, 6, 6, 7, 7, 8, 7, 7, 7, 6, 5, 5, 5, 6, 31, 28],
#  11: [11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 10, 11, 11, 11, 11, 11, 9, 11, 13, 11, 11, 10, 9, 10, 7, 4, 8],
#  13: [13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 12, 13, 13, 12, 12, 12, 12, 12, 13, 13, 13, 13, 13, 12, 12, 13, 10, 9, 9, 9, 9, 10, 14, 48, -1, -1],
#  16: [16, 16, 16, 16, 16, 16, 16, 16, 16, 17, 16, 14, 14, 14, 14, 14, 14, 13, 14, 14, 13, 13, 14, 14, 13, 14, 14, 14, 14, 14, 19, 20, 18, 16, 14, 18, 19, 17, 18, 19, 32, -1, 44]
haloids_elektra = traceback_to_dictionary(traceback_elektra)
haloids_rogue = traceback_to_dictionary(traceback_rogue)
haloids_storm = traceback_to_dictionary(traceback_storm)

# converting traceback to dictionary for the Justice League simulations
haloids_h148 = traceback_to_dictionary(traceback_h148)
haloids_h229 = traceback_to_dictionary(traceback_h229)
haloids_h242 = traceback_to_dictionary(traceback_h242)
haloids_h329 = traceback_to_dictionary(traceback_h329)


output = dict(
    filepaths=dict(
        storm=filepaths_storm,
        cptmarvel=filepaths_cptmarvel,
        elektra=filepaths_elektra,
        rogue=filepaths_rogue,
        h148=filepaths_h148,
        h229=filepaths_h229,
        h242=filepaths_h242,
        h329=filepaths_h329,
    ),
    haloids=dict(
        storm=haloids_storm,
        cptmarvel=haloids_cptmarvel,
        elektra=haloids_elektra,
        rogue=haloids_rogue,
        h148=haloids_h148,
        h229=haloids_h229,
        h242=haloids_h242,
        h329=haloids_h329,
    ),
)
with open("Data/simpaths_haloids.pickle", "wb") as f:
    pickle.dump(output, f)
