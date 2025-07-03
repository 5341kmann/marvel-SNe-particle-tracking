import os
import sys
import pickle
import pandas as pd
# Add the project root to sys.path so Python can find 'base.py'
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from base import *

### This file takes the snapshot numbers and main progenitor halo IDs and writes them to a file
### That file can then be read in via the get_stored_simpaths_haloids(sim,z0haloid) function defined in base.py

# Note: Snapshots 000113, 001025, and 000384 are corrupted or empty, so we skip them 
# snapnums_storm = ['004096'] # stom is inomplete, will complete later

snapnums_cptmarvel = ['004096', '003968', '003840', '003712', '003636', '003584', '003456', '003328', '003245', '003200', 
                      '003072', '002944', '002816', '002688', '002624', '002560', '002432', '002304', '002176', '002162',
                      '002048', '001920', '001813', '001792', '001664', '001543', '001536', '001408', '001331', '001280', 
                      '001162', '001152', '000896', '000818', '000768', '000672', '000640', '000512', '000482', '000291']

# excluded '000199', '000146', '000128' due to missing halo IDs in traceback


# snapnums_elektra = [] weird naming scheme, so we skip it for now
# snapnums_rouge = [] will do later

filepaths_cptmarvel = [SIM_FOLDER_PATH+'cptmarvel.cosmo25cmb/cptmarvel.cosmo25cmb.4096g5HbwK1BH/snapshots_200crit_cptmarvel/cptmarvel.cosmo25cmb.4096g5HbwK1BH.'+s for s in snapnums_cptmarvel]
traceback_cptmarvel = 'Data/cptmarvel.trace_back.hdf5'
# filepaths_storm = [SIM_FOLDER_PATH+'storm.cosmo25cmb/storm.cosmo25cmb.4096g5HbwK1BH/snapshots_200crit_storm/storm.cosmo25cmb.4096g5HbwK1BH.'+s for s in snapnums_storm]

haloids_storm = {}

# converting traceback to dictionary - for cptmarvel
df = pd.read_hdf(traceback_cptmarvel).apply(list, axis=1).to_dict()
haloids_cptmarvel = {k: [k] + v for k, v in df.items()}

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


haloids_elektra = {}
haloids_rouge = {}

output = dict(
    filepaths = dict(
        # storm = filepaths_storm,
        cptmarvel = filepaths_cptmarvel,
        # elektra = filepaths_elektra,
        # rogue = filepaths_rouge
    ),
    haloids = dict(
        # storm = haloids_storm,
        cptmarvel = haloids_cptmarvel,
        # elektra = haloids_elektra,
        # rogue = haloids_rouge
    )
)
with open('Data/simpaths_haloids.pickle','wb') as f:
    pickle.dump(output, f)