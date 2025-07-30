import os
import sys
import json
from tqdm.notebook import tqdm

basepath = os.path.dirname(os.path.dirname(os.getcwd()))
# print(os.getcwd())

sys.path.append(os.getcwd())

import analysis as an
for simulation in [an.MarvelSim.ELEKTRA, an.MarvelSim.ROGUE, an.MarvelSim.STORM]:
    
    z0halos = an.get_massive_haloids(simulation)
    sim_name = simulation.value
    print(f'Analyzing {sim_name} with {len(z0halos)} z=0 halos: {an.list_items(z0halos)}')
    # attempting to load in recompute settings
    if os.path.exists("Data/analysis/recompute.json"):
        with open("Data/analysis/recompute.json", "r") as f:
            recompute = json.load(f)
            print(
                f"Loaded recompute settings for {sim_name}: {an.list_items(recompute.keys())}"
            )
    else:
        recompute = {}
        print("recompute.json not found, using default settings.")
        for z0halo in tqdm(z0halos, desc="Loading z=0 halos"):
            key = f"{sim_name}_{z0halo}"
            recompute[key] = {
                "halo_mass": True,
                "stellar_mass": True,
                "expelled_mlf": True,
                "expelled_flux_mlf": {},
                "disk_mlf": True,
            }

            print(f"Recompute settings for {key}: {recompute[key]}")
    for z0halo in z0halos:
        key = f'{sim_name}_{z0halo}'
        if key not in recompute or not os.path.exists(f'Data/analysis/{key}.hdf5'):
            recompute[key] = {'halo_mass' : True,
                            'stellar_mass' : True,
                            'expelled_mlf' : True,
                            'expelled_flux_mlf' : {}, 
                            'disk_mlf' : True}
            print(f"Added recompute settings for {key}: {recompute[key]}")

    mlfs = an.compute_sim_mlfs(an.halo_mlf_by_tracked, 
                        recompute_dict=recompute, 
                        sim=simulation,
                        z0halos=z0halos, 
                        snap_num='4096')

    list_mlfs = list(mlfs['mlf'])

    masses = an.compute_sim_masses(
    sim=simulation,
    z0halos=z0halos,
    recompute_dict=recompute,
    snap_numbers = list(mlfs['snap_num'])
    )
    list_stellar_masses = list(masses['stellar_mass'])
    list_halo_masses = list(masses['halo_mass'])
