"""
This module provides functions for the preparation, processing, and manipulation of bulk particle data in the context of 
gas flow simulations influenced by supernovae. The primary focus is on retrieving gas particles that meet specific 
ejection and accretion criteria, and constructing data frames that encapsulate various attributes such as kinetic and 
potential energies, classifications based on spatial location (e.g., disk, halo, field), and the state of supernova 
heating.

Key functionalities include:
- Reading tracked particles from simulation data.
- Calculating properties of ejected, expelled, cooled, and accreted particles.
- Identifying discharged particles and their properties before and after discharge.

Note: Original Ejection/expulsion code adapted from Hollis Akins 2021; 
Github permalink: https://github.com/hollisakins/Justice_League_Code/blob/e049137edcfdc9838ebb3cf0fcaa4ee46e977cec/Analysis/RamPressure/analysis.py

Last revised: 26 Jul. 2025
"""

import pandas as pd
import numpy as np
from base import *
from analysis import JusticeLeagueSim, MarvelSim


def read_tracked_particles(sim, haloid):
    '''
    -> Reads in gas particles tracked across a number of simulation satellites and calculates/appends desired particle 
        properties for analysis.
    '''
    #--------------------------------#
    
    logger.debug(f'Loading tracked particles for {sim}-{haloid}...')
    
    key = f'{str(sim)}_{str(int(haloid))}'
    marvel_sim = True

    # importing tracked particles.
    if sim in [member.value for member in MarvelSim]:
        path = f'Data/tracked_particles.hdf5'
    elif sim in [member.value for member in JusticeLeagueSim]:
        path = JL_TRACKED_PARTICLES_PATH
        marvel_sim = False
    else:
        logger.error(f'Simulation {sim} not recognized. Please use a valid simulation name.')
        return None
    
    data = pd.read_hdf(path, key=key)
    logger.debug('Successfully loaded')
    
    times = np.unique(data.time)
    dt = times[1:]-times[:-1]
    dt = np.append(dt[0], dt)
    # wicked cool inverse array trick which matches the dt array to the original data
    dt = dt[np.unique(data.time, return_inverse=True)[1]]
    data['dt'] = dt
    
    if not marvel_sim:
        data['r_half'] = data['sat_r_half']
        data['r_gas'] = data['sat_r_gas']

        # dropping the satellite radius columns to avoid confusion.
        data = data.drop(columns=['sat_r_half', 'sat_r_gas'], errors='ignore')
    
    # calcuatiing r_gal to be the larger of either the stellar body (r_half) or the gas body (r_gas).
    # using the faster elementwise assignment of r_gal.
    data['r_gal'] = np.maximum(data['r_half'], data['r_gas'])
    # filling forward NaN values in 'r_gal' to ensure all particles have a valid r_gal value.
    data['r_gal'] = data.groupby('pid')['r_gal'].transform('ffill')
    
    thermo_disk = (np.array(data.temp) < 1.2e4) & (np.array(data.rho) > 0.1)
    
    # check for in case of processing satellite data.
    if 'in_galaxy' in data.columns:
        in_gal = np.array(data.in_galaxy) 
    else:
        # assumed to Justice League data, where 'in_sat' is used instead.
        in_gal = np.array(data.in_sat)
        # mending the column name to match the expected format.
        data['in_galaxy'] = data['in_sat']
        data = data.drop(columns=['in_sat'], errors='ignore')
    
    gal_disk = in_gal & thermo_disk
    gal_halo = in_gal & ~thermo_disk
    
    # basic location classifications.
    data['gal_disk'] = gal_disk
    data['gal_halo'] = gal_halo

    logger.debug('> Returning <tracked_particle> dataset <')
    return data


def calc_ejected_expelled(data, sim, haloid):
    '''
    -> Identifies gas particles meeting 'ejection' and 'expulsion' criteria, as well as those that have been cooled and
        reaccreted by their respective satellites.
    '''
    #--------------------------------#
#
    logger.debug(f'Now computing ejected/expelled particles for {sim}-{haloid}...')

    # Using vectorized operations for greater efficiency
    # 1. Sort the data by 'pid' and 'time' to ensure chronological order for each particle
    data = data.sort_values(['pid', 'time']).reset_index(drop=True)

    # 2. Creates previous state columns for each particle's history
    data['prev_gal_disk'] = data.groupby('pid')['gal_disk'].shift(1)
    data['prev_gal_halo'] = data.groupby('pid')['gal_halo'].shift(1)
    data['prev_in_galaxy'] = data.groupby('pid')['in_galaxy'].shift(1)
    logger.debug('Previous state columns created for gal_disk, gal_halo, and in_galaxy.')

    # 3. Identify all events at once using boolean masks
    # Ejection
    ejected_mask = (data['prev_gal_disk'] == True) & (data['gal_halo'] == True)
    logger.debug(f'Ejection mask created with {ejected_mask.sum()} particles ejected.')

    # Cooling
    cooled_mask = (data['prev_gal_halo'] == True) & (data['gal_disk'] == True)
    logger.debug(f'Cooling mask created with {cooled_mask.sum()} particles cooled.')
    # Expulsion
    expelled_mask = (data['prev_in_galaxy'] == True) & (data['in_galaxy'] == False)
    logger.debug(f'Expulsion mask created with {expelled_mask.sum()} particles expelled.')

    # Accretion
    accreted_mask = (data['prev_in_galaxy'] == False) & (data['in_galaxy'] == True)
    logger.debug(f'Accretion mask created with {accreted_mask.sum()} particles accreted.')

    # Added handling for NaN values from shift to ensure proper saving with format='table' ind HDF5
    for col in ['prev_gal_disk', 'prev_gal_halo', 'prev_in_galaxy']:
        if col in data.columns:
            data[col] = data[col].astype(float)
    logger.debug('Converted "prev_" boolean columns to float dtype for HDF5 compatibility.')

    # 4. Create the final DataFrames by applying the masks
    ejected = data[ejected_mask].copy()
    cooled = data[cooled_mask].copy()
    expelled = data[expelled_mask].copy()
    accreted = data[accreted_mask].copy()

    # 5. Adding origin/destination columns 
    if not expelled.empty:
        expelled['expulsion_origin'] = np.where(expelled['prev_gal_disk'], 'gal_disk', 'gal_halo')

    if not accreted.empty:
        accreted['accretion_destination'] = np.where(accreted['gal_disk'], 'gal_disk', 'gal_halo')
        
    # 6. Save the results
    key = f'{sim}_{haloid}'
    logger.debug(f'Saving results for {key}...')
    save_with_lock(ejected, 'Data/SNe/ejected_particles.hdf5', key)
    save_with_lock(cooled, 'Data/SNe/cooled_particles.hdf5', key)
    save_with_lock(expelled, 'Data/SNe/expelled_particles.hdf5', key)
    save_with_lock(accreted, 'Data/SNe/accreted_particles.hdf5', key)
    logger.debug('Ejected, cooled, expelled, and accreted particles calculated and saved successfully.')

def calc_disk_outflows(data, sim, haloid):
    '''
    -> Identifies particles that flow out of and into the galaxy disk.
    -> Specifically flags outflowing particles that were recently heated,
       likely by supernovae, prior to leaving the disk.
       ADAPTED FOR ISOLATED GALAXIES.
    '''
    logger.debug(f'Now identifying disk outflows for {sim}-{haloid}...')

    data = data.sort_values(['pid', 'time']).reset_index(drop=True)

    data['prev_gal_disk'] = data.groupby('pid')['gal_disk'].shift(1)

    # 3. Identify all outflow and inflow events simultaneously using boolean masks
    outflow_mask = (data['prev_gal_disk'] == True) & (data['gal_disk'] == False)
    logger.debug(f'Outflow mask created with {outflow_mask.sum()} particles outflowing.')
    inflow_mask = (data['prev_gal_disk'] == False) & (data['gal_disk'] == True)
    logger.debug(f'Inflow mask created with {inflow_mask.sum()} particles inflowing.')

    disk_outflows = data[outflow_mask].copy()
    disk_inflows = data[inflow_mask].copy()

    # Find the indices of the outflow events and subtract 1 to get the state just before.
    pre_outflow_indices = disk_outflows.index - 1
    pre_outflow_states = data.loc[pre_outflow_indices].copy()

    # ---  Identify Supernova-Heated Outflows ---
    logger.debug('Classifying supernova-heated subset of outflows...')
    if not disk_outflows.empty:
        # reseting indexes to ensure alignment
        disk_outflows = disk_outflows.reset_index(drop=True)
        pre_outflow_states = pre_outflow_states.reset_index(drop=True)

        #classifying SNe heated as having both
        # coolontime greater than the time at the previous step.
        # previous coolontime greater than the time at that timestep
        #  - capuring SNe activity occuring between two timestep while excluding activity begining prior to previous timestep
        condition1 = disk_outflows['coolontime'] > pre_outflow_states['time']
        condition2 = pre_outflow_states['coolontime'] < pre_outflow_states['time']
        
        # Combine the conditions to find particles heated between the two snapshots.
        is_sn_heated = condition1 & condition2

        # --- Assign the final boolean flag ---
        disk_outflows['snHeated'] = is_sn_heated

    # removing temporary columns used for calculations
    logger.debug('Dropping temporary columns used for calculations.')
    pre_outflow_states = pre_outflow_states.drop(columns=['prev_gal_disk'], errors='ignore')
    disk_outflows = disk_outflows.drop(columns=['prev_gal_disk'], errors='ignore')
    disk_inflows = disk_inflows.drop(columns=['prev_gal_disk'], errors='ignore')


    key = f'{sim}_{haloid}'

    # Use robust save_with_lock function for safety
    logger.debug(f'Saving results for {key}...')
    save_with_lock(pre_outflow_states, 'Data/SNe/pre_outflow_states.hdf5', key)
    save_with_lock(disk_outflows, 'Data/SNe/disk_outflows.hdf5', key)
    save_with_lock(disk_inflows, 'Data/SNe/disk_inflows.hdf5', key)

def calc_perm_expelled(key):

    # checks for the existence of the necessary hdf5 files.
    if not os.path.exists('Data/SNe/expelled_particles.hdf5'):
        logger.error('The expelled particles file has not yet been produced.\n\
                     Please populate inflows at Data/SNe/expelled_particles.hdf5 using calc_ejected_expelled()')
        return
    if not os.path.exists('Data/SNe/accreted_particles.hdf5'):
        logger.error('The accreated particles fil has not yet been produced.\n\
                     Please populate outflows at Data/SNe/accreted_particles.hdf5 using calc_ejected_expelled()')
        return

    # reads in the disk inflows and outflows, concatenating them into a single dataframe and adds event_type column classifying the events as either inflow or outflow.
    particle_flows = pd.concat([pd.read_hdf('Data/SNe/accreted_particles.hdf5', key=key).assign(event_type = 'inflow'), 
                              pd.read_hdf('Data/SNe/expelled_particles.hdf5', key=key).assign(event_type = 'outflow')])

    particle_flows = particle_flows.sort_values(['pid', 'time'])

    last_cross = particle_flows.drop_duplicates(subset=['pid'], keep='last')

    permanent_expelled = last_cross[last_cross['event_type'] == 'outflow']

    # saving the permanent expelled particles to an hdf5 file.
    save_with_lock(permanent_expelled, 'Data/SNe/permanent_expelled_particles.hdf5', key=key)