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
- Ensuring uniformity in data processing across different instances of computational analysis.

This module is essential for analyzing the effects of supernovae on gas dynamics within simulated environments.

Note: Ejection/expulsion code credit to Hollis Akins 2021; 
Github permalink: https://github.com/hollisakins/Justice_League_Code/blob/e049137edcfdc9838ebb3cf0fcaa4ee46e977cec/Analysis/RamPressure/analysis.py

Last revised: 3 Apr. 2022
"""
# This file gives functions for prepping, processing, and manipulation of bulk particle data. More specifically, the functions
# here retrieve gas particles satisfying certain ejection and accretion constraints and constructs respective data frames for 
# each. Said frames include arrtributes such as kinetic and potential energies, location based classifications (i.e., disk vs. 
# halo vs. field), state of supernova heating, and so on.
#
# Collecting these calculations here ensures uniformity in data processing for all instances in the computational analysis of 
# this project.
#
# ____________________________________________________________________________________________
# Ejection/Expulsion code credit to Hollis Akins 2021;
# Github permalink: https://github.com/hollisakins/Justice_League_Code/blob/ 
#                    e049137edcfdc9838ebb3cf0fcaa4ee46e977cec/Analysis/RamPressure/analysis.py
# ____________________________________________________________________________________________
# Last revised: 3 Apr. 2022

import pandas as pd
import numpy as np
import tqdm
from base import *



def read_tracked_particles(sim, haloid):
    '''
    -> Reads in gas particles tracked across a number of simulation satellites and calculates/appends desired particle 
        properties for analysis.
    '''
    #--------------------------------#
    
    logger.debug(f'Loading tracked particles for {sim}-{haloid}...')
    
    key = f'{str(sim)}_{str(int(haloid))}'

    # importing tracked particles.
    path1 = f'Data/tracked_particles.hdf5'
    data = pd.read_hdf(path1, key=key)
    logger.debug('Successfully loaded')
    
    times = np.unique(data.time)
    dt = times[1:]-times[:-1]
    dt = np.append(dt[0], dt)
    # wicked cool inverse array trick which matches the dt array to the original data
    dt = dt[np.unique(data.time, return_inverse=True)[1]]
    data['dt'] = dt
    
    
    
    # calcuatiing r_gal to be the larger of either the stellar body (r_half) or the gas body (r_gas).
    # using the faster elementwise assignment of r_gal.
    data['r_gal'] = np.maximum(data['r_half'], data['r_gas'])

    # filling forward NaN values in 'r_gal' to ensure all particles have a valid r_gal value.
    data['r_gal'] = data.groupby('pid')['r_gal'].transform('ffill')
    
    thermo_disk = (np.array(data.temp) < 1.2e4) & (np.array(data.rho) > 0.1)
    
    in_gal = np.array(data.in_galaxy)
    
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
    ejected = pd.DataFrame()
    cooled = pd.DataFrame()
    expelled = pd.DataFrame()
    accreted = pd.DataFrame()
    
    pids = np.unique(data.pid)
    # ------------ This was too slow -------------
    # # comparing particle at previous step with current step to determine ejection/expulsion.
    # # time sorted time-series analysis for more robust relative ordering for particle classifcation
    # for pid in tqdm.tqdm(pids, desc="Processing particles"):
    #     dat = data[data.pid==pid].sort_values('time').reset_index(drop=True)

    #     gal_disk = np.array(dat.gal_disk, dtype=bool)
    #     gal_halo = np.array(dat.gal_halo, dtype=bool)
    #     in_gal = np.array(dat.in_galaxy, dtype=bool)


    #     for i in range(1, len(dat)):
    #             # particle ejection
    #             if gal_disk[i-1] and gal_halo[i]:
    #                 ejected = pd.concat([ejected, dat.iloc[[i]]])

    #             # particle cooling
    #             if gal_halo[i-1] and gal_disk[i]:
    #                 cooled = pd.concat([cooled, dat.iloc[[i]]])

    #             # particle expulsion    
    #             if in_gal[i-1] and not in_gal[i]:
    #                 row_frame = data.iloc[[i]].copy()
    #                 # more descriptive orign classification
    #                 row_frame['expulsion_origin'] = 'gal_disk' if gal_disk[i-1] else 'gal_halo'
    #                 expelled = pd.concat([expelled, row_frame])
                
    #             # particle accretion  
    #             if not in_gal[i-1] and in_gal[i]:
    #                 row_frame = dat.iloc[[i]].copy()
    #                 row_frame['accreation_dest'] = 'gal_disk' if gal_disk[i] else 'gal_halo'
    #                 accreted = pd.concat([accreted, row_frame])
    
    # key = f'{sim}_{haloid}'

    # filepath = 'Data/SNe/ejected_particles.hdf5'
    # save_with_lock(ejected, filepath, key=key)
    
    # filepath = 'Data/SNe/cooled_particles.hdf5'
    # save_with_lock(cooled, filepath, key=key)

    # filepath = 'Data/SNe/expelled_particles.hdf5'
    # save_with_lock(expelled, filepath, key=key)
            
    # filepath ='Data/SNe/accreted_particles.hdf5'
    # save_with_lock(accreted, filepath, key=key)
    # --------------------------------

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

    # removing temporary columns which have been used for calculations 
    logger.debug('Dropping temporary columns used for calculations.')
    ejected = ejected.drop(columns=['prev_gal_disk', 'prev_gal_halo', 'prev_in_galaxy'])
    cooled = cooled.drop(columns=['prev_gal_disk', 'prev_gal_halo', 'prev_in_galaxy'])
    expelled = expelled.drop(columns=['prev_gal_disk', 'prev_gal_halo', 'prev_in_galaxy'])
    accreted = accreted.drop(columns=['prev_gal_disk', 'prev_gal_halo', 'prev_in_galaxy'])

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
        #if cooling is turned off, particle is SNe heated
        #checks state of particle prior to outflow
        is_sn_heated = np.array(pre_outflow_states.coolontime) > np.array(pre_outflow_states.time)
        disk_outflows = disk_outflows.reset_index(drop=True)  # Reset index for consistency
        disk_outflows['snHeated'] = is_sn_heated
        logger.debug(f'Outflows classified: {disk_outflows["snHeated"].sum()} particles are SN-heated.')

    # removing temporary columns used for calculations
    logger.debug('Dropping temporary columns used for calculations.')
    pre_outflow_states = pre_outflow_states.drop(columns=['prev_gal_disk'], errors='ignore')
    disk_outflows = disk_outflows.drop(columns=['prev_gal_disk'], errors='ignore')
    disk_inflows = disk_inflows.drop(columns=['prev_gal_disk'], errors='ignore')


    key = f'{sim}_{haloid}'
    # Use our robust save_with_lock function for safety
    logger.debug(f'Saving results for {key}...')
    save_with_lock(pre_outflow_states, 'Data/SNe/pre_outflow_states.hdf5', key)
    save_with_lock(disk_outflows, 'Data/SNe/disk_outflows.hdf5', key)
    save_with_lock(disk_inflows, 'Data/SNe/disk_inflows.hdf5', key)

def calc_hot_predischarged(sim, haloid, save=True, verbose=True):
    '''
    -> Identifies discharged gas particles that experienced supernova heating just prior to being 
        discharged. Properties prior to discharge are recorded in 'hot_predischarged'. 
    '''
    #--------------------------------#
    
    import tqdm
    data = read_tracked_particles(sim, haloid, verbose=verbose)
    
    logger.debug(f'Now compiling hot_predischarged particles for {sim}-{haloid}...')
    
    hot_predischarged = pd.DataFrame() # properties pre-discharge for heated gas.
    
    
    pids = np.unique(data.pid)
    for pid in tqdm.tqdm(pids):
        dat = data[data.pid==pid]

        gal_disk = np.array(dat.gal_disk, dtype=bool)
        in_gal = np.array(dat.in_galaxy, dtype=bool)
        outside_disk = ~gal_disk
        
        time = np.array(dat.time, dtype=float)
        coolontime = np.array(dat.coolontime, dtype=float)


        for i,t2 in enumerate(time[1:]):
                i += 1
                if gal_disk[i-1] and outside_disk[i] and (coolontime[i] > time[i-1]): #discharged and SN-heated
                    out = dat[time==time[i-1]].copy()
                    hot_predischarged = pd.concat([hot_predischarged, out])
                 
    # apply the calc_angles function along the rows of discharged particles.
    logger.debug('Calculating hot_predischarged angles;')
    hot_predischarged = hot_predischarged.apply(calc_angles, axis=1)
   
    
    if save:
        key = f'{sim}_{str(int(haloid))}'
        filepath = f'{rootPath}SNe-heated_Gas_Flow/SNe/hot_predischarged_particles.hdf5'
        logger.debug(f'Saving {key} pre-dsrg, SN-heated particles to {filepath}')
        hot_predischarged.to_hdf(filepath, key=key)
        
    logger.debug(f'> Returning (hot_predischarged) <')
    return hot_predischarged


def calc_reaccreted(sim, haloid, save=True, verbose=True):
    ''' 
    -> 'Advanced' computation of accreted gas particles denoted 'reaccreted'.
    -> Screening the 'accreted' df compiled by 'calc_discharge()' specifically for gas particles 
        previously discharged from their satellite's disk, and which are accreted (reaccreted) back onto 
            the disk at a later timestep. 
    -> (Only particles with an accretion event that has a directly preceeding discharge event are 
        compiled into 'reaccreted'.)
    '''
    #--------------------------------#
    
    import tqdm
    key = f'{sim}_{str(int(haloid))}'

    path = f'{rootPath}SNe-heated_Gas_Flow/SNe/discharged_particles.hdf5'
    discharged = pd.read_hdf(path, key=key)
    path = f'{rootPath}SNe-heated_Gas_Flow/SNe/accreted_particles.hdf5'
    accreted = pd.read_hdf(path, key=key)


    logger.debug(f'Now computing reaccreted particles for {sim}-{haloid}...')
    reaccreted = pd.DataFrame() # gas accreted following a discharge event.
  
    # defining attribute giving the length of time between discharge and accretion event for each gas particle:
    recycleTime = {'recycleTime': ''} 
    accreted = accreted.join(pd.DataFrame(columns=recycleTime)) 
    # ensuring that our new accreted dataframe inherits sne heating identified 'hot'.
    heating = {'snHeated': ''} 
    accreted = accreted.join(pd.DataFrame(columns=heating))

    pids = np.unique(discharged.pid) # quicker to use 'discharged' because fewer unique particles.
    for pid in tqdm.tqdm(pids):
        dis = discharged[discharged.pid==pid]
        acc = accreted[accreted.pid==pid]

        dTime = np.asarray(dis.time)
        aTime = np.asarray(acc.time)


        # removing initial accretion event if it does not correspond to a discharge; else no alterations.
        # check to ensure that our discharge event actually has an accretion:
        if (len(aTime) == 0) or (len(dTime) == 0): # if no instances of accretion or none of discharge, move on to next particle.
            continue
        if (aTime[0] < dTime[0]):
            aCache = acc[1:]

        else:
            aCache = acc

        if len(aCache) == 0: # if no instances of reaccretion, move on to next particle.
            continue


        dCache = dis[0:len(aCache)]

        import warnings
        pd.options.mode.chained_assignment = None #this will supress SettingWithCopyWarning
        with warnings.catch_warnings():
            warnings.simplefilter(action='ignore', category=FutureWarning)
            recyclingTime = np.array(aCache['time']).max() - np.array(dCache['time']).max()
            aCache.loc['recycleTime'] = recyclingTime

        heated = np.array(dCache['snHeated']).max()
        aCache.loc['snHeated'] = heated
        reaccreted = pd.concat([reaccreted, aCache])

        pd.options.mode.chained_assignment = "raise"

    if save:
        key = f'{sim}_{str(int(haloid))}'
        filepath = f'{rootPath}SNe-heated_Gas_Flow/SNe/reaccreted_particles.hdf5'
        logger.debug(f'Saving {key} reaccreted particle dataset to {filepath}')
        reaccreted.to_hdf(filepath, key=key)
        
    logger.debug(f'> Returning (reaccreted) dataset <')
    return reaccreted


def calc_expelled(key, save=True, verbose=True):
    #find expelled (permanently) gas
    #define gas particles as 'expelled' if they are not reaccreted onto the disk of their host satellite
    #after being discharged
    
    logger.debug(f'Now compiling expelled particles for {key}...')
    
    #load discharged, reaccreted gas particles
    predischarged,discharged = read_one_discharged(key)
    _,reaccreted = read_one_accreted(key)
    
    #create a dataframe for expelled gas
    expelled = pd.DataFrame()
    preexpelled = pd.DataFrame()
    
    #pid as np.array
    did = np.array(discharged.pid)
    rid = np.array(reaccreted.pid)
    
    #drop the duplicates of pids and get unique set
    dfunique = discharged.drop_duplicates(subset=['pid'], keep='last')
    
    #find the pid expelled
    index = np.argsort(rid)
    sorted_rid = rid[index]  # Sorted list of ids reaccreted

    from collections import Counter
    #find ids that discharged more than reaccreted (expelled at the last timestep)
    #also find pids that discharged but never be reaccreted
    subCount = Counter(did) - Counter(sorted_rid)
    dunique = np.array(list(subCount.items()))[:,0]
    expelled = dfunique[np.isin(dfunique.pid, dunique, assume_unique=True)]
    
    #do the same for predischarged
    pre_dfunique = predischarged.drop_duplicates(subset=['pid'], keep='last')
    preexpelled = pre_dfunique[np.isin(dfunique.pid, dunique, assume_unique=True)]
    if save:
        filepath = f'{rootPath}SNe-heated_Gas_Flow/SNe/expelled_particles.hdf5'
        logger.debug(f'Saving {key} expelled particles to {filepath}')
        expelled.to_hdf(filepath, key=key)
        
        filepath = f'{rootPath}SNe-heated_Gas_Flow/SNe/preexpelled_particles.hdf5'
        logger.debug(f'Saving {key} preexpelled particles to {filepath}')
        preexpelled.to_hdf(filepath, key=key)
        
    logger.debug(f'> Returning (predischarged, discharged, accreted) datasets <')
    return preexpelled, expelled


def calc_snHeated(particles):
    """
    This function detects the gas particles that are sn-heated by comparing the coolontime 
    with the time 1 timestep before. Since the detection of sn-heating doesn't depend on any positional argument,
    this includes gas particles being SN-heated at any point. However, in order to avoid counting gas particles
    that are outside the satellite, as those are likely not heated by the satellite, 
    we should constrain only within the satellite.
    """
    import tqdm
    #iterate detection process by pids
    pids = np.unique(particles['pid'])
    index = np.array([]) #initialize
    for pid in tqdm.tqdm(pids):
        data = particles[particles['pid']==pid]
        #create a structured array, containing index of dataframe, pid, time, and coolontime
        dtype = [('index', int), ('pid', int), ('time', float), ('coolontime', float)]
        structureArray=np.array(list(zip(data.index, *map(data.get, ['pid','time','coolontime']))), dtype=dtype)
        #limit to after being heated (avoid mistakingly take the row heated at the same timestep)
        heatedArray=structureArray[structureArray['time']>structureArray['coolontime']]
        #extract the list of unique coolontime, sorted by pid and time
        helper1, helper2 = np.unique(heatedArray['coolontime'], return_index = True)
        heatedunique = np.sort(heatedArray[helper2], order=['pid','time'])

        timebefore = heatedunique['time'][:-1]
        #find sn-heated list by comparing the time before and coolontime
        heatedLocal = heatedunique[1:][heatedunique['coolontime'][1:]>timebefore]
        indexLocal = heatedLocal['index'].astype(int)
        index = np.append(index, indexLocal)
    #based on detected indices, find the designated rows from original dataframe
    heated = particles[particles.index.isin(index)] 
    #drop if gas is outside satellite (e.g. host, other sat, IGM)
    heated = heated[heated['in_gal']==True]

    return heated


def calc_snGas(sim, haloid, save=True, verbose=True):
    '''
    -> I think this calculation is wrong because it's comparing coolontime with time before for all the particles discharged or not.
    -> July 2025: UPDATED: Compares coolontime with the current time, so that it only returns particles that are SN-heated at the current time step
    '''
    #--------------------------------#
    
    import tqdm
    data = read_tracked_particles(sim, haloid, verbose=verbose)
    
    logger.debug(f'Now compiling SN-heated gas for {sim}-{haloid}...')
        
    sngas = pd.DataFrame() # all gas in sims that experienced SN-heating.
    
    pids = np.unique(data.pid)
    for pid in tqdm.tqdm(pids):
        dat = data[data.pid==pid]

        time = np.array(dat.time, dtype=float)
        coolontime = np.array(dat.coolontime, dtype=float)
        
        for i,t in enumerate(time):
            if (coolontime[i] > time[i]):
                hot = dat[time==t].copy()
                sngas = pd.concat([sngas, hot])
    
    key = f'{sim}_{str(int(haloid))}'
    filepath = 'Data/SNe/sngas_particles.hdf5'
    logger.debug(f'Saving {key} SN-heated particles to {filepath}')
    sngas.to_hdf(filepath, key=key)

def read_all_ejected_expelled():
    '''
    -> Reads ejected, cooled, expelled, and accreted into workable dataframes for analysis in notebooks.
    '''
    #--------------------------------#
    
    ejected = pd.DataFrame()
    cooled = pd.DataFrame()
    expelled = pd.DataFrame()
    accreted = pd.DataFrame()
    keys = get_keys()
    for key in keys:
        if key in ['h148_3','h148_28','h242_12']: continue;
            
        ejected1 = pd.read_hdf(f'{rootPath}SNe-heated_Gas_Flow/SNe/ejected_particles.hdf5', key=key)
        ejected1['key'] = key
        ejected = pd.concat([ejected, ejected1])
        cooled1 = pd.read_hdf(f'{rootPath}SNe-heated_Gas_Flow/SNe/cooled_particles.hdf5', key=key)
        cooled1['key'] = key
        cooled = pd.concat([cooled, cooled1])
        expelled1 = pd.read_hdf(f'{rootPath}SNe-heated_Gas_Flow/SNe/expelled_particles.hdf5', key=key)
        expelled1['key'] = key
        expelled = pd.concat([expelled, expelled1])
        accreted1 = pd.read_hdf(f'{rootPath}SNe-heated_Gas_Flow/SNe/accreted_particles.hdf5', key=key)
        accreted1['key'] = key
        accreted = pd.concat([accreted, accreted1])

    logger.debug(f'> Returning (ejected, cooled, expelled, accreted) for all satellites <')
    return ejected, cooled, expelled, accreted


def read_all_discharged():
    '''
    -> Reads predischarged, discharged, accreted, and hot_predischarged into workable dataframes for
        analysis in notebooks.
    '''
    #--------------------------------#
    
    predischarged = pd.DataFrame()
    discharged = pd.DataFrame()
    # hot_predischarged= pd.DataFrame()
    
    keys = get_keys()

    for i,key in enumerate(keys):
        i += 1
        if key == 'h242_401':
            continue

        sim = key[:4]
        haloid = int(key[5:])
        predischarged1 = pd.read_hdf(f'{rootPath}SNe-heated_Gas_Flow/SNe/predischarged_particles.hdf5', key=key)
        predischarged1['key'] = key
        predischarged = pd.concat([predischarged, predischarged1])
        
        discharged1 = pd.read_hdf(f'{rootPath}SNe-heated_Gas_Flow/SNe/discharged_particles.hdf5', key=key)
        discharged1['key'] = key
        discharged = pd.concat([discharged, discharged1])
  
        # hot_predischarged1 = pd.read_hdf(f'{rootPath}SNe-heated_Gas_Flow/SNe/hot_predischarged_particles.hdf5', key=key)
        # hot_predischarged1['key'] = key
        # hot_predischarged = pd.concat([hot_predischarged, hot_predischarged1])
       
    logger.debug(f'> Returning (predischarged, discharged) for all satellites <')
    return predischarged, discharged


def read_accreted():
    '''
    -> Reads all accreted particles, reaccreted particles into workable dataframes for analysis.
    '''
    #--------------------------------#
    
    accreted = pd.DataFrame()
    reaccreted = pd.DataFrame()

    keys = get_keys()

    for i,key in enumerate(keys):
        i += 1

        if key == 'h242_401':
            continue

        sim = key[:4]
        haloid = int(key[5:])
        accreted1 = pd.read_hdf(f'{rootPath}SNe-heated_Gas_Flow/SNe/accreted_particles.hdf5', key=key)
        accreted1['key'] = key
        accreted = pd.concat([accreted, accreted1])
        
        reaccreted1 = pd.read_hdf(f'{rootPath}SNe-heated_Gas_Flow/SNe/reaccreted_particles.hdf5', key=key)
        reaccreted1['key'] = key
        reaccreted = pd.concat([reaccreted, reaccreted1])

    logger.debug(f'> Returning (accreted, reaccreted) for all satellites <')
    return accreted, reaccreted


def read_expelled():
    '''
    -> Reads all expelled particles, preexpelled particles into workable dataframes for analysis.
    '''
    #--------------------------------#
    
    expelled = pd.DataFrame()
    preexpelled = pd.DataFrame()

    keys = get_keys()

    for i,key in enumerate(keys):
        i += 1

        if key == 'h242_401':
            continue
        
        sim = key[:4]
        haloid = int(key[5:])
        expelled1 = pd.read_hdf(f'{rootPath}SNe-heated_Gas_Flow/SNe/expelled_particles.hdf5', key=key)
        expelled1['key'] = key
        expelled = pd.concat([expelled, expelled1])
        
        preexpelled1 = pd.read_hdf(f'{rootPath}SNe-heated_Gas_Flow/SNe/preexpelled_particles.hdf5', key=key)
        preexpelled1['key'] = key
        preexpelled = pd.concat([preexpelled, preexpelled1])

    logger.debug(f'> Returning (expelled, preexpelled) for all satellites <')
    return preexpelled, expelled

    
def read_sngas():
    '''
    -> Reads all gas particles in selected satellites ever SN-heated (irrespective of whether or not
        they were discharged.
    Starting Summer 23, this function was found to be not used for accurately calculating SN-heated gas. Therefore, use `sneHeated==True` instead.
    '''
    #--------------------------------#
    
    sntotal = pd.DataFrame()

    keys = get_keys()

    for i,key in enumerate(keys):
        i += 1
        sim = key[:4]
        haloid = int(key[5:])
        sntotal1 = pd.read_hdf(f'{rootPath}SNe-heated_Gas_Flow/SNe/sngas_particles.hdf5',
                               key=key)
        sntotal1['key'] = key
        sntotal = pd.concat([sntotal, sntotal1])

    logger.debug(f'> Returning (SN-heated gas) for all satellites <')
#     return sntotal


def read_one_discharged(key):
    '''
    -> Reads predischarged, discharged, accreted, and hot_predischarged into workable dataframes.
        Only run for one sat.
    '''
    #--------------------------------#
    
    predischarged = pd.DataFrame()
    discharged = pd.DataFrame()
    #hot_predischarged= pd.DataFrame()
    
    sim = key[:4]
    haloid = int(key[5:])
    predischarged = pd.read_hdf(f'{rootPath}SNe-heated_Gas_Flow/SNe/predischarged_particles.hdf5', key=key)

    discharged = pd.read_hdf(f'{rootPath}SNe-heated_Gas_Flow/SNe/discharged_particles.hdf5', key=key)

    #hot_predischarged = pd.read_hdf(f'{rootPath}SNe-heated_Gas_Flow/SNe/hot_predischarged_particles.hdf5', key=key)
       
    logger.debug(f'> Returning (predischarged, discharged, hot_predischarged) for satellite {key} <')
    return predischarged, discharged#, hot_predischarged


def read_one_accreted(key):
    '''
    -> Reads accreted particles, reaccreted particles into workable dataframes for analysis.
        Only run for one sat.
    '''
    #--------------------------------#
    
    accreted = pd.DataFrame()
    reaccreted = pd.DataFrame()

    sim = key[:4]
    haloid = int(key[5:])

    accreted = pd.read_hdf(f'{rootPath}SNe-heated_Gas_Flow/SNe/accreted_particles.hdf5', key=key)
    
    reaccreted = pd.read_hdf(f'{rootPath}SNe-heated_Gas_Flow/SNe/reaccreted_particles.hdf5', key=key)

    logger.debug(f'> Returning (accreted, reaccreted) for satellite {key} <')
    return accreted, reaccreted


def read_one_expelled(key):
    
    expelled = pd.DataFrame()
    preexpelled = pd.DataFrame()
    expelled = pd.read_hdf(f'{rootPath}SNe-heated_Gas_Flow/SNe/expelled_particles.hdf5', key=key)
    preexpelled = pd.read_hdf(f'{rootPath}SNe-heated_Gas_Flow/SNe/preexpelled_particles.hdf5', key=key)

    return preexpelled, expelled


logger.debug("compiler.py executed")
