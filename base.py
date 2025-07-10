import pynbody
import pandas as pd
import numpy as np
import pickle
import logging 
import sys
import os
import time
from config import *


logger = logging.getLogger('base.py') # set up logging for debugging

os.chdir(PROGRAM_PATH) # change to the program path, so that relative paths work

# set the config to prioritize the AHF catalog
pynbody.config['halo-class-priority'] =  [pynbody.halo.ahf.AHFCatalogue,
                                          pynbody.halo.subfind.SubfindCatalogue,
                                          pynbody.halo.hop.HOPCatalogue]

def save_with_lock(df, savepath, key):
    """
    Saves a DataFrame to an HDF5 file with a file lock to prevent race conditions.
    Allows for robust parallel operations where multiple processes might try to write to the same file.
    """
    lock_path = savepath + '.lock'

    # safely makedir of savepath if it does not exist
    os.makedirs(os.path.dirname(savepath), exist_ok=True)

    # checks
    while os.path.exists(lock_path):
        # Wait for a short, random interval to avoid overwhelming the file system
        time.sleep(0.1 + 0.1 * os.getpid() % 1) 
        
    try:
        # 2. Acquire the lock by creating the lock file
        with open(lock_path, 'w') as f:
            f.write(str(os.getpid())) # Write process ID for debugging

        # Use mode 'a'(append) so we don't overwrite the entire file each time
        logger.debug(f"Process {os.getpid()} acquired lock and is saving {key}...")

        with pd.HDFStore(savepath, mode='a') as store:
            if key in store:
                # If the key already exists, remove it.
                logger.debug(f'key {key} already exists, removing it before saving.')
                store.remove(key) 

            store.put(key, df, format='table') 
            logger.debug(f"Process {os.getpid()} finished saving.")

    finally:
        # remove the lock file after saving
        if os.path.exists(lock_path):
            os.remove(lock_path)
            

# define functions for basic data manipulation, importing, etc. used by everything
def get_stored_simpaths_haloids(sim,z0haloid):
    # get snapshot paths and haloids from stored file
    with open('Data/simpaths_haloids.pickle','rb') as f:
        d = pickle.load(f)
    try:
        filepaths = d['filepaths'][sim]
    except KeyError:
        logger.debug("sim must be one of 'storm','cptmarvel','elektra','rogue'")
        raise
    try:
        haloids = d['haloids'][sim][z0haloid]
    except KeyError:
        print('z0haloid not found, perhaps this is a halo that has no stars at z=0, and therefore isnt tracked')
        raise
    # print(f'Found {len(filepaths)} filepaths: {filepaths}')
    # print(f'Found {len(haloids)} haloids: {haloids}')
    return filepaths, haloids

# function to log uncaught exceptions
def handle_exception(exc_type, exc_value, exc_traceback):
    global logger
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return
    logger.error("Uncaught exception", exc_info=(exc_type, exc_value, exc_traceback))

# sets the custom exception handler as the default
sys.excepthook = handle_exception


def get_iords(sim, z0haloid, filepaths, haloids):
    # '''Get the particle indices (iords) for all gas particles that have been in the halo since snap_start.''''
    path = f'Data/iords/{sim}_{z0haloid}.pickle'
    if os.path.exists(path):
        logger.debug(f'Found iords file at {path}, loading these')
        logger.warning(f'If you have recently changed something, these iords might not be correct')
        with open(path,'rb') as infile:
            iords = pickle.load(infile)
    
    else:
        logger.debug(f'Could not find iords file, computing iords to track')
        iords = np.array([])
        for f,haloid in zip(filepaths,haloids):
            logger.debug(f'Capturing iords for {haloid} on snapshot {f[-4:]}')
            s = pynbody.load(f)
            s.physical_units()
            h = s.halos(halo_numbers='v1')
            halo = h[haloid]
            iord = np.array(halo.gas['iord'], dtype=int)
            iords = np.union1d(iords, iord)
        
        logger.debug(f'Saving iords file to {path}')
        with open(path,'wb') as outfile:
            pickle.dump(iords,outfile)

    return iords

def run_tracking(sim, z0haloid, filepaths,haloids):
    # now we need to start tracking, so we need to get the iords
    iords = get_iords(sim, z0haloid, filepaths, haloids)
    
    use_iords = True
    output = pd.DataFrame()
    logger.debug('Starting tracking')
    for f,haloid in zip(filepaths,haloids):
        s = pynbody.load(f)
        s.physical_units()
        h = s.halos(halo_numbers='v1')
        halo = h[haloid]
        snapnum = f[-4:]
        logger.debug(f'* Snapshot {snapnum}')

        if use_iords:
            iord = np.array(s.gas['iord'],dtype=float)
            gas_particles = s.gas[np.isin(iord,iords)]
            use_iords = False
        else:
            b = pynbody.bridge.OrderBridge(s_prev,s,allow_family_change=True)
            gas_particles = b(gas_particles_prev)

        # run analysis on gas particles!
        # this calls the analysis function i've defined, and concatenates the output from this snapshot to the output overall
        output = pd.concat([output, analysis(s,halo,gas_particles,h,haloid)])

        gas_particles_prev = gas_particles
        snapnum_prev = snapnum
        s_prev = s
    
    return output


def analysis(s,halo,gas_particles,h,haloid):
    output = pd.DataFrame()
    a = float(s.properties['a'])

    if len(gas_particles) != len(gas_particles.g):
        raise Exception('Some particles are no longer gas particles...')

    # calculate properties that are invariant to centering
    output['time'] = np.array([float(s.properties['time'].in_units('Gyr'))]*len(gas_particles))
    output['pid'] = np.array(gas_particles['iord'],dtype=int)
    output['rho'] = np.array(gas_particles.g['rho'].in_units('Msol kpc**-3'), dtype=float) * 4.077603812e-8 # multiply to convert to amu/cm^3
    output['temp'] = np.array(gas_particles.g['temp'].in_units('K'), dtype=float)
    output['mass'] = np.array(gas_particles.g['mass'].in_units('Msol'), dtype=float)
    output['coolontime'] = np.array(gas_particles.g['coolontime'].in_units('Gyr'),dtype=float)
    
    # calculate properties centered on the galaxy
    pynbody.analysis.halo.center(halo)
    x,y,z = gas_particles['x'],gas_particles['y'],gas_particles['z']
    Rvir = halo.properties['Rvir'] * a / HUBBLE
    output['r'] = np.array(np.sqrt(x**2 + y**2 + z**2), dtype=float)
    output['r_per_Rvir'] = output.r / Rvir
    output['x'] = x
    output['y'] = y
    output['z'] = z
    output['Rvir'] = np.array([Rvir]*len(x))
    output['a'] = np.array([a]*len(x))

    output['vx'] = np.array(gas_particles['vx'].in_units('km s**-1'),dtype=float)
    output['vy'] = np.array(gas_particles['vy'].in_units('km s**-1'),dtype=float)
    output['vz'] = np.array(gas_particles['vz'].in_units('km s**-1'),dtype=float)
    output['v'] = np.array(np.sqrt(output.vx**2 + output.vy**2 + output.vz**2))
    

    r_half, r_gas = np.nan, np.nan

    try:
        pynbody.analysis.angmom.faceon(halo)
        
        # Calculate stellar half-mass radius
        try:
            p_stars = pynbody.analysis.profile.Profile(halo.s, bins=np.power(10, np.linspace(-2, np.log10(0.2*Rvir), 200)))
            mass_profile = p_stars['mass_enc'] / p_stars['mass_enc'][-1]
            # uses linear interpolation to find the radius at which the cumulative mass is 50%
            # changed from the previous averaging of the two surrounding bins
            r_half = np.interp(0.5, mass_profile, p_stars['rbins'])
        except (ValueError, IndexError):
            # Keep r_half as np.nan
            pass

        # Calculate gas radius based on star formation density threshold
        try:
            p_gas = pynbody.analysis.profile.Profile(halo.g, bins=np.power(10, np.linspace(-2, np.log10(Rvir), 100)))
            sigma_th = 9e6 # Msol/kpc^2 Kennicutt-Schmidt density threshold
            gas_disk_radii = p_gas['rbins'][p_gas['density'] > sigma_th]
            if len(gas_disk_radii) > 0:
                r_gas = np.max(gas_disk_radii)
        except (ValueError, IndexError):
            # Keep r_gas as np.nan
            pass

    except Exception:
        # If faceon fails, skip radius calculation
        pass

    output['r_half'] = r_half
    output['r_gas'] = r_gas

    # --- Simplified Particle Classification ---
    # A particle is either in the galaxy's ISM or the IGM.
    in_galaxy = np.isin(output.pid, halo.g['iord'])
    
    output['in_galaxy'] = in_galaxy
    output['in_IGM'] = ~in_galaxy

    return output