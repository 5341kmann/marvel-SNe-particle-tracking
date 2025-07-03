import os, sys
# Add project roo to sys.path so Python can find 'base.py'
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from base import *


if __name__ == '__main__':
    sim = str(sys.argv[1])
    z0haloid = int(sys.argv[2])

    if not os.path.exists('./logs/'):
        os.mkdir('./logs/')
    logging.basicConfig(filename=f'./logs/{sim}_{z0haloid}.log', 
                        format='%(asctime)s :: %(name)s :: %(levelname)-8s :: %(message)s', 
                        datefmt='%m/%d/%Y %I:%M:%S %p',
                        level=logging.DEBUG)
    logger = logging.getLogger('PartTracker')

    logger.debug('--------------------------------------------------------------')
    logger.debug('Beginning particle tracking for {}-{}',format(sim, z0haloid))
    # in order: debug, info, warning, error

    logger.debug('Getting stored filepaths and haloids')
    filepaths, haloids, h1ids = get_stored_simpaths_haloids(sim,z0haloid)
    # filepaths starts with z=0 and goes to z=15 or so

    logger.debug('Getting starting snapshot (may take a while)')
    snap_start = get_snap_start(sim,z0haloid)
    print('Snap start: {}'.format(snap_start))
    logger.debug('Start on snapshot {}, {}'.format(snap_start, filepaths[snap_start-1][-4:]))
    
    # fix the case where the satellite doesn't have merger info prior to 
    if len(haloids) < snap_start:
        snap_start = len(haloids)
        raise Exception('Careful! You may have an error since the satellite doesnt have mergertree info out to the time where you want to start. This case is untested')
    
    if len(haloids) > snap_start:
        filepaths = np.flip(filepaths[:snap_start+1])
        haloids = np.flip(haloids[:snap_start+1])
        h1ids = np.flip(h1ids[:snap_start+1])

    if len(haloids) == snap_start:
        filepaths = np.flip(filepaths[:snap_start])
        haloids = np.flip(haloids[:snap_start])
        h1ids = np.flip(h1ids[:snap_start])   
        
    # filepaths and haloids now go the "right" way, i.e. starts from start_snap and goes until z=0
    assert len(filepaths) == len(haloids)
    assert len(haloids) == len(h1ids)

    # we save the data as an .hdf5 file since this is meant for large datasets, so that should work pretty good
    output = run_tracking(sim, z0haloid, filepaths, haloids, h1ids)

    savepath = 'Data/tracked_particles.hdf5'
    logger.debug(f'Saving output to {savepath}')
    output.to_hdf(savepath,key=f'{sim}_{z0haloid}')
