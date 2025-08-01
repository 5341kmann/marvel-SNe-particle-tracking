import os, sys
# Add project root to sys.path so Python can find 'base.py'
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
    logger.debug('Beginning particle tracking for {}-{}'.format(sim, z0haloid))
    # in order: debug, info, warning, error

    logger.debug('Getting stored filepaths and haloids')
    filepaths, haloids = get_stored_simpaths_haloids(sim,z0haloid)
    # filepaths starts with z=0 and goes to z=15 or so

    logger.debug('Start on snapshot: {}'.format(filepaths[0][-4:]))
        
    # # filepaths and haloids now go the "right" way, i.e. starts from start_snap and goes until z=0
    # assert len(filepaths) >= len(haloids)

    # we save the data as an .hdf5 file since this is meant for large datasets, so that should work pretty good
    output = run_tracking(sim, z0haloid, filepaths, haloids)

    savepath = 'Data/tracked_particles.hdf5'
    logger.debug(f'Saving output to {savepath}')
    save_with_lock(output, savepath, key=f'{sim}_{z0haloid}')
   
