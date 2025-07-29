import os, sys
# Add project roo to sys.path so Python can find 'base.py'
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from base import *
from compiler import *


if __name__ == '__main__':
    sim = str(sys.argv[1])
    z0haloid = int(sys.argv[2])

    if not os.path.exists('./logs/'):
        os.mkdir('./logs/')
    logging.basicConfig(filename=f'./logs/{sim}_{z0haloid}.log', 
                        format='%(asctime)s :: %(name)s :: %(levelname)-8s :: %(message)s', 
                        datefmt='%m/%d/%Y %I:%M:%S %p',
                        level=logging.DEBUG)
    logger = logging.getLogger('PartCalculator')

    logger.debug('--------------------------------------------------------------')
    logger.debug(f'Beginning particle calculations for {sim}_{z0haloid}')


    data = read_tracked_particles(sim, z0haloid)
    # logger.debug('Read tracked particles data')

    logger.debug('Calculating ejected and expelled particles')
    calc_ejected_expelled(data, sim, z0haloid)

    # logger.debug('Calculating disk outflows')
    # calc_disk_outflows(data, sim, z0haloid)

    # logger.debug('Calculating permantly expelled particles')
    # calc_perm_expelled(f'{sim}_{z0haloid}')
    # in order: debug, info, warning, error







    





                
