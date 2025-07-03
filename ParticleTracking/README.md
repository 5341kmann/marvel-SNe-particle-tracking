# Ram Pressure Stripping

In this directory lies code I've used for my Spring 2025 MAP, which has been adapted from Hollis Akin's 2021 Justice_League_Code
The goal of this project is to use particle tracking code to analyze...

## Base Code 
> `base.py`

This file stores basic data manipulation functions: getting the filepaths and haloids of the galaxies of interest, and more.


## Scripts

> `particletracking.py`
- This script runs particle tracking on a particular satellite in a particular simulation, which are specified at the command line (e.g. `python particletracking.py h329 11`). 
- The code tracks gas particles, starting at the snapshot where the satellite first crosses 2 Rvir from its host and ending at redshift 0. Gas particles that are tracked are *those that are in the satellite for at least one snapshot in this range*.
- Uses as input data: simulation snapshots, `Data/simpaths_haloids.pickle`.
- Produces as output data: `Data/iords/sim_haloid.pickle`, `Data/tracked_particles.hdf5`.
- Note: the bash script `runall.sh` can be used to run particle tracking on multiple satellites at once, to speed things along.  

## Data

The `particletracking.py` script draws data directly from the simulation snapshots. To speed up the process of analyzing these satellites over time, I have stored the simulation snapshot filepaths and main progenitor haloids for each redshift 0 satellite at `Data/simpaths_haloids.pickle`. The scripts in this directory utilize this pickle file to get satellite haloid information. 

Most of the datasets I've created as part of this project are in the form of `.hdf5` files. 
The HDF5 file format stores the output data efficiently and allows for data from all the galaixes to be stored in one file. The HDF5 file type can be read into python easily as a `pandas` DataFrame using `pandas.read_hdf(filepath, key=key)`, where the `key` option specifies the table you want. Each galaxy is saved with a key identifying its host and redshift 0 haloid, i.e. `key='h148_13'`. 

Due to its size, the produced `tracked_particles.hdf5` file is not stored on Github.
