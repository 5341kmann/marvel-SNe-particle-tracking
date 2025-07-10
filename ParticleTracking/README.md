# ParticleTracking

## Author
Grant Sackmann

## Description
ParticleTracking is a package designed for tracking particles in simulation data. This README provides instructions on how to use the key scripts within the package.
In progress — INCOMPLETE

## Installation
[Instructions on how to install the package, if applicable]

## Usage

### generate_haloid.py
The `generate_haloid.py` script is used to generate halo IDs for your data. Before running this script, ensure that you populate the `snapnums` variable with the appropriate snapshot numbers relevant to your simulation.

### particletracking.py

This script provides a command-line interface (CLI) for tracking particles in a given dataset. It utilizes various functions in base.py to track particle movement and generate relevant properties for each particle through timesteps.

#### Notes

To use the `particletracking.py` script, run the following command in your terminal:

python particletracking.py ['SIM'] [z0haloid]

for bulk calculations use runall.sh in project base dirtory

Make sure to update PROGRAM_PATH in config.py located in project base directory
if you are on quirm.math.grinnell.edu you should not have have to update SIM_FOLDER_PATH