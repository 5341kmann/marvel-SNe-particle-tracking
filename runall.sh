#!/bin/bash

# Summer 2025
# Grant Sackmann

# extracts PROGRAM_PATH from get_paths.py as describe in the config.py file
while IFS='=' read -r key value; do
  case "$key" in
    PROGRAM_PATH) PROGRAM_PATH="$value" 
    break;;
  esac
done < <(python get_paths.py) 

cd "$PROGRAM_PATH"ParticleTracking || { echo "Failed to change directory to PROGRAM_PATH: "$PROGRAM_PATH"ParticleTracking Please update config.py file with proper PROGRAM_PATH"; exit 1; }
echo "Changed directory to $PROGRAM_PATH"

date
python particletracking.py cptmarvel 1 &
python particletracking.py cptmarvel 2 &
python particletracking.py cptmarvel 3 &
python particletracking.py cptmarvel 5 &
wait
python particletracking.py cptmarvel 6 &
python particletracking.py cptmarvel 7 &
python particletracking.py cptmarvel 8 &
python particletracking.py cptmarvel 10 &
python particletracking.py cptmarvel 11 &
wait
python particletracking.py cptmarvel 13 &
python particletracking.py cptmarvel 16 &
# wait 

# python write_discharged.py
# include `wait` in between commands to do them in batches
