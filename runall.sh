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

# cd "$PROGRAM_PATH"ParticleTracking || { echo "Failed to change directory to PROGRAM_PATH: "$PROGRAM_PATH"ParticleTracking Please update config.py file with proper PROGRAM_PATH"; exit 1; }
# echo "Changed directory to $PROGRAM_PATH"

# date
# python particletracking.py cptmarvel 1 &
# python particletracking.py cptmarvel 2 &
# python particletracking.py cptmarvel 3 &
# python particletracking.py cptmarvel 5 &
# python particletracking.py cptmarvel 6 &
# wait
# python particletracking.py cptmarvel 7 &
# python particletracking.py cptmarvel 8 &
# python particletracking.py cptmarvel 10 &
# wait
# python particletracking.py cptmarvel 11 &
# python particletracking.py cptmarvel 13 &
# python particletracking.py cptmarvel 16 &
# wait 

cd "$PROGRAM_PATH"ParticlePropertiesCalc || { echo "Failed to change directory to PROGRAM_PATH: "$PROGRAM_PATH"ParticlePropertiesCalc Please update config.py file with proper PROGRAM_PATH"; exit 1; }
echo "Changed directory to $PROGRAM_PATH"

date
python propertiescalculating.py cptmarvel 1 &
python propertiescalculating.py cptmarvel 2 &
python propertiescalculating.py cptmarvel 3 &
python propertiescalculating.py cptmarvel 5 &
wait
python propertiescalculating.py cptmarvel 6 &
python propertiescalculating.py cptmarvel 7 &
python propertiescalculating.py cptmarvel 8 &
wait
python propertiescalculating.py cptmarvel 10 &
python propertiescalculating.py cptmarvel 11 &
python propertiescalculating.py cptmarvel 13 
python propertiescalculating.py cptmarvel 16 

# python write_discharged.py
# include `wait` in between commands to do them in batches
