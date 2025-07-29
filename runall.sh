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

# -- UNCOMMENT THIS SECTION TO TRACK PARTICLES --
cd "$PROGRAM_PATH"ParticleTracking || { echo "Failed to change directory to PROGRAM_PATH: "$PROGRAM_PATH"ParticleTracking Please update config.py file with proper PROGRAM_PATH"; exit 1; }
echo "Changed directory to $PROGRAM_PATH"
FILE_NAME="particletracking.py"
# -----------------------------------------------

# # --- UNCOMMENT THIS SECTION TO CALCULATE PROPERTIES OF PARTICLES ---
# cd "$PROGRAM_PATH"ParticlePropertiesCalc || { echo "Failed to change directory to PROGRAM_PATH: "$PROGRAM_PATH"ParticlePropertiesCalc Please update config.py file with proper PROGRAM_PATH"; exit 1; }
# echo "Changed directory to $PROGRAM_PATH"
# FILE_NAME="propertiescalculating.py"
# # -------------------------------------------------------------------
# date


# -----------------------MARVEL SIMULATIONS-----------------------
# CAPTAIN MARVEL
# HALO_NUMS=(1 2 3 5 6 7 8 10 11 13 16)
# SIMULATION_NAME="cptmarvel"
# for HALO_NUM in "${HALO_NUMS[@]}"; do
#   echo "Processing simulation: $SIMULATION_NAME, Z=0 halo: $HALO_NUM"
#   python propertiescalculating.py "$SIMULATION_NAME" "$HALO_NUM"
# done

# ELEKTRA
# HALO_NUMS=(1 2 3 4 5 6 9 10 12 15 67)
# SIMULATION_NAME="elektra"
# for HALO_NUM in "${HALO_NUMS[@]}"; do
#   echo "Processing simulation: $SIMULATION_NAME, Z=0 halo: $HALO_NUM"
#   python "$FILE_NAME" "$SIMULATION_NAME" "$HALO_NUM"
# done

# ROGUE
HALO_NUMS=(2 3 5 7 8 9 10 11 12 14 15 16 17 22 26 28 31 34 35 37 2659 3319 5601 6982)
SIMULATION_NAME="rogue"
for HALO_NUM in "${HALO_NUMS[@]}"; do
  echo "Processing simulation: $SIMULATION_NAME, Z=0 halo: $HALO_NUM"
  python "$FILE_NAME" "$SIMULATION_NAME" "$HALO_NUM"
done

# STORM
HALO_NUMS=(1 2 3 4 5 6 7 8 10 11 12 14 15 17 19 23 26 31 43 44 81 103 119 191 906 2163 2815 3385 3494 4908 4981)
for HALO_NUM in "${HALO_NUMS[@]}"; do
  echo "Processing simulation: $SIMULATION_NAME, Z=0 halo: $HALO_NUM"
  python "$FILE_NAME" "$SIMULATION_NAME" "$HALO_NUM"
done
# ------------------------------------------------------------------


# -----------------DC JUSTICE LEAGUE SIMULATIONS--------------------
# # SANDRA
# HALO_NUMS=(4 6 7 10 12 23 27 34 38 55 65 249 251 282)
# SIMULATION_NAME="h148"
# for HALO_NUM in "${HALO_NUMS[@]}"; do
#   echo "Processing simulation: $SIMULATION_NAME, Z=0 halo: $HALO_NUM"
#   python "$FILE_NAME" "$SIMULATION_NAME" "$HALO_NUM"
# done

# # RUTH
# HALO_NUMS=(14 18 20 22 49)
# SIMULATION_NAME="h229"
# for HALO_NUM in "${HALO_NUMS[@]}"; do
#   echo "Processing simulation: $SIMULATION_NAME, Z=0 halo: $HALO_NUM"
#   python "$FILE_NAME" "$SIMULATION_NAME" "$HALO_NUM"
# done

# # SONIA
# HALO_NUMS=(8 10 21 30 38 69 401)
# SIMULATION_NAME="h242"
# for HALO_NUM in "${HALO_NUMS[@]}"; do
#   echo "Processing simulation: $SIMULATION_NAME, Z=0 halo: $HALO_NUM"
#   python "$FILE_NAME" "$SIMULATION_NAME" "$HALO_NUM"
# done

# # ELENA
# HALO_NUMS=(7 29 117)
# SIMULATION_NAME="h329"
# for HALO_NUM in "${HALO_NUMS[@]}"; do
#   echo "Processing simulation: $SIMULATION_NAME, Z=0 halo: $HALO_NUM"
#   python "$FILE_NAME" "$SIMULATION_NAME" "$HALO_NUM"
# done
# ------------------------------------------------------------------