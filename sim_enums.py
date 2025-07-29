import os
from enum import Enum
from base import SIM_FOLDER_PATH
from utils import set_logger

# Class to handle Marvel and Justice League simulations

class MarvelSim(Enum):
    CAPTAINMARVEL = "cptmarvel"
    ELEKTRA = "elektra"
    ROGUE = "rogue"
    STORM = "storm"

    def get_path(self, snapnum) -> str:
        """
        Returns the path to the simulation data for the Marvel simulations.
        """
        
        logger = set_logger(self, "Get Simulation Path")

        if self.value in [ member.value for member in [MarvelSim.CAPTAINMARVEL, MarvelSim.ELEKTRA]]:
            return os.path.join(
                SIM_FOLDER_PATH,
                f"{self.value}.cosmo25cmb",
                f"{self.value}.cosmo25cmb.4096g5HbwK1BH",
                f"snapshots_200crit_{self.value}",
                f"{self.value}.cosmo25cmb.4096g5HbwK1BH.00{snapnum}",
            )
        elif self.value ==[member.value for member in [MarvelSim.ROGUE, MarvelSim.STORM]]:
            return os.path.join(
                SIM_FOLDER_PATH,
                f"{self.value}.cosmo25cmb",
                f"{self.value}.cosmo25cmb.4096g5HbwK1BH",
                f"{self.value}.cosmo25cmb.4096g5HbwK1BH.00{snapnum}",
                f"{self.value}.cosmo25cmb.4096g5HbwK1BH.00{snapnum}"
            )
        else:
            logger.debug(f"Unknown simulation: {self.name}. Cannot get path.")
            return ""

    def get_traceback_path(self) -> str:
        """
        Returns the path to the traceback data for the Marvel simulations.
        """
        logger = set_logger(self, "Get Traceback Path")

        if self.value in [member.value for member in MarvelSim if member != MarvelSim.CAPTAINMARVEL]:
            return os.path.join(
                SIM_FOLDER_PATH,
                f"{self.value}.cosmo25cmb",
                f"{self.value}.cosmo25cmb.4096g5HbwK1BH",
                f"{self.value}.cosmo25cmb.4096g5HbwK1BH.004096.M200.trace_back.hdf5",
            )
        elif self.value == MarvelSim.CAPTAINMARVEL.value:   
            return os.path.join(
                'Data',
                'cptmarvel.trace_back.hdf5'
            )
        else:
            logger.debug(f"Unknown simulation: {self.name}. Cannot get traceback path.")
            return ""
        

class JusticeLeagueSim(Enum):
    SANDRA = "h148"
    RUTH = "h229"
    SONIA = "h242"
    ELENA = "h329"

    def get_path(self, snapnum) -> str:

        logger = set_logger(self, "Get Simulation Path")
        """
        Returns the path to the simulation data for the Justice League simulations.
        """
        if self.value in [member.value for member in JusticeLeagueSim if member != JusticeLeagueSim.SANDRA]:
            return os.path.join(
                SIM_FOLDER_PATH,
                f"{self.value}.cosmo50PLK.3072g",
                f"{self.value}.cosmo50PLK.3072gst5HbwK1BH",
                f"snapshots_200crit_{self.value}",
                f"{self.value}.cosmo50PLK.3072gst5HbwK1BH.00{snapnum}",
            )
        elif self.value == JusticeLeagueSim.SANDRA.value:
            return os.path.join(
                SIM_FOLDER_PATH,
                f"{self.value}.cosmo50PLK.3072g",
                f"{self.value}.cosmo50PLK.3072g3HbwK1BH",
                f"snapshots_200crit_{self.value}",
                f"{self.value}.cosmo50PLK.3072g3HbwK1BH.00{snapnum}",
            )
        else:
            logger.debug(f"Unknown simulation: {self.name}. Cannot get path.")
            return ""

    def get_traceback_path(self) -> str:
        """
        Returns the path to the traceback data for the Justice League simulations.
        """

        logger = set_logger(self, "Get Traceback Path")

        if self.value in [member.value for member in JusticeLeagueSim if member.value != JusticeLeagueSim.SANDRA.value]:
            return os.path.join(
                SIM_FOLDER_PATH,
                f"{self.value}.cosmo50PLK.3072g",
                f"{self.value}.cosmo50PLK.3072gst5HbwK1BH",
                f"{self.value}.cosmo50PLK.3072gst5HbwK1BH.004096.M200.trace_back.hdf5",
            )
        elif self.value == JusticeLeagueSim.SANDRA.value:
            return os.path.join(
                SIM_FOLDER_PATH,
                f"{self.value}.cosmo50PLK.3072g",
                f"{self.value}.cosmo50PLK.3072g3HbwK1BH",
                f"{self.value}.cosmo50PLK.3072g3HbwK1BH.004096.M200.trace_back.hdf5",
            )
        else:
            logger.debug(f"Unknown simulation: {self.name}. Cannot get traceback path.")
            return ""