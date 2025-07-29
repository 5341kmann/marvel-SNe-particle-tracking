import os
import logging

def set_logger(sim, description, halo_num=None):
    """
    Set up a logger for the analysis module.

    Parameters:
    sim (Enum): The simulation enum (MarvelSim or JusticeLeagueSim).
    halo_num (int): The halo number.
    description (str): A description for the logger, typically function name.

    Returns:
    logger: A configured logger instance.
    """
    # 1. Get the logger instance by a unique name (or a consistent application name)
    # Using the module name is common, but for specific halo logs, a unique name is good.
    logger_name = f"{description}-{sim.value}_{halo_num}" if halo_num is not None else f"{description}-{sim.value}"
    logger_instance = logging.getLogger(logger_name)

    # 2. Set the overall logging level for this logger
    logger_instance.setLevel(logging.DEBUG) # Catch all messages (DEBUG, INFO, WARNING, ERROR, CRITICAL)

    # Prevent adding handlers multiple times on subsequent calls 
    if not logger_instance.handlers: 
        # Create a file handler
        log_file_path = os.path.join('logs', f"{sim.value}_{halo_num}.log") if halo_num is not None else os.path.join('logs', f"{sim.value}.log")
        file_handler = logging.FileHandler(log_file_path)
        file_handler.setLevel(logging.DEBUG) # All messages to file

        # Create a formatter and set it for the handler
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        file_handler.setFormatter(formatter)
        
        # Add the file handler to the logger
        logger_instance.addHandler(file_handler)

        # Add a StreamHandler for notebook output
        stream_handler = logging.StreamHandler()
        stream_handler.setLevel(logging.INFO) # Only INFO and above to console to avoid clutter
        stream_handler.setFormatter(formatter)
        logger_instance.addHandler(stream_handler)

    return logger_instance