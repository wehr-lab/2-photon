"""
This is a configuration file for network computation project. 
It contains settings for the project, including paths to data directories, model configurations, and other parameters.
"""

import os 
import sys

## Project data directories 
DATA_PATH = {}
DATA_PATH["universal"] = "/".join(os.path.dirname(__file__).split("/")[:-2]) + "/data/"

# print(DATA_PATH)
DATA_PATH["toneDecode"] = os.path.join(DATA_PATH["universal"], "tone_decode_data")


## Project code directories
CODE_PATH = "/".join(os.path.dirname(__file__).split("/")[:-1]) + "/src/"
print(CODE_PATH)