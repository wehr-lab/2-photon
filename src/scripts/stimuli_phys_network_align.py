import os 
import sys
from pathlib import Path

import numpy as np
import networkx as nx
import h5py 
import glob 

## import custom modules 
currentDir = Path(__file__).resolve().parent.parent
configDir = currentDir.parent / "config"
sys.path.append(str(configDir))

## define path variables 
from settings import DATA_PATH, CODE_PATH
sys.path.append(str(CODE_PATH))

import ephys.behavior as behavior
import ephys.cell as cell
import ephys.cell_ensemble as cell_ensemble
import ephys.events as events
import ephys.events_process as events_process
import ephys.session as session

## load .mat files and typecast to python classes 
DATA_PATH = Path(DATA_PATH)
print("DATA_PATH:", DATA_PATH)



print("run complete!")