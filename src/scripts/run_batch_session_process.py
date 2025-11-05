#!/usr/bin/env python3
"""
Batch processing script for arousal network analysis across multiple sessions.
"""

import sys
import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
import pickle
import warnings
import subprocess

warnings.filterwarnings("ignore")

# Setup paths
repoDir = Path(__file__).resolve().parent.parent.parent
# print(repoDir)
configDir = repoDir / "config"
sys.path.append(str(configDir))

from settings import DATA_PATH, CODE_PATH

sys.path.append(str(CODE_PATH["2-photon"]))
